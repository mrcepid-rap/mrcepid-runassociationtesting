import os
import pandas as pd
import gzip
import pandas.core.series
import dxpy

from ..tool_runners.glm_runner import GLMRunner
from ..association_pack import AssociationPack
from ..association_resources import *
from ..thread_utility import ThreadUtility

class ExtractVariants:

    def __init__(self, association_pack: AssociationPack, gene_ids: str):

        self._association_pack = association_pack

        # 1. Define our gene-list and make 'gene_info' objects of them (which are Pandas series classes)
        self._transcripts_table = pd.read_csv(gzip.open('transcripts.tsv.gz', 'rt'), sep = "\t")
        self._transcripts_table = self._transcripts_table.rename(columns={'#chrom':'chrom'})
        self._transcripts_table = self._transcripts_table.set_index('ENST')
        self._gene_ids = gene_ids.split(",")
        gene_infos = []
        for gene in self._gene_ids:
            gene_infos.append(self._get_gene_id(gene))

        # 2. Load per-chromosome genotype/variant information and filter to individuals we want to keep
        thread_utility = ThreadUtility(self._association_pack.threads,error_message='An extraction thread failed',incrementor=1,thread_factor=4)
        chromosomes = set()
        for gene_info in gene_infos:
            if gene_info['chrom'] not in chromosomes:
                for tarball_prefix in association_pack.tarball_prefixes:
                    thread_utility.launch_job(self._process_chromosome,
                                              tarball_prefix = tarball_prefix,
                                              chromosome = gene_info['chrom'],
                                              vep_dx=dxpy.DXFile(association_pack.bgen_dict[gene_info['chrom']]['vep']))
            chromosomes.add(gene_info['chrom'])
        thread_utility.collect_futures()

        # 3. Actually collect the information per-gene
        genes_to_run = []
        thread_utility = ThreadUtility(self._association_pack.threads,error_message='An extraction thread failed',incrementor=1,thread_factor=2)
        for gene_info in gene_infos:
            genes_to_run.append(gene_info.name)
            for tarball_prefix in association_pack.tarball_prefixes:
                thread_utility.launch_job(self._annotate_variants,
                                          tarball_prefix=tarball_prefix,
                                          gene_info=gene_info)
        self.outputs = []
        future_results = thread_utility.collect_futures()
        for result in future_results:
            self.outputs.extend(result)

        # 4. And run a linear model for all genes
        glm_run = GLMRunner(association_pack, genes_to_run=genes_to_run)
        self.outputs.extend(glm_run.outputs)
        os.rename('phenotypes_covariates.formatted.txt', association_pack.output_prefix + '.phenotypes_covariates.formatted.tsv')
        self.outputs.append(association_pack.output_prefix + '.phenotypes_covariates.formatted.tsv')

    def _get_gene_id(self, gene_id: str) -> pandas.core.series.Series:

        if 'ENST' in gene_id:
            print("gene_id – " + gene_id + " – looks like an ENST value... validating...")
            try:
                gene_info = self._transcripts_table.loc[gene_id]
                print("Found one matching ENST (%s - %s)... proceeding..." % (gene_id, gene_info['coord']))
            except KeyError:
                print("Did not find a transcript with ENST value %s... terminating..." % gene_id)
        else:
            print("gene_id – " + gene_id + " – does not look like an ENST value, searching for symbol instead...")
            found_rows = self._transcripts_table[self._transcripts_table['SYMBOL'] == gene_id]
            if len(found_rows) == 1:
                found_enst = found_rows.index[0]
                gene_info = self._transcripts_table.loc[found_enst]
                print("Found one matching ENST (%s - %s) for SYMBOL %s... proceeding..." % (found_enst, gene_info['coord'], gene_id))
            elif len(found_rows) > 1:
                gene_info = None
                print("Found %i ENST IDs (%s) for SYMBOL %s... Please re-run using exact ENST to ensure consistent results..." % (len(found_rows), ','.join(found_rows.index.to_list()), gene_id))
                raise Exception("Multiple ENST IDs")
            else:
                gene_info = None
                print("Did not find an associated ENST ID for SYMBOL %s... Please re-run after checking SYMBOL/ENST used..." % (gene_id))

        return gene_info

    def _process_chromosome(self, tarball_prefix: str, chromosome: str, vep_dx: dxpy.DXFile) -> None:

        dxpy.download_dxfile(vep_dx.get_id(), chromosome + ".filtered.vep.tsv.gz")

        # And filter the relevant SAIGE file to just the individuals we want so we can get actual MAC
        cmd = "bcftools view --threads 4 -S /test/SAMPLES_Include.txt -Ob -o /test/" + tarball_prefix + "." + chromosome + ".saige_input.bcf /test/" + tarball_prefix + "." + chromosome + ".SAIGE.bcf"
        run_cmd(cmd, True)

    def _annotate_variants(self, tarball_prefix: str, gene_info: pandas.core.series.Series) -> list:

        chromosome = gene_info['chrom']
        variant_index = pd.read_csv(gzip.open(chromosome + ".filtered.vep.tsv.gz", 'rt'), sep = "\t")

        # Need to get the variants from the SAIGE groupfile:
        with open(tarball_prefix + "." + chromosome + ".SAIGE.groupFile.txt") as saige_group_file:
            var_ids = []
            var_file = open(tarball_prefix + "." + gene_info['SYMBOL'] + '.variants.txt', 'w')
            for line in saige_group_file:
                data = line.rstrip().split("\t")
                if data[0] == gene_info.name:
                    for i in range(1,len(data)):
                        currID = data[i].replace('_',':').replace('/',':')
                        var_file.write(currID + "\n")
                        var_ids.append(currID)
                    break
            var_file.close()

        relevant_vars = variant_index[variant_index['varID'].isin(var_ids)]

        # Filter to the variants for this gene
        cmd = "bcftools view --threads 2 -i \'ID=@/test/" + tarball_prefix + "." + gene_info['SYMBOL'] + ".variants.txt\' -Ob -o /test/" + tarball_prefix + "." + gene_info['SYMBOL'] + ".variant_filtered.bcf /test/" + tarball_prefix + "." + chromosome + ".saige_input.bcf"
        run_cmd(cmd, True)
        cmd = "bcftools +fill-tags --threads 4 -Ob -o /test/" + tarball_prefix + "." + gene_info['SYMBOL'] + ".final.bcf /test/" + tarball_prefix + "." + gene_info['SYMBOL'] + ".variant_filtered.bcf"
        run_cmd(cmd, True)

        # Now get actual annotations back in:
        cmd = "bcftools query -f \'%ID\\t%MAF\\t%AC\\t%AC_Het\\t%AC_Hom\\n\' -o /test/" + tarball_prefix + "." + gene_info['SYMBOL'] + ".annotated_vars.txt /test/" + tarball_prefix + "." + gene_info['SYMBOL'] + ".final.bcf"
        run_cmd(cmd, True)
        # And get a list of individuals with a variant:
        cmd = "bcftools query -i \"GT=\'alt\'\" -f \'[%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\n]\' -o /test/" + tarball_prefix + "." + gene_info['SYMBOL'] + ".carriers.txt /test/" + tarball_prefix + "." + gene_info['SYMBOL'] + ".final.bcf"
        run_cmd(cmd, True)

        geno_table = pd.read_csv(tarball_prefix + "." + gene_info['SYMBOL'] + ".annotated_vars.txt",
                                 sep = "\t",
                                 names=['varID','MAF_tested','AC_tested','AC_tested_Het','AC_tested_Hom'])
        geno_table = pd.merge(relevant_vars, geno_table, on='varID', how="left")

        carriers_table = pd.read_csv(tarball_prefix + "." + gene_info['SYMBOL'] + ".carriers.txt",
                                     sep="\t",
                                     names=['CHROM','POS','varID','REF','ALT','IID','GT'])

        variant_file = tarball_prefix + "." + gene_info['SYMBOL'] + '.variant_table.tsv'
        carriers_file = tarball_prefix + "." + gene_info['SYMBOL'] + '.carriers_formated.tsv'
        geno_table.to_csv(path_or_buf=variant_file, index = False, sep="\t", na_rep='NA')
        carriers_table.to_csv(path_or_buf=carriers_file, index = False, sep="\t", na_rep='NA')

        return [variant_file, carriers_file]