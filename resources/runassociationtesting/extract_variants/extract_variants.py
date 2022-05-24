import csv
import warnings

import pandas as pd
import pandas.core.series

from ..tool_runners.glm_runner import GLMRunner
from ..association_pack import AssociationPack
from ..association_resources import *
from ..thread_utility import ThreadUtility

class ExtractVariants:

    def __init__(self, association_pack: AssociationPack):

        self._association_pack = association_pack
        transcripts_table = build_transcript_table()

        # 1. Define our gene-list and make 'gene_info' objects of them (which are Pandas series classes)
        gene_infos = []
        chromosomes = set()
        # If we are doing extraction based on individual SNPs, we need to make a 'fake' gene info but find all chromosomes
        if association_pack.is_snp_tar:
            print("Running extract variants in SNP mode...")
            gene_info = pd.Series({'chrom': 'SNP', 'SYMBOL': 'SNP'})
            gene_info.name = 'ENST00000000000'
            gene_infos.append(gene_info)
            sparse_matrix = csv.DictReader(open(association_pack.tarball_prefixes[0] + '.SNP.variants_table.STAAR.tsv', 'r'),
                                           delimiter="\t",
                                           quoting=csv.QUOTE_NONE)
            for row in sparse_matrix:
                chromosomes.add(str(row['chrom']))

            # And filter the relevant SAIGE file to just the individuals we want so we can get actual MAC
            cmd = "bcftools view --threads 4 -S /test/SAMPLES_Include.txt -Ob -o /test/" + association_pack.tarball_prefixes[0] + ".SNP.saige_input.bcf /test/" + association_pack.tarball_prefixes[0] + ".SNP.SAIGE.bcf"
            run_cmd(cmd, True)
        else:
            for gene in self._association_pack.gene_ids:
                gene_info = get_gene_id(gene, transcripts_table)
                gene_infos.append(gene_info)
                chromosomes.add(gene_info['chrom'])

        # 2. Load per-chromosome genotype/variant information
        print("Loading VEP annotations...")
        thread_utility = ThreadUtility(self._association_pack.threads,error_message='An extraction thread failed',incrementor=5,thread_factor=4)
        for chromosome in chromosomes:
            thread_utility.launch_job(self._download_vep,
                                      chromosome = chromosome,
                                      vep_dx=dxpy.DXFile(association_pack.bgen_dict[chromosome]['vep']))
        thread_utility.collect_futures()

        # 3. Filter relevant files to individuals we want to keep
        print("Filtering variant files to appropriate individuals...")
        thread_utility = ThreadUtility(self._association_pack.threads,error_message='An extraction thread failed',incrementor=20,thread_factor=4)
        for chromosome in set(['SNP']) if self._association_pack.is_snp_tar else chromosomes: # Just allows me to filter SNP vcfs if required...
            for tarball_prefix in association_pack.tarball_prefixes:
                thread_utility.launch_job(self._filter_individuals,
                                          tarball_prefix = tarball_prefix,
                                          chromosome = chromosome)
        thread_utility.collect_futures()

        # 4. Actually collect variant information per-gene
        print("Extracting variant information...")
        genes_to_run = [] # This just enables easy parallelisation with GLMRunner() – we use this to run all genes at once rather than 1 by 1
        thread_utility = ThreadUtility(self._association_pack.threads,error_message='An extraction thread failed',incrementor=20,thread_factor=2)
        for gene_info in gene_infos:
            genes_to_run.append(gene_info.name)
            for tarball_prefix in association_pack.tarball_prefixes:
                thread_utility.launch_job(self._annotate_variants,
                                          tarball_prefix=tarball_prefix,
                                          gene_info=gene_info,
                                          output_prefix=association_pack.output_prefix,
                                          chromosomes=chromosomes if self._association_pack.is_snp_tar else None)
        self.outputs = []
        future_results = thread_utility.collect_futures()
        for result in future_results:
            self.outputs.extend(result)

        # 5. And run a linear model for all genes
        print("Running linear models...")
        glm_run = GLMRunner(association_pack, genes_to_run=genes_to_run)
        self.outputs.extend(glm_run.outputs)
        os.rename('phenotypes_covariates.formatted.txt', association_pack.output_prefix + '.phenotypes_covariates.formatted.tsv')
        self.outputs.append(association_pack.output_prefix + '.phenotypes_covariates.formatted.tsv')

    @staticmethod
    def _download_vep(chromosome: str, vep_dx: dxpy.DXFile) -> None:

        dxpy.download_dxfile(vep_dx.get_id(), chromosome + ".filtered.vep.tsv.gz")

    @staticmethod
    def _filter_individuals(tarball_prefix: str, chromosome: str) -> None:
        # And filter the relevant SAIGE file to just the individuals we want so we can get actual MAC
        cmd = "bcftools view --threads 4 -S /test/SAMPLES_Include.txt -Ob -o /test/" + tarball_prefix + "." + chromosome + ".saige_input.bcf /test/" + tarball_prefix + "." + chromosome + ".SAIGE.bcf"
        run_cmd(cmd, True)

    @staticmethod
    def _annotate_variants(tarball_prefix: str, gene_info: pandas.core.series.Series, output_prefix: str, chromosomes: set) -> list:

        # The point here is that if we have a SNP tar, we may need to load in variant annotations for multiple chromosomes
        # If just a single gene (chromosomes = None), only need to load the data for the chromosome that gene is on
        if chromosomes is None:
            chromosome = gene_info['chrom']
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                variant_index = pd.read_csv(gzip.open(chromosome + ".filtered.vep.tsv.gz", 'rt'), sep = "\t")
        else:
            variant_index = []
            for chromosome in chromosomes:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    variant_index = pd.read_csv(gzip.open(chromosome + ".filtered.vep.tsv.gz", 'rt'), sep="\t")
            variant_index = pd.concat(variant_index)

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

        variant_file = output_prefix + "." + tarball_prefix + "." + gene_info['SYMBOL'] + '.variant_table.tsv'
        carriers_file = output_prefix + "." + tarball_prefix + "." + gene_info['SYMBOL'] + '.carriers_formatted.tsv'
        geno_table.to_csv(path_or_buf=variant_file, index = False, sep="\t", na_rep='NA')
        carriers_table.to_csv(path_or_buf=carriers_file, index = False, sep="\t", na_rep='NA')

        return [variant_file, carriers_file]