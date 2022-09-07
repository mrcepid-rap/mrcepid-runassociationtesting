import pandas.core.series

from runassociationtesting.association_pack import AssociationPack
from runassociationtesting.association_resources import *
from runassociationtesting.linear_model.proccess_model_output import merge_glm_staar_runs
from runassociationtesting.thread_utility import ThreadUtility


class ExtractVariants:

    def __init__(self, association_pack: AssociationPack):

        self._association_pack = association_pack

        # 1. Define our gene-list and make 'gene_info' objects of them (which are Pandas series classes)
        gene_infos = []
        chromosomes = set()
        # If we are doing extraction based on individual SNPs or a Gene list, we need to make a 'fake' gene info
        # but find all chromosomes those SNPS/Genes lie on
        if self._association_pack.is_non_standard_tar:
            gene_info, returned_chromosomes = process_snp_or_gene_tar(self._association_pack.is_snp_tar,
                                                                      self._association_pack.is_gene_tar,
                                                                      self._association_pack.tarball_prefixes[0])
            gene_infos.append(gene_info)
            chromosomes = returned_chromosomes
        else:
            for gene in self._association_pack.gene_ids:
                transcripts_table = build_transcript_table()
                gene_info = get_gene_id(gene, transcripts_table)
                gene_infos.append(gene_info)
                chromosomes.add(gene_info['chrom'])

        # 2. Download variant VEP annotations
        print("Loading VEP annotations...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='An extraction thread failed',
                                       incrementor=5,
                                       thread_factor=4)
        for chromosome in chromosomes:
            thread_utility.launch_job(self._download_vep,
                                      chromosome = chromosome)
        thread_utility.collect_futures()

        # 3. Filter relevant files to individuals we want to keep
        print("Filtering variant files to appropriate individuals...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='An extraction thread failed',
                                       incrementor=20,
                                       thread_factor=4)

        # if, elif, else simply depends on which type of tarball we are using
        if self._association_pack.is_snp_tar:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(self._filter_individuals,
                                          tarball_prefix=tarball_prefix,
                                          chromosome='SNP')
        elif self._association_pack.is_gene_tar:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(self._filter_individuals,
                                          tarball_prefix=tarball_prefix,
                                          chromosome='GENE')
        else:
            for chromosome in chromosomes:
                for tarball_prefix in self._association_pack.tarball_prefixes:
                    thread_utility.launch_job(self._filter_individuals,
                                              tarball_prefix=tarball_prefix,
                                              chromosome=chromosome)
        thread_utility.collect_futures()

        # 4. Actually collect variant information per-gene
        print("Extracting variant information...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='An extraction thread failed',
                                       incrementor=20,
                                       thread_factor=2)
        for gene_info in gene_infos:

            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(self._annotate_variants,
                                          tarball_prefix=tarball_prefix,
                                          gene_info=gene_info,
                                          chromosomes=chromosomes if self._association_pack.is_non_standard_tar else None)
        self.outputs = []
        future_results = thread_utility.collect_futures()
        for result in future_results:
            self.outputs.extend(result)

        # 5. And run a linear and STAAR model(s) for all genes
        print("Running linear models...")
        GLMRunner(self._association_pack, gene_infos=gene_infos)
        STAARRunner(self._association_pack, gene_infos=gene_infos)

        self.outputs.extend(merge_glm_staar_runs(self._association_pack)) # Function merges STAAR and GLM results together

        # 6. Finally add the phenotypes/covariates table to the outputs
        os.rename('phenotypes_covariates.formatted.txt', self._association_pack.output_prefix + '.phenotypes_covariates.formatted.tsv')
        self.outputs.append(self._association_pack.output_prefix + '.phenotypes_covariates.formatted.tsv')

    def _download_vep(self, chromosome: str) -> None:

        vep_dx = dxpy.DXFile(self._association_pack.bgen_dict[chromosome]['vep'])
        dxpy.download_dxfile(vep_dx.get_id(), chromosome + ".filtered.vep.tsv.gz")

    @staticmethod
    def _filter_individuals(tarball_prefix: str, chromosome: str) -> None:
        # And filter the relevant SAIGE file to just the individuals we want so we can get actual MAC
        cmd = "bcftools view --threads 4 -S /test/SAMPLES_Include.txt -Ob -o /test/" + tarball_prefix + "." + chromosome + ".saige_input.bcf /test/" + tarball_prefix + "." + chromosome + ".SAIGE.bcf"
        run_cmd(cmd, True)

    def _annotate_variants(self, tarball_prefix: str, gene_info: pandas.core.series.Series, chromosomes: set) -> list:

        # This is a bit confusing, so explaining in full.
        # We need to annotate EACH GENE separately EXCEPT when running a SNP/GENE list tarball, SO...
        # 1. If we have a SNP/GENE tar, we may need to load in variant annotations for multiple chromosomes, so we set 'chromosomes'
        # to a set of chromosomes that we extract from the SNP/GENE tar.
        # 2. If just a single gene or gene list (chromosomes = None), only need to load the data for the chromosome that specific Gene is on
        if chromosomes is None:
            chromosome = gene_info['chrom']
            variant_index = pd.read_csv(gzip.open(chromosome + ".filtered.vep.tsv.gz", 'rt'),
                                        sep = "\t",
                                        dtype={'SIFT': str, 'POLYPHEN': str})
        else:
            variant_index = []
            for chromosome in chromosomes:
                variant_index.append(pd.read_csv(gzip.open(chromosome + ".filtered.vep.tsv.gz", 'rt'),
                                                 sep="\t",
                                                 dtype={'SIFT': str,'POLYPHEN': str}))
            variant_index = pd.concat(variant_index)

        # Need to get the variants from the SAIGE groupfile:
        with open(tarball_prefix + "." + gene_info['chrom'] + ".SAIGE.groupFile.txt") as saige_group_file:
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
        cmd = "bcftools view --threads 2 -i \'ID=@/test/" + tarball_prefix + "." + gene_info['SYMBOL'] + ".variants.txt\' -Ob -o /test/" + tarball_prefix + "." + gene_info['SYMBOL'] + ".variant_filtered.bcf /test/" + tarball_prefix + "." + gene_info['chrom'] + ".saige_input.bcf"
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

        variant_file = self._association_pack.output_prefix + "." + tarball_prefix + "." + gene_info['SYMBOL'] + '.variant_table.tsv'
        carriers_file = self._association_pack.output_prefix + "." + tarball_prefix + "." + gene_info['SYMBOL'] + '.carriers_formatted.tsv'
        geno_table.to_csv(path_or_buf=variant_file, index = False, sep="\t", na_rep='NA')
        carriers_table.to_csv(path_or_buf=carriers_file, index = False, sep="\t", na_rep='NA')

        return [variant_file, carriers_file]
