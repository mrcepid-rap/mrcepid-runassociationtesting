from os.path import exists
from runassociationtesting.association_resources import *
from runassociationtesting.burden.tool_runners.tool_runner import ToolRunner
from runassociationtesting.thread_utility import ThreadUtility


class SAIGERunner(ToolRunner):

    def run_tool(self) -> None:

        # 1. Run SAIGE step one without parallelisation
        print("Running SAIGE step 1...")
        self._saige_step_one()

        # 2. Run SAIGE step two WITH parallelisation by chromosome
        print("Running SAIGE step 2...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A SAIGE thread failed', 
                                       incrementor=10, 
                                       thread_factor=1)
        
        for chromosome in get_chromosomes():
            for tarball_prefix in self._association_pack.tarball_prefixes:
                if exists(tarball_prefix + "." + chromosome + ".SAIGE.bcf"):
                    self._prep_group_file(tarball_prefix, chromosome)
                    thread_utility.launch_job(class_type=self._saige_step_two,
                                              tarball_prefix=tarball_prefix,
                                              chromosome=chromosome)
        future_results = thread_utility.collect_futures()

        # 3. Gather preliminary results
        print("Gathering SAIGE mask-based results...")
        completed_gene_tables = []
        log_file = open(self._output_prefix + '.SAIGE_step2.log', 'w')
        for result in future_results:
            tarball_prefix, finished_chromosome = result
            completed_gene_tables.append(self._process_saige_output(tarball_prefix, finished_chromosome))
            log_file.write(f'{tarball_prefix + "-" + finished_chromosome:{"-"}^{50}}')
            with open(tarball_prefix + "." + finished_chromosome + ".SAIGE_step2.log", 'r') as current_log:
                for line in current_log:
                    log_file.write(line)
                current_log.close()
        log_file.close()

        # 4. Run per-marker tests, if requested
        completed_marker_chromosomes = []
        if self._association_pack.run_marker_tests:
            print("Running per-marker tests...")
            thread_utility = ThreadUtility(self._association_pack.threads,
                                           error_message='A SAIGE marker thread failed',
                                           incrementor=1,
                                           thread_factor=4)

            for chromosome in get_chromosomes():
                thread_utility.launch_job(class_type=self._saige_marker_run,
                                          chromosome=chromosome)
                completed_marker_chromosomes.append(chromosome)

            future_results = thread_utility.collect_futures()
            markers_log_file = open(self._output_prefix + '.SAIGE_markers.log', 'w')
            for result in future_results:
                finished_chromosome = result
                markers_log_file.write(f'{finished_chromosome + ".log":{"-"}^{50}}')
                with open(finished_chromosome + ".SAIGE_markers.log", 'r') as current_log:
                    for line in current_log:
                        markers_log_file.write(line)
                    current_log.close()
            markers_log_file.close()

        # 5. Process final results
        print("Processing final SAIGE output...")
        self._outputs.extend(self._annotate_saige_output(completed_gene_tables, completed_marker_chromosomes))

    # Run rare variant association testing using SAIGE-GENE
    def _saige_step_one(self) -> None:

        # See the README.md for more information on these parameters
        # Just to note – I previously tried to implement the method that includes variance ratio estimation. However,
        # there are too few low MAC variants in the genotype files to perform this step accurately. The SAIGE
        # documentation includes this step, but I am very unsure how it works...
        cmd = 'step1_fitNULLGLMM.R ' \
                    '--phenoFile=/test/phenotypes_covariates.formatted.txt ' \
                    f'--phenoCol={self._association_pack.pheno_names[0]} ' \
                    '--isCovariateTransform=FALSE ' \
                    '--sampleIDColinphenoFile=IID ' \
                    f'--outputPrefix=/test/{self._association_pack.pheno_names[0]}.SAIGE_OUT ' \
                    '--sparseGRMFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx ' \
                    '--sparseGRMSampleIDFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt ' \
                    f'--nThreads={str(self._association_pack.threads)} ' \
                    '--LOCO=FALSE ' \
                    '--skipModelFitting=FALSE ' \
                    '--useSparseGRMtoFitNULL=TRUE ' \
                    '--skipVarianceRatioEstimation=TRUE '
        if self._association_pack.is_binary:
            cmd = cmd + '--traitType=binary '
        else:
            cmd = cmd + '--traitType=quantitative '

        # Manage covariates
        if self._association_pack.sex == 2:
            default_covars = ['PC' + str(x) for x in range(1, 11)] + ['age', 'age_squared', 'sex', 'wes_batch']
        else:
            default_covars = ['PC' + str(x) for x in range(1, 11)] + ['age', 'age_squared', 'wes_batch']
        all_covariates = [','.join(default_covars)]  # A list to manage appending additional covariates
        if len(self._association_pack.found_quantitative_covariates) > 0:
            all_covariates.append(','.join(self._association_pack.found_quantitative_covariates))
        if len(self._association_pack.found_categorical_covariates) > 0:
            cat_covars_join = ','.join(self._association_pack.found_categorical_covariates)
            all_covariates.append(cat_covars_join)
            cmd = cmd + '--qCovarColList=wes_batch,' + cat_covars_join + ' '
        else:
            cmd = cmd + '--qCovarColList=wes_batch '
        cmd = cmd + '--covarColList=' + ','.join(all_covariates)

        run_cmd(cmd, True, self._output_prefix + ".SAIGE_step1.log", print_cmd=True)

    # This exists for a very stupid reason – they _heavily_ modified the groupFile for v1.0 and I haven't gone back
    # to change how this file is made in 'collapse variants'
    @staticmethod
    def _prep_group_file(tarball_prefix: str, chromosome: str) -> None:

        with open(tarball_prefix + '.' + chromosome + '.SAIGE.groupFile.txt', 'r') as group_file:
            modified_group = open(tarball_prefix + '.' + chromosome + '.SAIGE_v1.0.groupFile.txt', 'w')
            for line in group_file:
                data = line.rstrip().split('\t')
                mod_data = [data[0], 'var']
                found = set()
                for var in data[1:]:
                    var = var.translate(str.maketrans('_/', '::'))
                    if var not in found:
                        mod_data.append(var)
                        found.add(var)
                modified_group.writelines(' '.join(mod_data) + "\n")
                mod_annote = [data[0], 'anno']
                for i in range(2, len(mod_data)):
                    mod_annote.append('foo')
                modified_group.write(' '.join(mod_annote) + "\n")
            group_file.close()
            modified_group.close()

    # This is a helper function to parallelise SAIGE step 2 by chromosome
    # This returns the tarball_prefix and chromosome number to make it easier to generate output
    def _saige_step_two(self, tarball_prefix: str, chromosome: str) -> tuple:

        cmd = f'bcftools view --threads 1 -S /test/SAMPLES_Include.txt -Ob ' \
              f'-o /test/{tarball_prefix}.{chromosome}.saige_input.bcf ' \
              f'/test/{tarball_prefix}.{chromosome}.SAIGE.bcf'
        run_cmd(cmd, True)
        cmd = f'bcftools index --threads 1 /test/{tarball_prefix}.{chromosome}.saige_input.bcf'
        run_cmd(cmd, True)
        
        # See the README.md for more information on these parameters
        cmd = 'step2_SPAtests.R ' \
              f'--vcfFile=/test/{tarball_prefix}.{chromosome}.saige_input.bcf ' \
              '--vcfField=GT ' \
              f'--GMMATmodelFile=/test/{self._association_pack.pheno_names[0]}.SAIGE_OUT.rda ' \
              '--sparseGRMFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx ' \
              '--sparseGRMSampleIDFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt ' \
              '--LOCO=FALSE ' \
              f'--SAIGEOutputFile=/test/{tarball_prefix}.{chromosome}.SAIGE_OUT.SAIGE.gene.txt ' \
              f'--groupFile=/test/{tarball_prefix}.{chromosome}.SAIGE_v1.0.groupFile.txt ' \
              '--is_output_moreDetails=TRUE ' \
              '--maxMAF_in_groupTest=0.5 ' \
              '--maxMissing=1 ' \
              f'--chrom={chromosome} ' \
              '--annotation_in_groupTest=foo '

        if self._association_pack.is_binary:
            cmd = cmd + '--is_Firth_beta=TRUE'

        run_cmd(cmd, True, tarball_prefix + "." + chromosome + ".SAIGE_step2.log")

        return tarball_prefix, chromosome

    def _saige_marker_run(self, chromosome: str) -> str:

        process_bgen_file(self._association_pack.bgen_dict[chromosome], chromosome)

        cmd = 'step2_SPAtests.R ' \
              f'--bgenFile=/test/{chromosome}.markers.bgen ' \
              f'--bgenFileIndex=/test/{chromosome}.markers.bgen.bgi ' \
              f'--sampleFile=/test/{chromosome}.markers.bolt.sample ' \
              f'--GMMATmodelFile=/test/{self._association_pack.pheno_names[0]}.SAIGE_OUT.rda ' \
              '--sparseGRMFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx ' \
              '--sparseGRMSampleIDFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt ' \
              f'--SAIGEOutputFile=/test/{chromosome}.SAIGE_OUT.SAIGE.markers.txt ' \
              '--LOCO=FALSE ' \
              '--is_output_moreDetails=TRUE ' \
              '--maxMissing=1 '
        if self._association_pack.is_binary:
            cmd = cmd + '--is_Firth_beta=TRUE'
        run_cmd(cmd, True, chromosome + ".SAIGE_markers.log")

        return chromosome

    def _process_saige_output(self, tarball_prefix: str, chromosome: str) -> pandas.DataFrame:

        # Load the raw table
        saige_table = pd.read_csv(tarball_prefix + "." + chromosome + ".SAIGE_OUT.SAIGE.gene.txt", sep='\t')
        saige_table = saige_table.rename(columns={'Region': 'ENST'})
        saige_table = saige_table.drop(columns=['Group', 'max_MAF'])

        # Get column names for Mask/MAF information if possible
        saige_table = self._define_field_names_from_tarball_prefix(tarball_prefix, saige_table)

        return saige_table

    def _annotate_saige_output(self, completed_gene_tables: list, completed_marker_chromosomes: list) -> list:

        saige_table = pd.concat(completed_gene_tables)

        # Now process the gene table into a useable format:
        # First read in the transcripts file
        transcripts_table = build_transcript_table()

        # Now merge the transcripts table into the gene table to add annotation and the write
        saige_table = pd.merge(transcripts_table, saige_table, on='ENST', how="left")
        with open(self._output_prefix + '.genes.SAIGE.stats.tsv', 'w') as gene_out:

            # Sort just in case
            saige_table = saige_table.sort_values(by=['chrom', 'start', 'end'])

            saige_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')
            gene_out.close()

            # And bgzip and tabix...
            cmd = "bgzip /test/" + self._output_prefix + '.genes.SAIGE.stats.tsv'
            run_cmd(cmd, True)
            cmd = "tabix -S 1 -s 2 -b 3 -e 4 /test/" + self._output_prefix + '.genes.SAIGE.stats.tsv.gz'
            run_cmd(cmd, True)

        outputs = [self._output_prefix + '.SAIGE_step1.log',
                   self._output_prefix + '.SAIGE_step2.log',
                   self._output_prefix + '.genes.SAIGE.stats.tsv.gz',
                   self._output_prefix + '.genes.SAIGE.stats.tsv.gz.tbi']

        if self._association_pack.run_marker_tests:

            variant_index = []
            saige_table_marker = []
            # Open all chromosome indicies and load them into a list and append them together
            for chromosome in completed_marker_chromosomes:
                variant_index.append(pd.read_csv(gzip.open(f'filtered_bgen/{chromosome}.filtered.vep.tsv.gz', 'rt'),
                                                 sep="\t",
                                                 dtype={'SIFT': str, 'POLYPHEN': str}))
                saige_table_marker.append(pd.read_csv(chromosome + ".SAIGE_OUT.SAIGE.markers.txt", sep="\t"))

            variant_index = pd.concat(variant_index)
            variant_index = variant_index.set_index('varID')

            saige_table_marker = pd.concat(saige_table_marker)

            # For markers, we can use the SNP ID column to get what we need
            saige_table_marker = saige_table_marker.rename(columns={'MarkerID': 'varID',
                                                                    'AC_Allele2': 'SAIGE_AC',
                                                                    'AF_Allele2': 'SAIGE_MAF'})
            saige_table_marker = saige_table_marker.drop(columns=['CHR', 'POS', 'Allele1', 'Allele2', 'MissingRate'])
            saige_table_marker = pd.merge(variant_index, saige_table_marker, on='varID', how="left")
            with open(self._output_prefix + '.markers.SAIGE.stats.tsv', 'w') as marker_out:
                # Sort by chrom/pos just to be sure...
                saige_table_marker = saige_table_marker.sort_values(by=['CHROM', 'POS'])

                saige_table_marker.to_csv(path_or_buf=marker_out, index=False, sep="\t", na_rep='NA')
                marker_out.close()

                # And bgzip and tabix...
                cmd = "bgzip /test/" + self._output_prefix + '.markers.SAIGE.stats.tsv'
                run_cmd(cmd, True)
                cmd = "tabix -S 1 -s 2 -b 3 -e 3 /test/" + self._output_prefix + '.markers.SAIGE.stats.tsv.gz'
                run_cmd(cmd, True)

            outputs.append([self._output_prefix + '.markers.SAIGE.stats.tsv.gz',
                            self._output_prefix + '.markers.SAIGE.stats.tsv.gz.tbi',
                            self._output_prefix + '.SAIGE_markers.log'])

        return outputs
