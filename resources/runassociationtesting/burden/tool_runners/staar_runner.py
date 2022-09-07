from os.path import exists
from typing import List

from runassociationtesting.association_resources import *
from runassociationtesting.burden.tool_runners.tool_runner import ToolRunner
from runassociationtesting.thread_utility import ThreadUtility


class STAARRunner(ToolRunner):

    def run_tool(self) -> None:

        # self._gene_infos = gene_infos
        chromosomes = set()

        # 1. Run the STAAR NULL model
        print("Running STAAR Null Model(s)...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A STAAR thread failed',
                                       incrementor=10,
                                       thread_factor=1)

        for phenoname in self._association_pack.pheno_names:
            thread_utility.launch_job(self._staar_null,
                                      phenoname=phenoname)
        thread_utility.collect_futures()

        # 2. Run the actual per-gene association tests
        print("Running STAAR masks * chromosomes...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A STAAR thread failed',
                                       incrementor=10,
                                       thread_factor=1)

        for phenoname in self._association_pack.pheno_names:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                # First 'if' is standard mode when running a regular burden test
                for chromosome in get_chromosomes(is_snp_tar=self._association_pack.is_snp_tar,
                                                  is_gene_tar=self._association_pack.is_gene_tar):
                    if exists(tarball_prefix + "." + chromosome + ".STAAR.matrix.rds"):
                        thread_utility.launch_job(self._staar_genes,
                                                  tarball_prefix=tarball_prefix,
                                                  chromosome=chromosome,
                                                  phenoname=phenoname)

        future_results = thread_utility.collect_futures()

        # 3. Print a preliminary STAAR output
        print("Finalising STAAR outputs...")
        completed_staar_files = []
        # And gather the resulting futures
        for result in future_results:
            tarball_prefix, finished_chromosome, phenoname = result
            completed_staar_files.append(f'{tarball_prefix}.{phenoname}.{finished_chromosome}.STAAR_results.tsv')

        # 4. Annotate and print final STAAR output
        self._outputs.extend(self._process_staar_outputs(completed_staar_files))

    # Generate the NULL model for STAAR
    def _staar_null(self, phenoname: str) -> None:

        # I have made a custom script in order to generate the STAAR Null model:
        # located in /usr/bin/runSTAAR_Null.R
        # This generates an RDS output file containing the NULL model
        # See the README.md for more information on these parameters
        cmd = "Rscript /prog/runSTAAR_Null.R " \
              "/test/phenotypes_covariates.formatted.txt " + \
              phenoname + " " + \
              str(self._association_pack.is_binary)
        if len(self._association_pack.found_quantitative_covariates) > 0:
            cmd += " " + ','.join(self._association_pack.found_quantitative_covariates)
        else:
            cmd += " NULL"
        if len(self._association_pack.found_categorical_covariates) > 0:
            cmd += " " + ','.join(self._association_pack.found_categorical_covariates)
        else:
            cmd += " NULL"
        run_cmd(cmd, True)

    # Run rare variant association testing using STAAR
    # Returns the finished chromosome to aid in output file creation
    @staticmethod
    def _staar_genes(tarball_prefix: str, chromosome: str, phenoname: str) -> tuple:

        # I have made a custom script in order to run the actual per-gene model in STAAR:
        # located in /usr/bin/runSTAAR_Genes.R
        # This generates a text output file of p.values
        # See the README.md for more information on these parameters
        cmd = f'Rscript /prog/runSTAAR_Genes.R ' \
                f'/test/{tarball_prefix}.{chromosome}.STAAR.matrix.rds ' \
                f'/test/{tarball_prefix}.{chromosome}.variants_table.STAAR.tsv ' \
                f'/test/{phenoname}.STAAR_null.rds ' + \
                f'{phenoname} ' \
                f'{tarball_prefix} ' \
                f'{chromosome} ' \
                f'none'

        run_cmd(cmd, True)

        return tarball_prefix, chromosome, phenoname

    # Helper method to just write STAAR outputs
    @staticmethod
    def _write_staar_csv(file_name: str, completed_staar_files: List[str]):
        with open(file_name, 'w') as staar_output:
            output_csv = csv.DictWriter(staar_output,
                                        delimiter="\t",
                                        fieldnames=['SNP', 'n.samps', 'pheno', 'relatedness.correction',
                                                    'staar.O.p', 'staar.SKAT.p', 'staar.burden.p',
                                                    'staar.ACAT.p', 'n_var', 'cMAC'])
            output_csv.writeheader()
            for file in completed_staar_files:
                with open(file, 'r') as curr_file_reader:
                    curr_file_csv = csv.DictReader(curr_file_reader, delimiter="\t")
                    for gene in curr_file_csv:
                        output_csv.writerow(gene)
                curr_file_reader.close()
            staar_output.close()

    # Process STAAR output files
    def _process_staar_outputs(self, completed_staar_files: List[str]) -> list:

        if self._association_pack.is_non_standard_tar:

            if self._association_pack.is_gene_tar:
                prefix = 'GENE'
            elif self._association_pack.is_snp_tar:
                prefix = 'SNP'
            else:
                raise dxpy.AppError('Somehow output is neither GENE or SNP from STAAR!')

            self._write_staar_csv(self._output_prefix + '.' + prefix + '.STAAR.stats.tsv', completed_staar_files)
            outputs = [self._output_prefix + '.' + prefix + '.STAAR.stats.tsv']
        else:
            # Here we are concatenating a temp file of each tsv from completed_staar_files:
            self._write_staar_csv(self._output_prefix + '.genes.STAAR.stats.temp', completed_staar_files)

            # Now read in the concatenated STAAR stats file:
            staar_table = pd.read_csv(open(self._output_prefix + '.genes.STAAR.stats.temp', 'r'), sep="\t")

            # Now process the gene table into a useable format:
            # First read in the transcripts file
            transcripts_table = build_transcript_table()

            # Test what columns we have in the 'SNP' field so we can name them...
            field_names = self._define_field_names_from_pandas(staar_table.iloc[0])
            staar_table[field_names] = staar_table['SNP'].str.split("-", expand=True)
            staar_table = staar_table.drop(columns=['SNP'])  # And drop the SNP column now that we have processed it

            # Now merge the transcripts table into the gene table to add annotation and the write
            staar_table = pd.merge(transcripts_table, staar_table, on='ENST', how="left")
            with open(self._output_prefix + '.genes.STAAR.stats.tsv', 'w') as gene_out:

                # Sort just in case
                staar_table = staar_table.sort_values(by=['chrom', 'start', 'end'])

                staar_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')
                gene_out.close()

                # And bgzip and tabix...
                cmd = "bgzip /test/" + self._output_prefix + '.genes.STAAR.stats.tsv'
                run_cmd(cmd, True)
                cmd = "tabix -S 1 -s 2 -b 3 -e 4 /test/" + self._output_prefix + '.genes.STAAR.stats.tsv.gz'
                run_cmd(cmd, True)

            outputs = [self._output_prefix + '.genes.STAAR.stats.tsv.gz',
                       self._output_prefix + '.genes.STAAR.stats.tsv.gz.tbi']

        return outputs
