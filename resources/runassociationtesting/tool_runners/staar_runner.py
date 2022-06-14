import csv

from os.path import exists
from ..association_pack import AssociationPack
from ..association_resources import *
from ..thread_utility import ThreadUtility


class STAARRunner:

    def __init__(self, association_pack: AssociationPack, gene_infos: list = None):

        self._association_pack = association_pack
        self._gene_infos = gene_infos
        self._chromosomes = set()

        # Create a file of genes if genes_to_run is !none
        if self._gene_infos is not None:
            with open('staar.gene_list', 'w') as gene_list_file:
                for gene_info in self._gene_infos:
                    gene_list_file.write(gene_info.name + '\n')
                    if gene_info['chrom'] not in self._chromosomes:
                        self._chromosomes.add(gene_info['chrom'])
                gene_list_file.close()

        # 1. Run the STAAR NULL model
        print("Running STAAR Null Model(s)...")
        thread_utility = ThreadUtility(self._association_pack.threads, error_message='A STAAR thread failed',incrementor=10, thread_factor=1)
        for phenoname in self._association_pack.pheno_names:
            thread_utility.launch_job(self._staar_null,
                                      phenoname=phenoname)
        thread_utility.collect_futures()

        # 2. Run the actual per-gene association tests
        print("Running STAAR masks * chromosomes...")
        thread_utility = ThreadUtility(self._association_pack.threads, error_message='A STAAR thread failed', incrementor=10, thread_factor=1)

        for phenoname in self._association_pack.pheno_names:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                # First 'if' is standard mode when running a regular burden test
                if self._gene_infos is None:
                    for chromosome in get_chromosomes(is_snp_tar=self._association_pack.is_snp_tar,is_gene_tar=self._association_pack.is_gene_tar):
                        if exists(tarball_prefix + "." + chromosome + ".STAAR.matrix.rds"):
                            thread_utility.launch_job(self._staar_genes,
                                                      tarball_prefix = tarball_prefix,
                                                      chromosome = chromosome,
                                                      phenoname = phenoname)
                # Second 'if' is when we extract information for specific genes
                else:
                    for chromosome in self._chromosomes:
                        if exists(tarball_prefix + "." + chromosome + ".STAAR.matrix.rds"):
                            thread_utility.launch_job(self._staar_genes,
                                                      tarball_prefix = tarball_prefix,
                                                      chromosome = chromosome,
                                                      phenoname = phenoname)

        future_results = thread_utility.collect_futures()

        # 3. Print a preliminary STAAR output
        print("Finalising STAAR outputs...")
        self._completed_staar_files = []
        # And gather the resulting futures
        for result in future_results:
            tarball_prefix, finished_chromosome, phenoname = result
            self._completed_staar_files.append(tarball_prefix + "." + phenoname + "." + finished_chromosome + ".STAAR_results.tsv")

        # 4. Annotate and print final STAAR output
        self.outputs = self._process_staar_outputs()

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
    def _staar_genes(self, tarball_prefix: str, chromosome: str, phenoname: str) -> tuple:

        # I have made a custom script in order to run the actual per-gene model in STAAR:
        # located in /usr/bin/runSTAAR_Genes.R
        # This generates a text output file of p.values
        # See the README.md for more information on these parameters
        cmd = "Rscript /prog/runSTAAR_Genes.R " \
                "/test/" + tarball_prefix + "." + chromosome + ".STAAR.matrix.rds " \
                "/test/" + tarball_prefix + "." + chromosome + ".variants_table.STAAR.tsv " \
                "/test/" + phenoname + ".STAAR_null.rds " + \
                phenoname + " " + \
                tarball_prefix + " " + \
                chromosome + " "
        # If we only want to analyse a subset of genes, define here:
        if self._gene_infos is None:
            cmd = cmd + "none"
        else:
            cmd = cmd + "/test/staar.gene_list"

        run_cmd(cmd, True)

        return tarball_prefix, chromosome, phenoname

    # Process STAAR output files
    def _process_staar_outputs(self) -> list:

        if self._association_pack.is_non_standard_tar:

            if self._association_pack.is_gene_tar:
                prefix = 'GENE'
            elif self._association_pack.is_snp_tar:
                prefix = 'SNP'

            with open(self._association_pack.output_prefix + '.' + prefix + '.STAAR.stats.tsv', 'w') as staar_output:
                output_csv = csv.DictWriter(staar_output,
                                            delimiter = "\t",
                                            fieldnames=['SNP','n.samps','pheno','relatedness.correction','staar.O.p','staar.SKAT.p','staar.burden.p',
                                                        'staar.ACAT.p','n_var','cMAC'])
                output_csv.writeheader()
                # In this section (e.g. SNP/GENE tars) this loop technically doesn't need to exist, but keeping here
                # for similarity to section below for regular tars
                for file in self._completed_staar_files:
                    with open(file, 'r') as curr_file_reader:
                        curr_file_csv = csv.DictReader(curr_file_reader, delimiter = "\t")
                        for gene in curr_file_csv:
                            output_csv.writerow(gene)
                    curr_file_reader.close()
                staar_output.close()

            outputs = [self._association_pack.output_prefix + '.' + prefix + '.STAAR.stats.tsv']
        else:
            # Here we are concatenating a temp file of each tsv from completed_staar_files:
            with open(self._association_pack.output_prefix + '.genes.STAAR.stats.temp', 'w') as staar_output:
                output_csv = csv.DictWriter(staar_output,
                                            delimiter = "\t",
                                            fieldnames=['SNP','n.samps','pheno','relatedness.correction','staar.O.p','staar.SKAT.p','staar.burden.p',
                                                        'staar.ACAT.p','n_var','cMAC'])
                output_csv.writeheader()
                for file in self._completed_staar_files:
                    with open(file, 'r') as curr_file_reader:
                        curr_file_csv = csv.DictReader(curr_file_reader, delimiter = "\t")
                        for gene in curr_file_csv:
                            output_csv.writerow(gene)
                    curr_file_reader.close()
                staar_output.close()

            # Now read in the concatenated STAAR stats file:
            staar_table = pd.read_csv(open(self._association_pack.output_prefix + '.genes.STAAR.stats.temp', 'r'), sep = "\t")

            # Now process the gene table into a useable format:
            # First read in the transcripts file
            transcripts_table = pd.read_csv(gzip.open('transcripts.tsv.gz', 'rt'), sep = "\t")
            transcripts_table = transcripts_table.rename(columns={'#chrom':'chrom'})
            transcripts_table = transcripts_table.set_index('ENST')
            transcripts_table = transcripts_table[transcripts_table['fail'] == False]
            transcripts_table = transcripts_table.drop(columns=['syn.count','fail.cat','fail'])

            # Limit to genes we care about if running only a subset:
            if self._gene_infos is not None:
                valid_genes = []
                for gene_info in self._gene_infos:
                    valid_genes.append(gene_info.name)
                transcripts_table = transcripts_table.loc[valid_genes]

            # Test what columns we have in the 'SNP' field so we can name them...
            field_one = staar_table.iloc[0]
            field_one = field_one['SNP'].split("-")
            field_names = ['ENST']
            if len(field_one) == 2: # This is the bare minimum, always name first column ENST, and second column 'var1'
                field_names.append('var1')
            elif len(field_one) == 3: # This could be the standard naming format... check that column [2] is MAF/AC
                if 'MAF' in field_one[2] or 'AC' in field_one[2]:
                    field_names.extend(['MASK','MAF'])
                else: # This means we didn't hit on MAF in column [2] and a different naming convention is used...
                    field_names.extend(['var1','var2'])
            else:
                for i in range(2,len(field_one) + 1):
                    field_names.append('var%i' % i)

            # Process the 'SNP' column into separate fields and remove
            staar_table[field_names] = staar_table['SNP'].str.split("-",expand = True)
            staar_table = staar_table.drop(columns=['SNP'])

            # Now merge the transcripts table into the gene table to add annotation and the write
            staar_table = pd.merge(transcripts_table, staar_table, on='ENST', how="left")
            with open(self._association_pack.output_prefix + '.genes.STAAR.stats.tsv', 'w') as gene_out:

                # Sort just in case
                staar_table = staar_table.sort_values(by=['chrom','start','end'])

                staar_table.to_csv(path_or_buf=gene_out, index = False, sep="\t", na_rep='NA')
                gene_out.close()

                # And bgzip and tabix...
                cmd = "bgzip /test/" + self._association_pack.output_prefix + '.genes.STAAR.stats.tsv'
                run_cmd(cmd, True)
                cmd = "tabix -S 1 -s 2 -b 3 -e 4 /test/" + self._association_pack.output_prefix + '.genes.STAAR.stats.tsv.gz'
                run_cmd(cmd, True)

            outputs = [self._association_pack.output_prefix + '.genes.STAAR.stats.tsv.gz',
                       self._association_pack.output_prefix + '.genes.STAAR.stats.tsv.gz.tbi']

        return outputs
