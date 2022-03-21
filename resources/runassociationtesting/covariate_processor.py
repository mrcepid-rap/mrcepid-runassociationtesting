import os
import csv

from association_resources import *
from ingest_data import IngestData
from association_pack import AssociationPack

class CovariateProcessor:

    def __init__(self, ingested_data: IngestData, phenoname: str, categorical_covariates: str, quantitative_covariates: str,
                 sex: int, is_binary: bool, run_marker_tests: bool, output_prefix: str, mode: str):

        self._sex = sex
        self._mode = mode
        self._get_num_threads()

        self._genetics_samples = set()
        self._select_individuals(ingested_data.inclusion_found,
                                 ingested_data.exclusion_found)

        self._phenotypes = {}
        self._pheno_names = []
        self._process_phenotype(phenoname)

        # Process additional covariates (check if requested in the function):
        self._add_covars = {}
        self._found_quantitative_covariates = []
        self._found_categorical_covariates = []
        self._process_additional_covariates(ingested_data.additional_covariates_found,
                                            categorical_covariates,
                                            quantitative_covariates)

        self._create_covariate_file(ingested_data.additional_covariates_found)

        self.association_pack = AssociationPack(tarball_prefixes=ingested_data.tarball_prefixes,
                                                bgen_dict=ingested_data.bgen_dict,
                                                is_binary=is_binary,
                                                sex=sex,
                                                threads=self._threads,
                                                run_marker_tests=run_marker_tests,
                                                output_prefix=output_prefix,
                                                pheno_name=self._pheno_name,
                                                found_quantitative_covariates=self._found_quantitative_covariates,
                                                found_categorical_covariates=self._found_categorical_covariates)

    # Get number of cores available to the instance:
    def _get_num_threads(self) -> None:
        self._threads = os.cpu_count()
        print("{0:65}: {val}".format("Number of threads available", val = self._threads))

    # Three steps here:
    # 1. Get individuals we plan to include
    # 2. Exclude individuals not wanted in the analysis
    # 3. Get individuals that are POSSIBLE to include (they actually have WES) and only keep 'include' samples
    def _select_individuals(self, inclusion_found: bool, exclusion_found: bool) -> set:

        # Get a list of individuals that we ARE going to use
        include_samples = set()
        if inclusion_found is True:
            inclusion_file = open('INCLUSION.lst', 'r')
            for indv in inclusion_file:
                indv = indv.rstrip()
                include_samples.add(indv)

        # Get a list of individuals that we ARE NOT going to use
        exclude_samples = set()
        if exclusion_found is True:
            exclude_file = open('EXCLUSION.lst', 'r')
            for indv in exclude_file:
                indv = indv.rstrip()
                exclude_samples.add(indv)

        # Get individuals with genetic data
        # Remember! the genetic data has already been filtered to individuals with WES data.
        genetics_fam_file = open('genetics/UKBB_450K_Autosomes_QCd.fam', 'r')
        for line in genetics_fam_file:
            line = line.rstrip()
            fields = line.split()
            eid = fields[0]
            if inclusion_found == False and exclusion_found == False:
                self._genetics_samples.add(eid)
            elif inclusion_found == False and exclusion_found == True:
                if eid not in exclude_samples:
                    self._genetics_samples.add(eid)
            elif inclusion_found == True and exclusion_found == False:
                if eid in include_samples:
                    self._genetics_samples.add(eid)
            else:
                if eid in include_samples and eid not in exclude_samples:
                    self._genetics_samples.add(eid)

        print("{0:65}: {val}".format("Total samples after inclusion/exclusion lists applied", val = len(self._genetics_samples)))

    # This is a helper function for 'create_covariate_file()' that processes the phenotype file
    def _process_phenotype(self, phenoname: str) -> None:
        # Need to go through phenofile first and injest into a dictionary and get the name of the phenofield:
        dialect = csv.Sniffer().sniff(open('model_phenotypes.pheno', 'r').readline(), delimiters=[' ','\t'])
        pheno_reader = csv.DictReader(open('model_phenotypes.pheno', 'r'), delimiter=dialect.delimiter, skipinitialspace=True)
        field_names = pheno_reader.fieldnames
        if len(field_names) != 3:
            if self._mode == 'phewas':

            else:
                if phenoname is None:
                    raise RuntimeError("Pheno file has more than three columns and phenoname is not set!")
                else:
                    if phenoname in field_names:
                        self._pheno_names = phenoname
                    else:
                        raise RuntimeError("phenoname was not found in the provided phenofile!")
        else:
            if phenoname is None:
                for field in field_names:
                    if field != "FID" and field != "IID":
                        self._pheno_names = field
            else:
                if phenoname in field_names:
                    self._pheno_names = phenoname
                else:
                    raise RuntimeError("phenoname was not found in the provided phenofile!")

        if "FID" not in field_names and "IID" not in field_names:
            raise RuntimeError("Pheno file does not contain FID/IID fields!")

        for indv in pheno_reader:
            # Will spit out an error if a given sample does not have data
            if indv[self._pheno_name] is None:
                raise dxpy.AppError("Phenotype file has blank lines!")
            # Exclude individuals that have missing data (NA/NAN)
            elif indv[self._pheno_name].lower() != "na" and indv[self._pheno_name].lower() != "nan" and indv[self._pheno_name].lower() != "":
                self._phenotypes[indv['FID']] = indv[self._pheno_name]

    # This is a helper function for 'create_covariate_file()' that processes requested additional phenotypes
    def _process_additional_covariates(self, additional_covariates_found: bool, categorical_covariates: str, quantitative_covariates: str) -> tuple:

        if additional_covariates_found:
            dialect = csv.Sniffer().sniff(open('additional_covariates.covariates', 'r').readline(), delimiters=[' ','\t'])
            additional_covar_reader = csv.DictReader(open('additional_covariates.covariates', 'r'), delimiter=dialect.delimiter, skipinitialspace=True)
            field_names = list.copy(additional_covar_reader.fieldnames)

            # make sure the sample ID field is here and remove it from 'field_names' to help with iteration
            if 'FID' not in field_names and 'IID' not in field_names:
                raise dxpy.AppError("FID & IID column not found in provided covariates file!")
            else:
                field_names.remove('FID')
                field_names.remove('IID')

            # Now process & check the categorical/quantitative covariates lists and match it to field_names:
            if categorical_covariates is not None:
                categorical_covariates = categorical_covariates.split(',')
                for covar in categorical_covariates:
                    if covar in field_names:
                        self._found_categorical_covariates.append(covar)
                    else:
                        print("Provided categorical covariate %s not found in additional covariates file..." % covar)

            if quantitative_covariates is not None:
                quantitative_covariates = quantitative_covariates.split(',')
                for covar in quantitative_covariates:
                    if covar in field_names:
                        self._found_quantitative_covariates.append(covar)
                    else:
                        print("Provided quantitative covariate %s not found in additional covariates file..." % covar)

            # Throw an error if user provided covariates but none were found
            if (len(self._found_categorical_covariates) + len(self._found_quantitative_covariates)) == 0:
                raise dxpy.AppError('Additional covariate file provided but no additional covariates found based on covariate names provided...')

            for sample in additional_covar_reader:
                # First check no NAs/Blanks exist
                all_covars_found = True
                sample_dict = {}
                for field_name in (self._found_quantitative_covariates + self._found_categorical_covariates):
                    if sample[field_name].lower() == "na" or sample[field_name].lower() == "nan" or sample[field_name].lower() == "":
                        all_covars_found = False
                    else:
                        sample_dict[field_name] = sample[field_name]

                if all_covars_found == True:
                    self._add_covars[sample['IID']] = sample_dict

    # Do covariate processing and sample inclusion/exclusion
    def _create_covariate_file(self, additional_covariates_found: bool) -> dict:

        # Read the base covariates into this code that we want to analyse:
        # Formatting is weird to fit with other printing below...
        print("{0:65}: {val}".format("Phenotype", val = self._pheno_name))
        print("{0:65}: {val}".format("Default covariates included in model", val = ''))
        print("{pad:^5}{:60}: {val}".format("Quantitative", val = 'age, age^2, PC1..PC10',pad=' '))
        if (self._sex == 2):
            print("{pad:^5}{:60}: {val}".format("Categorical", val = 'sex, WES_batch',pad=' '))
        else:
            print("{pad:^5}{:60}: {val}".format("Categorical", val = 'WES_batch',pad=' '))
        if additional_covariates_found:
            print("{0:65}: {val}".format("Number of individuals with non-null additional covariates", val = len(self._add_covars)))
            print("{0:65}: {val}".format("Additional covariates included in model", val = ''))
            if len(self._found_quantitative_covariates) > 0:
                print("{pad:^5}{:60}: {val}".format("Quantitative", val = ', '.join(self._found_quantitative_covariates) ,pad=' '))
            else:
                print("{pad:^5}{:60}: {val}".format("Quantitative", val = 'None', pad=' '))
            if len(self._found_categorical_covariates) > 0:
                print("{pad:^5}{:60}: {val}".format("Categorical", val = ', '.join(self._found_categorical_covariates) ,pad=' '))
            else:
                print("{pad:^5}{:60}: {val}".format("Categorical", val = 'None' ,pad=' '))
        else:
            print("No additional covariates provided/found beyond defaults...")

        base_covar_reader = csv.DictReader(open('base_covariates.covariates', 'r'), delimiter="\t")
        indv_written = 0 # Just to count the number of samples we will analyse
        formatted_combo_file = open('phenotypes_covariates.formatted.txt', 'w') # SAIGE needs a combo file

        write_fields = ["FID", "IID"]
        write_fields = write_fields + ["PC%s" % (x) for x in range(1,41)]
        write_fields = write_fields + ["age", "age_squared", "sex", "wes_batch"]
        write_fields = write_fields + [self._pheno_name]
        # This doesn't matter to python if we didn't find additional covariates. A list of len() == 0 does not lengthen
        # the target list (e.g. 'write_fields')
        write_fields = write_fields + self._found_quantitative_covariates + self._found_categorical_covariates

        combo_writer = csv.DictWriter(formatted_combo_file,
                                      fieldnames = write_fields,
                                      quoting = csv.QUOTE_NONE,
                                      delimiter = " ",
                                      extrasaction='ignore')
        combo_writer.writeheader()

        # Need a list of included individuals ONLY:
        include_samples = open('SAMPLES_Include.txt', 'w')
        num_all_samples = 0
        na_pheno_samples = 0 # for checking number of individuals missing phenotype information
        for indv in base_covar_reader:
            if indv['22001-0.0'] != "NA": # need to exclude blank row individuals, eid is normally the only thing that shows up, so filter on sex
                indv_writer = {'FID': indv['eid'],
                               'IID': indv['eid']}
                for PC in range(1,41):
                    old_PC = "22009-0.%s" % (PC)
                    new_pc = "PC%s" % (PC)
                    indv_writer[new_pc] = indv[old_PC]
                indv_writer['age'] = int(indv['21003-0.0'])
                indv_writer['age_squared'] = indv_writer['age']**2
                indv_writer['sex'] = int(indv['22001-0.0'])
                indv_writer['wes_batch'] = indv['wes.batch']

                # Check if we found additional covariates and make sure this sample has non-null values
                write_sample = False
                if len(self._add_covars) > 0:
                    if indv['eid'] in self._add_covars:
                        write_sample = True
                        for covariate in self._add_covars[indv['eid']]:
                            indv_writer[covariate] = self._add_covars[indv['eid']][covariate]
                else:
                    write_sample = True

                # exclude based on sex-specific analysis if required:
                if indv['eid'] in self._genetics_samples:
                    num_all_samples += 1
                    if indv['eid'] in self._phenotypes and write_sample:
                        indv_writer[self._pheno_name] = self._phenotypes[indv['eid']]
                        if self._sex == 2:
                            indv_written += 1
                            combo_writer.writerow(indv_writer)
                            include_samples.write(indv['eid'] + "\n")
                        elif self._sex == indv_writer['sex']:
                            indv_written += 1
                            combo_writer.writerow(indv_writer)
                            include_samples.write(indv['eid'] + "\n")
                    else:
                        na_pheno_samples += 1

        formatted_combo_file.close()
        include_samples.close()

        # Generate a plink file to use that only has included individuals:
        cmd = "plink2 " \
              "--bfile /test/genetics/UKBB_450K_Autosomes_QCd --make-bed --keep-fam /test/SAMPLES_Include.txt " \
              "--out /test/genetics/UKBB_450K_Autosomes_QCd_WBA"
        run_cmd(cmd, True)

        # I have to do this to recover the sample information from plink
        cmd = "docker run -v /home/dnanexus/:/test/ egardner413/mrcepid-associationtesting plink2 " \
              "--bfile /test/genetics/UKBB_450K_Autosomes_QCd_WBA " \
              "--validate | grep samples"
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()

        # Print to ensure that total number of individuals is consistent between genetic and covariate/phenotype data
        print("{0:65}: {val}".format("Samples with covariates after include/exclude lists applied", val = num_all_samples))
        print("{0:65}: {val}".format("Number of individuals with NaN/NA phenotype information", val = na_pheno_samples))
        print("{0:65}: {val}".format("Number of individuals written to covariate/pheno file", val = indv_written))
        print("{0:65}: {val}".format("Plink individuals written", val = stdout.decode('utf-8').rstrip(' loaded from\n')))