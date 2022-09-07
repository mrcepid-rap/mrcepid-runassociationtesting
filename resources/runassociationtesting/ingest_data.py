from abc import ABC
from typing import List, Tuple, Set

from runassociationtesting.association_pack import AssociationPack, ProgramArgs
from runassociationtesting.association_resources import *


# This class is slightly different that in other applets I have designed; it handles ALL inputs rather
# than just external dependencies.
# Note that it is an abstract class (ABC) so that different modules can access the methods within this class, but
# add additional imports as required.
class IngestData(ABC):

    def __init__(self, parsed_options: ProgramArgs):

        self._parsed_options = parsed_options

        # Work our way through all the resources we need
        threads = self._get_num_threads()
        is_binary = parsed_options.is_binary
        sex = parsed_options.sex
        self._ingest_docker_file()
        self._ingest_transcript_index(parsed_options.transcript_index)
        pheno_files = self._ingest_phenofile(parsed_options.phenofile)
        additional_covariates_found = self._ingest_covariates(parsed_options.base_covariates,
                                                              parsed_options.covarfile)
        inclusion_found, exclusion_found = self._define_exclusion_lists(parsed_options.inclusion_list,
                                                                        parsed_options.exclusion_list)

        # Once all data is ingested, process the covariates/phenotypes into a single file of
        # individuals that we want to analyse
        # Set samples to use
        genetics_samples = self._select_individuals(inclusion_found, exclusion_found)

        # Define phenotype information
        phenotypes, pheno_names = self._process_phenotype(parsed_options.phenoname, pheno_files)

        # Process additional covariates (check if requested in the function)
        found_categorical_covariates, found_quantitative_covariates, add_covars = \
            self._process_additional_covariates(additional_covariates_found,
                                                parsed_options.categorical_covariates,
                                                parsed_options.quantitative_covariates)

        # Create the joint pheno/covariate file for testing
        self._create_covariate_file(genetics_samples=genetics_samples,
                                    phenotypes=phenotypes,
                                    pheno_names=pheno_names,
                                    additional_covariates_found=additional_covariates_found,
                                    found_quantitative_covariates=found_quantitative_covariates,
                                    found_categorical_covariates=found_categorical_covariates,
                                    add_covars=add_covars,
                                    sex=sex)

        # And build an object that will contain all the information we need to run some specified analysis
        self._association_pack = AssociationPack(pheno_files=pheno_files,
                                                 inclusion_found=inclusion_found,
                                                 exclusion_found=exclusion_found,
                                                 additional_covariates_found=additional_covariates_found,
                                                 is_binary=is_binary,
                                                 sex=sex,
                                                 threads=threads,
                                                 pheno_names=pheno_names,
                                                 found_quantitative_covariates=found_quantitative_covariates,
                                                 found_categorical_covariates=found_categorical_covariates)

    def get_association_pack(self):
        return self._association_pack

    def set_association_pack(self, objects: AssociationPack):
        self._association_pack = objects

    # Get number of cores available to the instance:
    @staticmethod
    def _get_num_threads() -> int:
        threads = os.cpu_count()
        print(f'{"Number of threads available":{65}}: {threads}')
        return threads

    # Bring our docker image into our environment so that we can run commands we need:
    @staticmethod
    def _ingest_docker_file() -> None:
        cmd = "docker pull egardner413/mrcepid-associationtesting:latest"
        run_cmd(cmd)

    # Get transcripts for easy annotation:
    @staticmethod
    def _ingest_transcript_index(transcript_index: dxpy.DXFile) -> None:
        dxpy.download_dxfile(transcript_index.get_id(), 'transcripts.tsv.gz')

    @staticmethod
    def _ingest_phenofile(phenofile: List[dxpy.DXFile]) -> List[str]:

        # There can be multiple phenotype files in some modes, so we iterate through the list of them here
        phenofiles = []
        for pheno in phenofile:
            curr_pheno_name = pheno.describe()['name']

            dxpy.download_dxfile(pheno.get_id(), curr_pheno_name)
            phenofiles.append(curr_pheno_name)

        return phenofiles

    # Get covariate/phenotype data:
    @staticmethod
    def _ingest_covariates(base_covariates: dxpy.DXFile, covarfile: dxpy.DXFile) -> bool:
        dxpy.download_dxfile(base_covariates.get_id(), 'base_covariates.covariates')

        # Check if additional covariates were provided:
        additional_covariates_found = False
        if covarfile is not None:
            dxpy.download_dxfile(covarfile.get_id(), 'additional_covariates.covariates')
            additional_covariates_found = True

        return additional_covariates_found

    # Get inclusion/exclusion sample lists
    @staticmethod
    def _define_exclusion_lists(inclusion_list: dxpy.DXFile, exclusion_list: dxpy.DXFile) -> Tuple[bool, bool]:
        inclusion_found, exclusion_found = False, False
        if inclusion_list is not None:
            dxpy.download_dxfile(inclusion_list.get_id(), 'INCLUSION.lst')
            inclusion_found = True
        if exclusion_list is not None:
            dxpy.download_dxfile(exclusion_list.get_id(), 'EXCLUSION.lst')
            exclusion_found = True

        return inclusion_found, exclusion_found

    # Define individuals based on exclusion/inclusion lists
    @staticmethod
    def _select_individuals(inclusion_found, exclusion_found) -> Set[str]:

        # A set of all possible samples to include in this analysis
        genetics_samples = set()

        # Three steps:
        # 1. Get a list of individuals that we ARE going to use
        include_samples = set()
        if inclusion_found is True:
            inclusion_file = open('INCLUSION.lst', 'r')
            for indv in inclusion_file:
                indv = indv.rstrip()
                include_samples.add(indv)

        # 2. Get a list of individuals that we ARE NOT going to use
        exclude_samples = set()
        if exclusion_found is True:
            exclude_file = open('EXCLUSION.lst', 'r')
            for indv in exclude_file:
                indv = indv.rstrip()
                exclude_samples.add(indv)

        # 3. Get individuals that are POSSIBLE to include (they actually have WES) and only keep 'include' samples
        # Remember! the genetic data has already been filtered to individuals with WES data.
        genetics_fam_file = open('genetics/UKBB_470K_Autosomes_QCd.fam', 'r')
        for line in genetics_fam_file:
            line = line.rstrip()
            fields = line.split()
            eid = fields[0]
            if inclusion_found is False and exclusion_found is False:
                genetics_samples.add(eid)
            elif inclusion_found is False and exclusion_found is True:
                if eid not in exclude_samples:
                    genetics_samples.add(eid)
            elif inclusion_found is True and exclusion_found is False:
                if eid in include_samples:
                    genetics_samples.add(eid)
            else:
                if eid in include_samples and eid not in exclude_samples:
                    genetics_samples.add(eid)

        print(f'{"Total samples after inclusion/exclusion lists applied":{65}}: {len(genetics_samples)}')
        return genetics_samples

    # This is a helper function for 'create_covariate_file()' that processes the phenotype file
    @staticmethod
    def _process_phenotype(pheno_name: str, pheno_files: List[str]) -> Tuple[dict, List[str]]:

        # Define ways to store...
        phenotypes = {}  # A dictionary of phenotypes for every given individual
        pheno_names = []  # A list of all possible pheno_name strings

        # Since we check for phewas mode above, this should only 'iterate' through multiple pheno_files if running in
        # phewas mode...
        # This for loop simply checks for phenotype names from the pheno file(s) and ingests into a list:
        for phenofile in pheno_files:
            dialect = csv.Sniffer().sniff(open(phenofile, 'r').readline(), delimiters=[' ', '\t'])
            pheno_reader = csv.DictReader(open(phenofile, 'r'), delimiter=dialect.delimiter, skipinitialspace=True)
            field_names = pheno_reader.fieldnames

            # Check to make sure we have a proper identifier
            if "FID" not in field_names and "IID" not in field_names:
                raise RuntimeError("Pheno file does not contain FID/IID fields!")

            # And ingest individual phenofields...
            curr_pheno_names = []

            # If phenoname not provided, then ingest all possible phenotypes
            if pheno_name is None:
                for field in field_names:
                    if field != "FID" and field != "IID":
                        curr_pheno_names.append(field)
            # Else just try and munge the one phenoname from the file
            else:
                if pheno_name in field_names:
                    curr_pheno_names.append(pheno_name)
                else:
                    raise RuntimeError("phenoname was not found in the provided phenofile!")

            # And then iterate through every sample in the pheno_file and add the information to our phenotypes
            # dictionary
            for indv in pheno_reader:
                # Will spit out an error if a given sample does not have data
                for pheno in curr_pheno_names:
                    if pheno not in phenotypes:
                        phenotypes[pheno] = {}
                    if indv[pheno] is None:
                        raise dxpy.AppError("Phenotype file has blank lines!")
                    # Exclude individuals that have missing data (NA/NAN)
                    elif indv[pheno].lower() != "na" and indv[pheno].lower() != "nan" and indv[pheno].lower() != "":
                        phenotypes[pheno][indv['FID']] = indv[pheno]

            # Finally add all pheno_name added from this file to the final list of pheno_name to run.
            pheno_names.extend(curr_pheno_names)

        return phenotypes, pheno_names

    # This is a helper function for 'create_covariate_file()' that processes requested additional phenotypes
    @staticmethod
    def _process_additional_covariates(additional_covariates_found: bool, categorical_covariates: str,
                                       quantitative_covariates: str) -> Tuple[List[str], List[str], dict]:

        found_categorical_covariates = []
        found_quantitative_covariates = []
        add_covars = {}

        if additional_covariates_found:
            dialect = csv.Sniffer().sniff(open('additional_covariates.covariates', 'r').readline(),
                                          delimiters=[' ', '\t'])
            additional_covar_reader = csv.DictReader(open('additional_covariates.covariates', 'r'),
                                                     delimiter=dialect.delimiter,
                                                     skipinitialspace=True)
            field_names = list.copy(additional_covar_reader.fieldnames)

            # make sure the sample ID field is here and remove it from 'field_names' to help with iteration
            if 'FID' not in field_names and 'IID' not in field_names:
                raise dxpy.AppError('FID & IID column not found in provided covariates file!')
            else:
                field_names.remove('FID')
                field_names.remove('IID')

            # Now process & check the categorical/quantitative covariates lists and match it to field_names:
            if categorical_covariates is not None:
                categorical_covariates = categorical_covariates.split(',')
                for covar in categorical_covariates:
                    if covar in field_names:
                        found_categorical_covariates.append(covar)
                    else:
                        print(f'Provided categorical covariate {covar} not found in additional covariates file...')

            if quantitative_covariates is not None:
                quantitative_covariates = quantitative_covariates.split(',')
                for covar in quantitative_covariates:
                    if covar in field_names:
                        found_quantitative_covariates.append(covar)
                    else:
                        print(f'Provided quantitative covariate {covar} not found in additional covariates file...')

            # Throw an error if user provided covariates but none were found
            if (len(found_categorical_covariates) + len(found_quantitative_covariates)) == 0:
                raise dxpy.AppError('Additional covariate file provided but no additional covariates found based on'
                                    ' covariate names provided...')

            for sample in additional_covar_reader:
                # First check no NAs/Blanks exist
                all_covars_found = True
                sample_dict = {}
                for field_name in (found_quantitative_covariates + found_categorical_covariates):
                    if sample[field_name].lower() == "na" or\
                            sample[field_name].lower() == "nan" or\
                            sample[field_name].lower() == "":
                        all_covars_found = False
                    else:
                        sample_dict[field_name] = sample[field_name]

                if all_covars_found:
                    add_covars[sample['IID']] = sample_dict

        return found_categorical_covariates, found_quantitative_covariates, add_covars

    # Do covariate processing and sample inclusion/exclusion
    @staticmethod
    def _create_covariate_file(genetics_samples: Set[str], phenotypes: dict, pheno_names: List[str],
                               additional_covariates_found: bool,
                               found_quantitative_covariates: List[str], found_categorical_covariates: List[str],
                               add_covars: dict, sex: int, ) -> None:

        # Read the base covariates into this code that we want to analyse:
        # Formatting is weird to fit with other printing below...
        print(f'{"Phenotype(s)":{65}}: {", ".join(pheno_names)}')
        print(f'{"Default covariates included in model":{65}}:')
        print(f'{" ":^{5}}{"Quantitative":{60}}: {"age, age^2, PC1..PC10"}')
        print(f'{" ":^{5}}{"Categorical":{60}}: {"sex, WES_batch" if sex == 2 else "WES_batch"}')
        if additional_covariates_found:
            print(f'{"Number of individuals with non-null additional covariates":{65}}: {len(add_covars)}')
            print(f'{"Additional covariates included in model":{65}}:')
            print(f'{" ":^{5}}{"Quantitative":{60}}: {", ".join(found_quantitative_covariates) if len(found_quantitative_covariates) > 0 else "None"}')
            print(f'{" ":^{5}}{"Categorical":{60}}: {", ".join(found_categorical_covariates) if len(found_categorical_covariates) > 0 else "None"}')
        else:
            print(f'No additional covariates provided/found beyond defaults...')

        base_covar_reader = csv.DictReader(open('base_covariates.covariates', 'r'), delimiter="\t")
        indv_written = 0  # Just to count the number of samples we will analyse
        formatted_combo_file = open('phenotypes_covariates.formatted.txt', 'w', newline='\n')  # SAIGE needs combo file

        write_fields = ["FID", "IID"]
        write_fields = write_fields + ["PC%s" % x for x in range(1, 41)]
        write_fields = write_fields + ["age", "age_squared", "sex", "wes_batch"]
        write_fields = write_fields + pheno_names
        # This doesn't matter to python if we didn't find additional covariates. A list of len() == 0 does not lengthen
        # the target list (e.g. 'write_fields')
        write_fields = write_fields + found_quantitative_covariates + found_categorical_covariates

        combo_writer = csv.DictWriter(formatted_combo_file,
                                      fieldnames=write_fields,
                                      quoting=csv.QUOTE_NONE,
                                      delimiter=" ",
                                      extrasaction='ignore',
                                      lineterminator='\n')
        combo_writer.writeheader()

        # Need a list of included individuals ONLY:
        # We write both a file for every other tool and a file for regenie at the same time just in case
        # the user requests a REGENIE run
        include_samples = open('SAMPLES_Include.txt', 'w')
        num_all_samples = 0
        na_pheno_samples = 0  # for checking number of individuals missing phenotype information
        for indv in base_covar_reader:
            # need to exclude blank row individuals, eid is normally the only thing that shows up, so filter on sex
            if indv['22001-0.0'] != "NA" and indv['eid'] in genetics_samples:
                indv_writer = {'FID': indv['eid'],
                               'IID': indv['eid']}
                for PC in range(1, 41):
                    old_pc = "22009-0.%s" % PC
                    new_pc = "PC%s" % PC
                    indv_writer[new_pc] = indv[old_pc]
                indv_writer['age'] = int(indv['21003-0.0'])
                indv_writer['age_squared'] = indv_writer['age']**2
                indv_writer['sex'] = int(indv['22001-0.0'])
                indv_writer['wes_batch'] = indv['wes.batch']
                num_all_samples += 1

                # Check if we found additional covariates and make sure this sample has non-null values
                found_covars = False
                if len(add_covars) > 0:
                    if indv['eid'] in add_covars:
                        found_covars = True
                        for covariate in add_covars[indv['eid']]:
                            indv_writer[covariate] = add_covars[indv['eid']][covariate]
                else:
                    found_covars = True

                found_phenos = False
                if len(pheno_names) == 1:
                    pheno = pheno_names[0]
                    if indv['eid'] in phenotypes[pheno]:
                        found_phenos = True
                        indv_writer[pheno] = phenotypes[pheno][indv['eid']]
                    else:
                        na_pheno_samples += 1
                else:
                    found_phenos = True  # Always true because we can exclude NAs later when running phewas
                    for pheno in pheno_names:
                        if indv['eid'] in phenotypes[pheno]:
                            indv_writer[pheno] = phenotypes[pheno][indv['eid']]
                        else:
                            indv_writer[pheno] = 'NA'

                # exclude based on sex-specific analysis if required:
                if found_covars and found_phenos:
                    if sex == 2:
                        indv_written += 1
                        combo_writer.writerow(indv_writer)
                        include_samples.write(indv['eid'] + "\n")
                    elif sex == indv_writer['sex']:
                        indv_written += 1
                        combo_writer.writerow(indv_writer)
                        include_samples.write(indv['eid'] + "\n")

        formatted_combo_file.close()
        include_samples.close()

        # Generate a plink file to use that only has included individuals:
        cmd = "plink2 " \
              "--bfile /test/genetics/UKBB_470K_Autosomes_QCd --make-bed --keep-fam /test/SAMPLES_Include.txt " \
              "--out /test/genetics/UKBB_470K_Autosomes_QCd_WBA"
        run_cmd(cmd, True)

        # I have to do this to recover the sample information from plink
        cmd = "docker run -v /home/dnanexus/:/test/ egardner413/mrcepid-associationtesting plink2 " \
              "--bfile /test/genetics/UKBB_470K_Autosomes_QCd_WBA " \
              "--validate | grep samples"
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()

        # Print to ensure that total number of individuals is consistent between genetic and covariate/phenotype data
        print(f'{"Samples with covariates after include/exclude lists applied":{65}}: {num_all_samples}')

        if len(pheno_names) == 1:
            print(f'{"Number of individuals with NaN/NA phenotype information":{65}}: {na_pheno_samples}')
        print(f'{"Number of individuals written to covariate/pheno file":{65}}: {indv_written}')
        print(f'{"Plink individuals written":{65}}: {stdout.decode("utf-8").rstrip(" loaded from")}')
