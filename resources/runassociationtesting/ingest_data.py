import os
import csv
import dxpy

from abc import ABC
from typing import Set, Tuple, List, Any, Dict
from pathlib import Path

from general_utilities.association_resources import run_cmd, download_dxfile_by_name
from general_utilities.mrc_logger import MRCLogger

from runassociationtesting.association_pack import AssociationPack, ProgramArgs


class IngestData(ABC):
    """Download and process required files and data to enable module functionality.

    This is an abstract class (ABC) so that different modules can access the methods within this class,
    but add additional imports as required in their own implementations. This interface does all of the processing of
    files required by all modules that use the associationtesting framework. This includes:

    1. The phenotype file
    2. The covariate file(s) – base and additional (if provided)
    3. Sample inclusion / exclusion lists
    4. Gene transcript dictionary

    :param parsed_options: A subclass of ProgramArgs containing options parsed by argparse
    """

    def __init__(self, parsed_options: ProgramArgs):

        # Initiate logger – This can also be used by a class which implements this Interface
        self._logger = MRCLogger(__name__).get_logger()

        self._parsed_options = parsed_options

        # Grab the Docker image, so we can run tools not on the DNANexus platform by default
        self._ingest_docker_file()

        # Work our way through all the resources we need
        # Gene transcript dictionary
        self._ingest_transcript_index(parsed_options.transcript_index)

        # Phenotype files
        phenotypes = self._ingest_phenofile(parsed_options.phenofile, parsed_options.phenoname)
        pheno_names = list(phenotypes.keys())  # A list of pheno_names since we don't store phenotypes at runtime

        # Base and Additional covariates
        additional_covariates_found = self._ingest_covariates(parsed_options.base_covariates,
                                                              parsed_options.covarfile)

        # Sample inclusion / exclusion lists
        inclusion_found, exclusion_found = self._define_exclusion_lists(parsed_options.inclusion_list,
                                                                        parsed_options.exclusion_list)

        # Once all data is ingested, process the covariates/phenotypes into a single file of individuals that we want
        # to analyse
        genetics_samples = self._select_individuals(inclusion_found, exclusion_found)

        # Process additional covariates (check if requested in the function)
        found_categorical_covariates, found_quantitative_covariates, add_covars = \
            self._process_additional_covariates(additional_covariates_found,
                                                parsed_options.categorical_covariates,
                                                parsed_options.quantitative_covariates)

        # After all sample/phenotype/covariate processing, create the joint pheno/covariate file for testing
        self._create_covariate_file(genetics_samples=genetics_samples,
                                    phenotypes=phenotypes,
                                    pheno_names=pheno_names,
                                    additional_covariates_found=additional_covariates_found,
                                    found_categorical_covariates=found_categorical_covariates,
                                    found_quantitative_covariates=found_quantitative_covariates,
                                    add_covars=add_covars,
                                    sex=parsed_options.sex)

        # And build an object that will contain all the information we need to run some specified analysis
        self._association_pack = AssociationPack(is_binary=parsed_options.is_binary,
                                                 sex=parsed_options.sex,
                                                 threads=self._get_num_threads(),
                                                 pheno_names=pheno_names,
                                                 found_quantitative_covariates=found_quantitative_covariates,
                                                 found_categorical_covariates=found_categorical_covariates)

    def get_association_pack(self) -> AssociationPack:
        """Getter for self._association_pack

        :return: self._association_pack, which contains stored parameters necessary for running the current module
        """
        return self._association_pack

    def set_association_pack(self, objects: AssociationPack) -> None:
        """Setter for self._association_pack

        :param objects: An AssociationPack or subclass implementing AssociationPack
        :return: None
        """
        self._association_pack = objects

    def _get_num_threads(self) -> int:
        """Get number of cores available to the instance using :func:`os.cpu_count()`

        :return: The number of CPUs on this instance
        """

        threads = os.cpu_count()
        self._logger.info(f'{"Number of threads available":{65}}: {threads}')

        return threads

    @staticmethod
    def _ingest_docker_file() -> None:
        """Download the default Docker image so that we can run tools not on the DNANexus platform.

        :return: None
        """

        cmd = "docker pull egardner413/mrcepid-burdentesting:latest"
        run_cmd(cmd, is_docker=False)

    @staticmethod
    def _ingest_transcript_index(transcript_index: dxpy.DXFile) -> None:
        """Get transcripts for gene annotation

        The provided file with always be placed at `$HOME/transcripts.tsv.gz`.

        :param transcript_index: A transcript index in .tsv format. For more information on the structure of the file,
            see the README.
        :return: None
        """

        dxpy.download_dxfile(transcript_index.get_id(), 'transcripts.tsv.gz')

    def _ingest_phenofile(self, pheno_files: List[dxpy.DXFile], pheno_name: str) -> Dict[str, Dict[str, Any]]:
        """Download provided phenotype files and attempt to retrieve phenotypes from these file(s)

        There can be multiple phenotype files in some modes, so this method will iterate through the provided list
        and download all of them and store the respective file name in a List. This method with then process a (
        possible list of) phenotype file(s) (pheno_files) of the format:

                FID IID pheno1 pheno2 ... pheno3

        where a pheno_file can be either space or tab-delimited. This file **cannot** have blank phenotypes,
        but can have NA/NAN (in any case). This method will then parse a provided phenotype file and either:

        1. Extract all putative phenotypes (other than FID/IID) from pheno_file if pheno_name is NOT provided
        2. Extract a single provided pheno_name from a given pheno_file, else throw an error that the name wasn't found.

        Note that option (2) does not work when multiple pheno_files are provided and will raise an exception. Raw
        values from the provided pheno_file(s) are then extracted and placed into a dict() of phenotypes with a
        format of:

                {'pheno_name': {'1234567': 0.1,
                                '2345678': 0.3}}

        :param pheno_name: Phenotype name that MUST be in the header (csv.fieldnames) of the provided file.
        :param pheno_files: List of DXFile(s) pointing to provided pheno_files.
        :return: A dict with keys of pheno_name(s) found and values a dict with keys of IIDs and values of the
            given phenotype
        """

        # A dictionary of phenotypes for every individual in the phenofile
        phenotypes = {}

        for dx_pheno_file in pheno_files:
            pheno_file = download_dxfile_by_name(dx_pheno_file)

            total_missing_dict = {}
            total_samples = 0

            # Now process the downloaded pheno_file and extract phenotype names/raw phenotypes
            dialect = csv.Sniffer().sniff(pheno_file.open('r').readline(), delimiters=[' ', '\t'])
            with pheno_file.open('r') as pheno_reader:
                pheno_csv = csv.DictReader(pheno_reader, delimiter=dialect.delimiter)
                field_names = pheno_csv.fieldnames

                # Check to make sure we have a proper identifier
                if "FID" not in field_names or "IID" not in field_names:
                    raise ValueError("Pheno file does not contain FID/IID fields!")

                # And ingest individual phenofields...
                curr_pheno_names = []

                # If phenoname not provided, then ingest all possible phenotypes
                if pheno_name is None:
                    for field in field_names:
                        if field != "FID" and field != "IID":
                            curr_pheno_names.append(field)
                            total_missing_dict[field] = 0
                # Else just try and munge the one phenoname from the file
                else:
                    if pheno_name in field_names:
                        curr_pheno_names.append(pheno_name)
                        total_missing_dict[pheno_name] = 0
                    else:
                        raise ValueError("phenoname was not found in the provided phenofile!")

                # And then iterate through every sample in the pheno_file and add the information to our phenotypes
                # dictionary
                for indv in pheno_csv:
                    total_samples += 1
                    for pheno in curr_pheno_names:
                        if pheno not in phenotypes:
                            phenotypes[pheno] = {}

                        # Will spit out an error if a given sample is not formatted properly
                        if indv[pheno] is None:
                            raise ValueError("Phenotype file has blank lines!")

                        # Exclude individuals that have missing data (NA/NAN) from the dictionary (these will be
                        # properly added to the final covariate file as Null values later in the data ingestion
                        # process).
                        elif indv[pheno].lower() != "na" and indv[pheno].lower() != "nan" and indv[pheno].lower() != "":
                            phenotypes[pheno][indv['FID']] = indv[pheno]
                        else:
                            total_missing_dict[pheno] += 1

            for pheno in total_missing_dict:
                print_string = f'Total samples in file {pheno_file.name}'
                self._logger.info(f'{print_string:{65}}: {total_samples}')
                print_string = f'Phenotype "{pheno}" samples missing'
                prop_missing = (total_missing_dict[pheno] / total_samples) * 100
                self._logger.info(f'{" ":^{5}}{print_string:{60}}: {total_missing_dict[pheno]} ({prop_missing:0.2f}%)')

        return phenotypes

    @staticmethod
    def _ingest_covariates(base_covariates: dxpy.DXFile, covarfile: dxpy.DXFile) -> bool:
        """Download covariate data

        Base covariates will always be placed at `$HOME/base_covariates.covariates`. Additional covariates (if provided)
        will always be placed at `$HOME/additional_covariates.covariates`.

        :param base_covariates: A DXFile pointing to the base_covariates file on the RAP
        :param covarfile: A DXFile (possibly None if not provided) pointing to the base_covariates file on the RAP
        :return: A boolean that is true if additional covariates beyond the base covariates were provided
        """

        dxpy.download_dxfile(base_covariates.get_id(), 'base_covariates.covariates')

        # Check if additional covariates were provided:
        additional_covariates_found = False
        if covarfile is not None:
            dxpy.download_dxfile(covarfile.get_id(), 'additional_covariates.covariates')
            additional_covariates_found = True

        return additional_covariates_found

    @staticmethod
    def _define_exclusion_lists(inclusion_list: dxpy.DXFile, exclusion_list: dxpy.DXFile) -> Tuple[bool, bool]:
        """Get inclusion/exclusion sample lists

        If provided, inclusion and exclusion lists will be downloaded to `$HOME/INCLUSION.lst` and
        `$HOME/EXCLUSION.list`, respectively.

        :param inclusion_list: A DXFile (possibly None if not provided) pointing to a sample inclusion list file on
            the RAP
        :param exclusion_list: A DXFile (possibly None if not provided) pointing to a sample inclusion list file on
            the RAP
        :return: A Tuple with two booleans, describing if an inclusion or exclusion list were found, respectively
        """

        inclusion_found, exclusion_found = False, False
        if inclusion_list is not None:
            dxpy.download_dxfile(inclusion_list.get_id(), 'INCLUSION.lst')
            inclusion_found = True
        if exclusion_list is not None:
            dxpy.download_dxfile(exclusion_list.get_id(), 'EXCLUSION.lst')
            exclusion_found = True

        return inclusion_found, exclusion_found

    def _select_individuals(self, inclusion_found, exclusion_found) -> Set[str]:
        """Define individuals based on exclusion/inclusion lists.

        Three steps to this:

        1. Read inclusion samples
        2. Read exclusion samples
        3. Read samples with WES & Genetic data (by default encoded as the 'genetics_qc_pass' in the base_covariates file)
        4. Mash these three sources together to get a `set` of individuals we will include in downstream analysis

        The logic to step (4) is:

        * If inclusion + exclusion: Any individual to include in our study MUST be in the inclusion list and NOT in the exclusion list if both are provided.
        * If inclusion: a sample MUST be in this file to be included.
        * If exclusion: a sample CANNOT be in the file and be included.
        * If neither: All samples in base_covariates with genetics_qc_pass = PASS are included

        :param inclusion_found: Was an inclusion list found?
        :param exclusion_found: Was an exclusion list found?
        :return: A set() containing eids of the samples for all downstream analysis
        """

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

        # 3. Get individuals that are POSSIBLE to include (they actually have WES AND pass genetic data filtering) and
        # only keep 'include' samples.
        # Remember! the genetic data has already been filtered to individuals with WES data.

        # 4. Generate A set of all possible samples to include in this analysis
        genetics_samples = set()

        with Path('base_covariates.covariates').open('r') as base_covariates_file:
            base_covar_reader = csv.DictReader(base_covariates_file, delimiter="\t")
            for indv in base_covar_reader:
                eid = indv['eid']
                genetics_status = int(indv['genetics_qc_pass'])

                # extract the variable indicating whether an individual has passed genetic data QC:
                if genetics_status == 1:
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
            base_covariates_file.close()

        self._logger.info(f'{"Total samples after inclusion/exclusion lists applied":{65}}: {len(genetics_samples)}')
        return genetics_samples

    def _process_additional_covariates(self, additional_covariates_found: bool, categorical_covariates: List[str],
                                       quantitative_covariates: List[str]
                                       ) -> Tuple[List[str], List[str], Dict[str, Dict[str, Any]]]:
        """Identify additional covariates beyond base covariates and ensure all asked for additional covariate(s) are
        in the provided covariate file.

        Process a covariate file of the form:

                FID IID catCovar1 catCovar2 quantCovar3

        where the provided file can be either space or tab-delimited. This method is relatively, simple, and will
        iterate through the provided covariate file, ensuring that covariates that are asked for are in the provided
        file. This method will also report statistics on total / proportion of individuals missing each covariate.

        :param additional_covariates_found: Was an additional covariates file provided on the command-line?
        :param categorical_covariates: List of names of additional categorical covariates
        :param quantitative_covariates: List of names of additional quantitative covariates
        :return: A Tuple containing two lists (of categorical and quantitative covariate names, respectively) and a Dict
            with keys of sample IDs and values a Dict with keys of covariate names and values of covariate values.
        """

        self._logger.info('Processing additional covariates...')

        found_categorical_covariates = []
        found_quantitative_covariates = []
        add_covars = {}

        if additional_covariates_found:

            # Keep track of sample totals:
            total_missing_dict = {}

            additional_covariates_file = Path('additional_covariates.covariates')

            # Determine delimiter of provided file (space or tab)
            dialect = csv.Sniffer().sniff(additional_covariates_file.open('r').readline(), delimiters=[' ', '\t'])

            with additional_covariates_file.open('r') as additional_covariates_reader:
                additional_covar_csv = csv.DictReader(additional_covariates_reader,
                                                      delimiter=dialect.delimiter)
                field_names = list.copy(additional_covar_csv.fieldnames)

                # make sure the sample ID field is here and remove it from 'field_names' to help with iteration
                if 'FID' not in field_names or 'IID' not in field_names:
                    raise ValueError('FID & IID column not found in provided covariates file!')
                else:
                    field_names.remove('FID')
                    field_names.remove('IID')

                # Now process & check the categorical/quantitative covariates lists and match it to field_names:
                if categorical_covariates is not None:
                    for covar in categorical_covariates:
                        if covar in field_names:
                            found_categorical_covariates.append(covar)
                            total_missing_dict[covar] = 0
                        else:
                            self._logger.warning(f'Provided categorical covariate {covar} not found in '
                                                 f'additional covariates file...')

                if quantitative_covariates is not None:
                    for covar in quantitative_covariates:
                        if covar in field_names:
                            found_quantitative_covariates.append(covar)
                            total_missing_dict[covar] = 0
                        else:
                            self._logger.warning(f'Provided quantitative covariate {covar} not found in '
                                                 f'additional covariates file...')

                # Throw an error if user provided covariates but none were found
                if (len(found_categorical_covariates) + len(found_quantitative_covariates)) == 0:
                    raise RuntimeError('Additional covariate file provided but no additional covariates found based on '
                                       'covariate names provided...')

                total_samples = 0
                for sample in additional_covar_csv:
                    # First check no NAs/Blanks exist for a given sample. i.e., a sample must have ALL asked for
                    # covariates to be included in the final analysis
                    all_covars_found = True
                    sample_dict = {}
                    total_samples += 1
                    for covar_name in (found_quantitative_covariates + found_categorical_covariates):
                        if sample[covar_name] is None:
                            all_covars_found = False
                            total_missing_dict[covar_name] += 1
                        elif sample[covar_name].lower() in ['na', 'nan', '']:
                            all_covars_found = False
                            total_missing_dict[covar_name] += 1
                        else:
                            sample_dict[covar_name] = sample[covar_name]

                    if all_covars_found:
                        add_covars[sample['IID']] = sample_dict

                # Report total number of samples with missing additional covariates
                for covar_name in total_missing_dict:
                    print_string = f'Add covar "{covar_name}" samples missing'
                    prop_missing = (total_missing_dict[covar_name] / total_samples) * 100
                    self._logger.info(f'{print_string:{65}}: {total_missing_dict[covar_name]} ({prop_missing:0.2f}%)')

        return found_categorical_covariates, found_quantitative_covariates, add_covars

    def _create_covariate_file(self, genetics_samples: Set[str], phenotypes: Dict[str, Dict[str, Any]],
                               pheno_names: List[str], additional_covariates_found: bool,
                               found_categorical_covariates: List[str], found_quantitative_covariates: List[str],
                               add_covars: Dict[str, Dict[str, Any]], sex: int) -> None:

        """Print final covariate + phenotype file while implementing sample inclusion/exclusion

        While long, this method is relatively simple. Simple checking that all covariates and phenotypes for a given
        individual are present, and then writing that individual to the final covariates file. If multiple phenotypes
        are provided, individuals with missing phenotypes are coded as NA. Otherwise, individuals are excluded
        from the final file. For covariates, individuals must have ALL asked for and base covariates to to be written
        to the final pheno + covariate file (always located at `phenotypes_covariates.formatted.txt`).

        :param genetics_samples: List of valid sample IDs as determined by the combination of inclusion / exclusion
            lists and base covariate file(s)
        :param phenotypes: Dictionary of phenotypes with per-individual values
        :param pheno_names: List of pheno_names in phenotypes
        :param additional_covariates_found: Was an additional covariates file provided on the command-line?
        :param found_categorical_covariates: List of names of additional categorical covariates
        :param found_quantitative_covariates: List of names of additional quantitative covariates
        :param add_covars: Dictionary of additional covariates with per-individual values
        :param sex: Sex to restrict analysis to (either 0 [female], 1 [male], or 2 [both])
        :return: None
        """

        # Print some statistics about what we found in previous ingestion classes:
        self._logger.info(f'{"Phenotype(s)":{65}}: {", ".join(pheno_names)}')
        self._logger.info(f'{"Default covariates included in model":{65}}:')
        self._logger.info(f'{" ":^{5}}{"Quantitative":{60}}: {"age, age^2, PC1..PC10"}')
        self._logger.info(f'{" ":^{5}}{"Categorical":{60}}: {"sex, WES_batch" if sex == 2 else "WES_batch"}')
        if additional_covariates_found:
            self._logger.info(f'{"Number of individuals with non-null additional covariates":{65}}: {len(add_covars)}')
            self._logger.info(f'{"Additional covariates included in model":{65}}:')
            self._logger.info(f'{" ":^{5}}{"Quantitative":{60}}: {", ".join(found_quantitative_covariates) if len(found_quantitative_covariates) > 0 else "None"}')
            self._logger.info(f'{" ":^{5}}{"Categorical":{60}}: {", ".join(found_categorical_covariates) if len(found_categorical_covariates) > 0 else "None"}')
        else:
            self._logger.info(f'No additional covariates provided/found beyond defaults...')

        base_covariates_file = Path('base_covariates.covariates')
        final_covariates_file = Path('phenotypes_covariates.formatted.txt')

        with base_covariates_file.open('r') as base_covar_reader,\
                final_covariates_file.open('w', newline='\n') as final_covariates_writer:

            base_covar_csv = csv.DictReader(base_covar_reader, delimiter="\t")

            # Build fieldnames for final_covariates_file
            write_fields = ['FID', 'IID']
            write_fields = write_fields + ['PC%s' % x for x in range(1, 41)]
            write_fields = write_fields + ['age', 'age_squared', 'sex', 'wes_batch', 'array_batch']
            write_fields = write_fields + pheno_names
            # This doesn't matter to python if we didn't find additional covariates. A list of len() == 0 does not
            # lengthen the target list (e.g. 'write_fields')
            write_fields = write_fields + found_quantitative_covariates + found_categorical_covariates

            # Open the final CSV writer
            combo_writer = csv.DictWriter(final_covariates_writer,
                                          fieldnames=write_fields,
                                          quoting=csv.QUOTE_NONE,
                                          delimiter=' ',
                                          extrasaction='ignore',
                                          lineterminator='\n')
            indv_written = 0  # Just to count the number of samples we will analyse
            combo_writer.writeheader()

            # Need a list of included individuals ONLY:
            # We write both a file for every other tool and a file for regenie at the same time just in case
            # the user requests a REGENIE run
            with Path('SAMPLES_Include.txt').open('w') as include_samples:
                num_all_samples = 0
                for indv in base_covar_csv:
                    # need to exclude blank row individuals, eid is normally the only thing that shows up, so filter
                    # on sex
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
                        indv_writer['array_batch'] = indv['22000-0.0']
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

        # Print to ensure that total number of individuals is consistent between genetic and covariate/phenotype data
        self._logger.info(f'{"Samples with covariates after include/exclude lists applied":{65}}: {num_all_samples}')
        self._logger.info(f'{"Number of individuals written to covariate/pheno file":{65}}: {indv_written}')
