import tarfile
from os.path import exists

from abc import ABC
from association_resources import *


# This class is slightly different that in other applets I have designed; it handles ALL inputs rather
# than just external dependencies.
# Note that it is an abstract class (ABC) so that different modules can access the methods within this class, but
# add additional imports as required.
class IngestData(ABC):

    def __init__(self, parsed_options: dict, required_options: set):

        self._parsed_options = parsed_options

        # Set reasonable defaults
        self.additional_covariates_found = False
        self.inclusion_found = False
        self.exclusion_found = False
        self.phenofiles = []
        self.tarball_prefixes = []
        self.bgen_dict = {}

        self.default_options = {'transcript_index', 'bed_file', 'fam_file', 'bim_file', 'low_MAC_list', 'sparse_grm',
                                'sparse_grm_sample', 'phenofile', 'base_covariates', 'covarfile',
                                'association_tarballs', 'bgen_index', 'inclusion_list', 'exclusion_list'}

        self._check_required_options(required_options)

        # Work our way through all the resources we need
        self._ingest_docker_file()
        self._ingest_transcript_index(parsed_options['transcript_index'])
        self._ingest_genetic_data(parsed_options['bed_file'],
                                  parsed_options['fam_file'],
                                  parsed_options['bim_file'],
                                  parsed_options['low_MAC_list'],
                                  parsed_options['sparse_grm'],
                                  parsed_options['sparse_grm_sample'])
        self._ingest_phenofile(parsed_options['phenofile'])
        self._ingest_covariates(parsed_options['base_covariates'], parsed_options['covarfile'])
        self._ingest_tarballs(parsed_options['association_tarballs'])
        self._ingest_bgen(parsed_options['bgen_index'])
        self._define_exclusion_lists(parsed_options['inclusion_list'],
                                     parsed_options['exclusion_list'])

    def _check_required_options(self, required_options: set):

        for option in required_options:
            if option not in self._parsed_options:
                raise dxpy.AppError(f'Option {option} not found in input_yaml. Please re-run with this option!')

    # Bring our docker image into our environment so that we can run commands we need:
    @staticmethod
    def _ingest_docker_file() -> None:
        cmd = "docker pull egardner413/mrcepid-associationtesting:latest"
        run_cmd(cmd)

    # Get transcripts for easy annotation:
    @staticmethod
    def _ingest_transcript_index(transcript_index: dict) -> None:
        dxpy.download_dxfile(dxpy.DXFile(transcript_index).get_id(), 'transcripts.tsv.gz')

    @staticmethod
    def _ingest_genetic_data(bed_file: dict, fam_file: dict, bim_file: dict, low_MAC_list: dict,
                              sparse_grm: dict, sparse_grm_sample: dict) -> None:
        # Now grab all genetic data that I have in the folder /project_resources/genetics/
        os.mkdir("genetics/")  # This is for legacy reasons to make sure all tests work...
        dxpy.download_dxfile(dxpy.DXFile(bed_file).get_id(), 'genetics/UKBB_470K_Autosomes_QCd.bed')
        dxpy.download_dxfile(dxpy.DXFile(bim_file).get_id(), 'genetics/UKBB_470K_Autosomes_QCd.bim')
        dxpy.download_dxfile(dxpy.DXFile(fam_file).get_id(), 'genetics/UKBB_470K_Autosomes_QCd.fam')
        dxpy.download_dxfile(dxpy.DXFile(low_MAC_list).get_id(), 'genetics/UKBB_470K_Autosomes_QCd.low_MAC.snplist')
        # This is the sparse matrix
        dxpy.download_dxfile(dxpy.DXFile(sparse_grm).get_id(),
                             'genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx')
        dxpy.download_dxfile(dxpy.DXFile(sparse_grm_sample).get_id(),
                             'genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt')

    def _ingest_phenofile(self, phenofile: list) -> None:

        # There can be multiple phenotype files in phewas mode, so we iterate through the list of them here
        for pheno in phenofile:
            curr_pheno = dxpy.DXFile(pheno)
            curr_pheno_name = curr_pheno.describe()['name']

            dxpy.download_dxfile(curr_pheno.get_id(), curr_pheno_name)
            self.phenofiles.append(curr_pheno_name)

    # Get covariate/phenotype data:
    def _ingest_covariates(self, base_covariates: dict, covarfile: dict) -> None:
        dxpy.download_dxfile(dxpy.DXFile(base_covariates).get_id(), 'base_covariates.covariates')

        # Check if additional covariates were provided:
        if covarfile is not None:
            covarfile = dxpy.DXFile(covarfile)
            dxpy.download_dxfile(covarfile, 'additional_covariates.covariates')
            self.additional_covariates_found = True

    # Need to grab the tarball file for associations...
    # This was generated by the applet mrcepid-collapsevariants
    # Ingest the list file into this AWS instance
    def _ingest_tarballs(self, association_tarballs: dict) -> None:

        dxtarballs = dxpy.DXFile(association_tarballs)
        self.is_snp_tar = False
        self.is_gene_tar = False
        if '.tar.gz' in dxtarballs.describe()['name']:
            # likely to be a single tarball, download, check, and extract:
            tarball_name = dxtarballs.describe()['name']
            dxpy.download_dxfile(dxtarballs, tarball_name)
            if tarfile.is_tarfile(tarball_name):
                tarball_prefix = tarball_name.replace(".tar.gz", "")
                self.tarball_prefixes.append(tarball_prefix)
                tar = tarfile.open(tarball_name, "r:gz")
                tar.extractall()
                if exists(tarball_prefix + ".SNP.BOLT.bgen"):
                    self.is_snp_tar = True
                elif exists(tarball_prefix + ".GENE.BOLT.bgen"):
                    self.is_gene_tar = True
            else:
                raise dxpy.AppError("Provided association tarball (" + dxtarballs.describe()['id'] + ") is not a tar.gz file")
        else:
            # Likely to be a list of tarballs, download and extract...
            dxpy.download_dxfile(dxtarballs, "tarball_list.txt")
            with open("tarball_list.txt", "r") as tarball_reader:
                for association_tarball in tarball_reader:
                    association_tarball = association_tarball.rstrip()
                    tarball = dxpy.DXFile(association_tarball)
                    tarball_name = tarball.describe()['name']
                    dxpy.download_dxfile(tarball, tarball_name)

                    # Need to get the prefix on the tarball to access resources within:
                    # All files within SHOULD have the same prefix as this file
                    tarball_prefix = tarball_name.rstrip('.tar.gz')
                    self.tarball_prefixes.append(tarball_prefix)
                    tar = tarfile.open(tarball_name, "r:gz")
                    tar.extractall()
                    if exists(tarball_prefix + ".SNP.BOLT.bgen"):
                        raise dxpy.AppError("Cannot run masks from a SNP list (" + dxtarballs.describe()['id'] + ") when running tarballs as batch...")
                    elif exists(tarball_prefix + ".GENE.BOLT.bgen"):
                        raise dxpy.AppError("Cannot run masks from a GENE list (" + dxtarballs.describe()['id'] + ") when running tarballs as batch...")

    # Grab the entire WES variant data in bgen format
    def _ingest_bgen(self, bgen_index: dict) -> None:

        # Ingest the INDEX of bgen files:
        bgen_index = dxpy.DXFile(bgen_index)
        dxpy.download_dxfile(bgen_index.get_id(), "bgen_locs.tsv")
        # and load it into a dict:
        os.mkdir("filtered_bgen/")  # For downloading later...
        bgen_index_csv = csv.DictReader(open("bgen_locs.tsv", "r"), delimiter="\t")
        for line in bgen_index_csv:
            self.bgen_dict[line['chrom']] = {'index': line['bgen_index_dxid'],
                                             'sample': line['sample_dxid'],
                                             'bgen': line['bgen_dxid'],
                                             'vep': line['vep_dxid']}

    # Get inclusion/exclusion sample lists
    def _define_exclusion_lists(self, inclusion_list, exclusion_list) -> None:
        if inclusion_list is not None:
            inclusion_list = dxpy.DXFile(inclusion_list)
            dxpy.download_dxfile(inclusion_list, 'INCLUSION.lst')
            self.inclusion_found = True
        if exclusion_list is not None:
            exclusion_list = dxpy.DXFile(exclusion_list)
            dxpy.download_dxfile(exclusion_list, 'EXCLUSION.lst')
            self.exclusion_found = True
