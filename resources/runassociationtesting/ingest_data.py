import csv
import os
import tarfile

import dxpy

from association_resources import *


# This class is slightly different that in other applets I have designed; it handles ALL inputs rather
# than just external dependencies
class IngestData:

    def __init__(self, association_tarballs: dict, phenofile: dict, covarfile: dict, inclusion_list: dict,
                 exclusion_list: dict, bgen_index: dict, transcript_index: dict, base_covariates: dict,
                 bed_file: dict, fam_file: dict, bim_file: dict, low_MAC_list: dict,
                 sparse_grm: dict, sparse_grm_sample: dict):

        # Set reasonable defaults
        self.additional_covariates_found = False
        self.inclusion_found = False
        self.exclusion_found = False
        self.tarball_prefixes = []
        self.bgen_dict = {}

        self._ingest_docker_file()
        self._ingest_transcript_index(transcript_index)
        self._ingest_genetic_data(bed_file, fam_file, bim_file, low_MAC_list, sparse_grm, sparse_grm_sample)
        self._ingest_covariates(base_covariates, phenofile, covarfile)
        self._ingest_tarballs(association_tarballs)
        self._ingest_bgen(bgen_index)
        self._define_exclusion_lists(inclusion_list, exclusion_list)

    # Bring our docker image into our environment so that we can run commands we need:
    @staticmethod
    def _ingest_docker_file() -> None:
        cmd = "docker pull egardner413/mrcepid-associationtesting:latest"
        run_cmd(cmd)

    # Get transcripts for easy annotation:
    @staticmethod
    def _ingest_transcript_index(transcript_index: dict):
        dxpy.download_dxfile(dxpy.DXFile(transcript_index).get_id(), 'transcripts.tsv.gz')

    @staticmethod
    def _ingest_genetic_data(bed_file: dict, fam_file: dict, bim_file: dict, low_MAC_list: dict,
                              sparse_grm: dict, sparse_grm_sample: dict):
        # Now grab all genetic data that I have in the folder /project_resources/genetics/
        os.mkdir("genetics/")  # This is for legacy reasons to make sure all tests work...
        dxpy.download_dxfile(dxpy.DXFile(bed_file).get_id(), 'genetics/UKBB_450K_Autosomes_QCd.bed')
        dxpy.download_dxfile(dxpy.DXFile(bim_file).get_id(), 'genetics/UKBB_450K_Autosomes_QCd.bim')
        dxpy.download_dxfile(dxpy.DXFile(fam_file).get_id(), 'genetics/UKBB_450K_Autosomes_QCd.fam')
        dxpy.download_dxfile(dxpy.DXFile(low_MAC_list).get_id(), 'genetics/UKBB_450K_Autosomes_QCd.low_MAC.snplist')
        # This is the sparse matrix
        dxpy.download_dxfile(dxpy.DXFile(sparse_grm).get_id(),
                             'genetics/fixed_rel.sorted.mtx')
        dxpy.download_dxfile(dxpy.DXFile(sparse_grm_sample).get_id(),
                             'genetics/fixed_rel.sorted.mtx.sampleIDs.txt')

    # Get covariate/phenotype data:
    def _ingest_covariates(self, base_covariates: dict, phenofile: dict, covarfile: dict):
        dxpy.download_dxfile(dxpy.DXFile(base_covariates).get_id(), 'base_covariates.covariates')
        dxpy.download_dxfile(dxpy.DXFile(phenofile).get_id(), 'model_phenotypes.pheno')
        # Check if additional covariates were provided:
        if covarfile is not None:
            covarfile = dxpy.DXFile(covarfile)
            dxpy.download_dxfile(covarfile, 'additional_covariates.covariates')
            self.additional_covariates_found = True

    # Need to grab the tarball file for associations...
    # This was generated by the applet mrcepid-collapsevariants
    # Ingest the list file into this AWS instance
    def _ingest_tarballs(self, association_tarballs: dict):

        dxtarballs = dxpy.DXFile(association_tarballs)
        if '.tar.gz' in dxtarballs.describe()['name']:
            # likely to be a single tarball, download, check, and extract:
            tarball_name = dxtarballs.describe()['name']
            dxpy.download_dxfile(dxtarballs, tarball_name)
            if tarfile.is_tarfile(tarball_name):
                tarball_prefix = tarball_name.rstrip('.tar.gz')
                self.tarball_prefixes.append(tarball_prefix)
                tar = tarfile.open(tarball_name, "r:gz")
                tar.extractall()
            else:
                dxpy.AppError("Provided association tarball (" + dxtarballs.describe()['id'] + ") is not a tar.gz file")
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

    # Grab the entire WES variant data in bgen format
    def _ingest_bgen(self, bgen_index: dict):

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
    def _define_exclusion_lists(self, inclusion_list, exclusion_list):
        if inclusion_list is not None:
            inclusion_list = dxpy.DXFile(inclusion_list)
            dxpy.download_dxfile(inclusion_list, 'INCLUSION.lst')
            self.inclusion_found = True
        if exclusion_list is not None:
            exclusion_list = dxpy.DXFile(exclusion_list)
            dxpy.download_dxfile(exclusion_list, 'EXCLUSION.lst')
            self.exclusion_found = True
