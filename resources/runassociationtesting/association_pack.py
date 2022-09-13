import dxpy

from abc import ABC
from dataclasses import dataclass
from typing import List


@dataclass
class ProgramArgs(ABC):
    phenofile: List[dxpy.DXFile]
    phenoname: str
    covarfile: dxpy.DXFile
    categorical_covariates: List[str]
    quantitative_covariates: List[str]
    is_binary: bool
    sex: int
    exclusion_list: dxpy.DXFile
    inclusion_list: dxpy.DXFile
    transcript_index: dxpy.DXFile
    base_covariates: dxpy.DXFile


class AssociationPack(ABC):

    def __init__(self, pheno_files: List[str],
                 inclusion_found: bool, exclusion_found: bool, additional_covariates_found: bool,
                 is_binary: bool, sex: int, threads: int, pheno_names: List[str],
                 found_quantitative_covariates: List[str], found_categorical_covariates: List[str]):

        self.pheno_files = pheno_files
        self.exclusion_found = exclusion_found
        self.inclusion_found = inclusion_found
        self.additional_covariates_found = additional_covariates_found
        self.is_binary = is_binary
        self.sex = sex
        self.threads = threads
        self.pheno_names = pheno_names
        self.found_quantitative_covariates = found_quantitative_covariates
        self.found_categorical_covariates = found_categorical_covariates
