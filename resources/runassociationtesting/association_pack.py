import dxpy

from abc import ABC
from dataclasses import dataclass
from typing import List


@dataclass
class ProgramArgs(ABC):
    """A @dataclass that stores information on arguments that are passed to optparse

    This class is for ease of programming and has no actual functionality for any of the processing occurring in this
    workflow. It allows for coders to know possible input parameters and to get proper typehints when adding new
    functionality. For more information see either _load_module_options() in ModuleLoader (module_loader.py) or the
    implemented version of this method in individual modules subclasses that implement ModuleLoader.
    """

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
    """An interface that stores information necessary for genetic association.

    All modules that can be run through RunAssociationTesting MUST implement this interface to store information
    # required for genetic association. The parameters stored in this interface are those required for all such modules
    (i.e., these parameters are required no matter what analysis is being performed) and are processed and implemented
    in the Additional subclasses can implement this interface to add additional parameters.

    :param is_binary: Is the phenotype binary?
    :param sex: Sex to perform sex-stratified analysis for. Possible values are 0 = Female, 1 = Male, 2 = Both.
    :param threads: Number of threads available to the current machine running this job.
    :param pheno_names: str names of phenotypes requested to be run for this analysis.
    :param found_quantitative_covariates: str names of any found quantitative covariates.
    :param found_categorical_covariates: str names of any found categorical covariates.
    """

    def __init__(self, is_binary: bool, sex: int, threads: int, pheno_names: List[str],
                 found_quantitative_covariates: List[str], found_categorical_covariates: List[str]):

        self.is_binary = is_binary
        self.sex = sex
        self.threads = threads
        self.pheno_names = pheno_names
        self.found_quantitative_covariates = found_quantitative_covariates
        self.found_categorical_covariates = found_categorical_covariates
