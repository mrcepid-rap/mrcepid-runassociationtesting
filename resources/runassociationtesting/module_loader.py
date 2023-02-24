import re
import dxpy
import argparse
import importlib

from importlib import util
from typing import List, Type, Optional
from abc import ABC, abstractmethod

from runassociationtesting.association_pack import AssociationPack, ProgramArgs


class ModuleLoader(ABC):

    def __init__(self, output_prefix: str, input_args: str):
        self._outputs = []
        self.output_prefix = output_prefix

        self._input_args = input_args

        # Load the top-level required options
        self._parser = argparse.ArgumentParser()
        self._load_general_options()

        # Set all possible required options as a union of top-level options required by all tools and options specific
        # to a given tool. Then generate a union of these actual options to load.
        self._load_module_options()
        self.parsed_options = self._parse_options()

        # Download required files and process covariates (this will call the subclass _ingest_data,
        # NOT this abstractclass' _ingest_data)
        self.association_pack = self._ingest_data(self.parsed_options)

    def get_outputs(self) -> List[str]:
        return self._outputs

    def set_outputs(self, outputs: List[str]):
        self._outputs = outputs

    # A class that defines a 'dxfile' type for argparse. Allows for a 'None' input default for optional files
    @staticmethod
    def dxfile_input(input_str: str) -> Optional[dxpy.DXFile]:
        try:
            if input_str == 'None':
                return None
            else:
                dxfile = dxpy.DXFile(dxid=input_str)
                dxfile.describe()
                return dxfile
        except dxpy.exceptions.DXError:  # This just checks if the format of the input is correct
            raise TypeError(f'The input for parameter – {input_str} – '
                            f'does not look like a valid DNANexus file ID.')
        except dxpy.exceptions.ResourceNotFound:
            raise TypeError(f'The input for parameter – {input_str} – '
                            f'does not exist on the DNANexus platform.')

    @staticmethod
    def comma_str(input_str: str) -> List[str]:
        input_list = input_str.split(',')
        return input_list

    # Since I don't use the standard argparse (because I want to ensure options are specific to submodules)
    # I have written a separate 'splitter' to ensure options are properly parsed into a list properly
    @staticmethod
    def _split_options(input_args: str) -> List[str]:

        arg_indicies = [-1]
        last_index = -99
        for i in range(0, len(input_args)):
            curr_index = input_args.find(' -', i)
            if (curr_index - 1) != last_index and curr_index != last_index and curr_index != -1:
                arg_indicies.append(curr_index)
                last_index = curr_index

        split_args = []
        for i in range(0, len(arg_indicies)):
            if i == len(arg_indicies) - 1:
                start = arg_indicies[i] + 1
                end = len(input_args)
            else:
                start = arg_indicies[i] + 1
                end = arg_indicies[i+1]
            split_args.append(input_args[start:end])

        parsed_args = []
        for arg in split_args:

            opt_search = re.match('^(-{1,2}[\\w\\d-]+)\\s*([\\S\\s]+)?', arg)

            if opt_search:
                parsed_args.append(opt_search.group(1))
                if opt_search.group(2):
                    group = opt_search.group(2)
                    # Strip any leading/lagging quotes from anything we parse:
                    str_len = len(group)
                    for match in re.finditer('[\'\"]', group):
                        if match.start() == 0:
                            group = group[1:len(group)]
                        elif match.end() == str_len:
                            group = group[0:len(group)-1]
                    parsed_args.append(group)

            else:
                raise ValueError(f'Incorrectly formatted argument string {arg}')

            # last match might not have an end group and len() in the 'while' statement cannot have 'None' as input type
            input_args = '' if input_args is None else input_args
        return parsed_args

    def _load_general_options(self) -> None:

        example_dxfile = 'file-123...'

        self._parser.add_argument('--phenofile',
                                  help="Path/hash to phenotype file(s). If providing multiple files, provide "
                                       "them space separated, one after another.",
                                  type=self.dxfile_input, dest='phenofile', required=True, nargs='+',
                                  metavar=example_dxfile)
        self._parser.add_argument('--phenoname',
                                  help="A single phenotype name to run association tests for. This MUST EXACTLY match "
                                       "the header of the provided phenofile to function properly. Only required if "
                                       "providing a pheno file with multiple phenotype columns (other than IID/FID).",
                                  type=str, dest='phenoname', required=False, default=None)
        self._parser.add_argument('--covarfile',
                                  help="Path/hash to additional covariates file.",
                                  type=self.dxfile_input, dest='covarfile', required=False, default=None,
                                  metavar=example_dxfile)
        self._parser.add_argument('--categorical_covariates',
                                  help="A space-separated list (e.g. covar1,covar2,covar3) of categorical "
                                       "(e.g. WES batch) in <covarfile>. Names MUST match column header.",
                                  type=str, dest='categorical_covariates', required=False, default=None,
                                  nargs='*', metavar=('CAT_COVAR1', 'CAT_COVAR2'))
        self._parser.add_argument('--quantitative_covariates',
                                  help="A space-separated list (e.g. covar1,covar2,covar3) of categorical (e.g. PC1) "
                                       "in <covarfile>. Names MUST match column header.",
                                  type=str, dest='quantitative_covariates', required=False, default=None,
                                  nargs='*', metavar=('QUANT_COVAR1', 'QUANT_COVAR2'))
        self._parser.add_argument('--is_binary',
                                  help="Is the phenotype binary?",
                                  dest='is_binary', action='store_true')
        self._parser.add_argument('--sex',
                                  help="Run only one sex or both sexes (0 = female, 1 = male, 2 = both) [2]?",
                                  type=int, dest='sex', required=False, default=2, choices=[0, 1, 2])
        self._parser.add_argument('--exclusion_list',
                                  help="File of individual IDs to exclude from analyses. Do not use the command to "
                                       "exclude individuals by sex!",
                                  type=self.dxfile_input, dest='exclusion_list', required=False, default=None,
                                  metavar=example_dxfile)
        self._parser.add_argument('--inclusion_list',
                                  help="File of individual IDs to include in analyses.",
                                  type=self.dxfile_input, dest='inclusion_list', required=False, default=None,
                                  metavar=example_dxfile)
        self._parser.add_argument('--transcript_index',
                                  help="transcripts for annotation of association stats",
                                  type=self.dxfile_input, dest='transcript_index', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--base_covariates',
                                  help="A file of standard covariates to include in the model. This file follows a "
                                       "very specific format outlined in the README and should not be changed unless "
                                       "absolutely necessary",
                                  type=self.dxfile_input, dest='base_covariates', required=True,
                                  metavar=example_dxfile)

    @abstractmethod
    def start_module(self) -> None:
        pass

    @abstractmethod
    def _load_module_options(self) -> None:
        pass

    @abstractmethod
    def _parse_options(self) -> ProgramArgs:
        pass

    @abstractmethod
    def _ingest_data(self, parsed_options: ProgramArgs) -> AssociationPack:
        pass


# This method takes a possible 'mode' and tries to import it into the current projects namespace.
# This is to enable 'modules' that are not part of the standard distribution to be run without being part of the
# project at runtime
def conditional_import(package: str) -> Type[ModuleLoader]:

    # Return a reference to the module itself
    loader = util.find_spec(package)
    if loader is not None:
        loader_package = importlib.import_module(package)
        return loader_package.LoadModule
    else:
        raise ModuleNotFoundError
