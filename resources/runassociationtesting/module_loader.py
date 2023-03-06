import re
import dxpy
import argparse

from importlib import util, import_module
from typing import List, Type, Optional
from abc import ABC, abstractmethod

from general_utilities.mrc_logger import MRCLogger
from runassociationtesting.association_pack import AssociationPack, ProgramArgs


class ModuleLoader(ABC):
    """An interface for loading individual modules and their parameters into the DNANexus environment.

    ModuleLoader is the primary entrypoint into individual modules. To ensure proper loading of individual modules, all
    modules **must**:

     1. Have a python file in the top-level package of their namespace named `loader.py` that **must**
     2. Implement a subclass of ModuleLoader that **must**
     3. be named LoadModule

    This is to ensure that :func:`conditional_import` can find the appropriate packages to load the module and that
    :func:`main` in runassociationtesting can use the exact same name to run the module.

    This interface will always load a set of default options shared by all modules, and process them with ArgParse. It
    will then force the module that implements this interface to process these
    arguments and load any individual arguments (via :func:`_` that it requests in its own implementation.

    :param output_prefix: A prefix to name the output tarball returned by this method.
    :param input_args: A string containing options that conform to optparse specification. These are processed
        downstream in the ModuleLoader class.
    """

    def __init__(self, output_prefix: str, input_args: str):

        # Initiate logger – This can also be used by a class which implements this Interface
        self._logger = MRCLogger(__name__).get_logger()

        # Set empty outputs array to hold anything a module may produce and returnable to the main() method calling this
        # class
        self._outputs = []
        self.output_prefix = output_prefix

        # Load the top-level required options
        self._input_args = input_args
        self._parser = argparse.ArgumentParser()
        self._load_general_options()

        # Set all possible required options as a union of top-level options required by all tools and options specific
        # to a given tool. Then generate a union of these actual options to load.
        self._load_module_options()
        self.parsed_options = self._parse_options()

        # Download required files and process covariates (this will call the version implemented by the module
        # subclass _ingest_data, NOT this abstractclass' _ingest_data)
        self.association_pack = self._ingest_data(self.parsed_options)

    def get_outputs(self) -> List[str]:
        """Getter for the list of output filepaths in str format"""

        return self._outputs

    def set_outputs(self, outputs: List[str]):
        """Setter for the list of output filepaths in str format"""

        self._outputs = outputs

    @staticmethod
    def dxfile_input(input_str: str) -> Optional[dxpy.DXFile]:
        """A method that defines a 'dxfile' type for argparse. Allows for a 'None' input default for optional files

        This method is to allow for a dxfile 'type' when defining arguments for argparse. This method is passed as an
        argument to :func:`argparse.ArgumentParser().add_argument()` method as type=self.dxfile_input.

        :param input_str: A DXFile ID in the form file-1234567890ABCDEFG
        :return: A 'nullable' dxpy.DXFile
        """
        try:
            if input_str == 'None':
                return None
            else:
                dxfile = dxpy.DXFile(dxid=input_str)
                dxfile.describe()  # This will trigger the Exceptions caught below if not actually a DXFile / not found
                return dxfile
        except dxpy.exceptions.DXError:  # This just checks if the format of the input is correct
            raise TypeError(f'The input for parameter – {input_str} – '
                            f'does not look like a valid DNANexus file ID.')
        except dxpy.exceptions.ResourceNotFound:
            raise TypeError(f'The input for parameter – {input_str} – '
                            f'does not exist on the DNANexus platform.')

    @staticmethod
    def comma_str(input_str: str) -> List[str]:
        """A method that defines a 'comma_str' type for argparse. Allows for comma-seperated strings of arguments rather
        than space-delimited.

        This method is passed as an argument to :func:`argparse.ArgumentParser().add_argument()` method as
        type=self.comma_str.

        :param input_str: A string putatively separated by commas
        :return: A list of strings representing the comma-split input_str
        """

        input_list = input_str.split(',')
        return input_list


    @staticmethod
    def _split_options(input_args: str) -> List[str]:
        """A custom argument parser that allows for the use of the association_pack.ProgramArgs class

        Since I don't use the standard argparse (because I want to ensure options are specific to submodules)
        I have written a separate 'splitter' to ensure options are parsed into a list properly. This allows the options
        to be passed to argparse.parse_args() as a list, which then returns a Tuple with named arguments to ProgramArgs
        and all of it's implementing classes. This in turn allows access of parsed args by name using autocomplete in
        pycharm like::

            parsed_options.option1

        :param input_args: The string of options provided to mrcepid-runassociationtesting -iinput_args in argparse
            'format'
        :return: A List of options and their arguments in a format readable by argparse.parse_args()
        """

        # First find all locations of ' -', which will tell us where all flags in the string are
        arg_indicies = [-1]
        last_index = -99
        for i in range(0, len(input_args)):
            curr_index = input_args.find(' -', i)
            if (curr_index - 1) != last_index and curr_index != last_index and curr_index != -1:
                arg_indicies.append(curr_index)
                last_index = curr_index

        # Then find the start of all flags and their arguments to generate coordinates within the string that we can
        # substring on:
        # --arg1 test --boolArg2 --arg3 fun on a bun --arg4 "/A bad file/path/foo.txt"
        # becomes:
        # ['--arg1 test', '--boolArg2', '--arg3 fun on a bun', '--arg4 "/A bad file/path/foo.txt"']
        split_args = []
        for i in range(0, len(arg_indicies)):
            if i == len(arg_indicies) - 1:
                start = arg_indicies[i] + 1
                end = len(input_args)
            else:
                start = arg_indicies[i] + 1
                end = arg_indicies[i+1]
            split_args.append(input_args[start:end])

        # And finally split the flag from the argument(s) potentially further splitting the argument by whitespace if
        # not surrounded by quotes
        parsed_args = []
        for arg in split_args:

            opt_search = re.match('^(-{1,2}[\\w\\d-]+)\\s*([\\S\\s]+)?', arg)

            if opt_search:
                parsed_args.append(opt_search.group(1))
                if opt_search.group(2):
                    group = opt_search.group(2)

                    # Do post-processing for quoted arguments:
                    # Strip any leading/lagging quotes from anything we parse:
                    str_len = len(group)
                    quote_match = 0
                    for match in re.finditer('[\'\"]', group):
                        if match.start() == 0:
                            quote_match += 1
                            group = group[1:len(group)]
                        elif match.end() == str_len:
                            quote_match += 1
                            group = group[0:len(group)-1]

                    # If quote_match = 2, then the parameter was wrapped in quotes and we don't strip by whitespace,
                    # add the argument directly as-is
                    if quote_match == 2:
                        parsed_args.append(group)

                    # If not 2, then we need to strip by whitespace and add each argument separately
                    else:
                        split_group = group.split()
                        for split in split_group:
                            parsed_args.append(split)

            else:
                raise ValueError(f'Incorrectly formatted argument string {arg}')

            # last match might not have an end group and len() in the 'while' statement cannot have 'None' as input type
            input_args = '' if input_args is None else input_args

        return parsed_args

    def _load_general_options(self) -> None:
        """A sort of DataClass that just adds arguments to argparse. The options here are defined for all
        implementing sub-classes."""

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
        """Triggers module-specific loading required for the functionality of that module.

        As an abstractmethod, all modules must call this method to ensure something other than optionparsing and data
        loading is performed. What this method actually does is completely up to the user.
        """
        pass

    @abstractmethod
    def _load_module_options(self) -> None:
        """Load module-specific options.

        This method is CALLED in the constructor of the interface, but its functionality is DEFINED in the module
        currently being 'loaded'. This method effectively works as an 'add-on' to the :func:`_load_general_options`
        function and loads module-specific options into the same argparse instance as this interface.
        """
        pass

    @abstractmethod
    def _parse_options(self) -> ProgramArgs:
        """Parse input_args into individual options that can be accessed as a named Tuple.

        This method is CALLED in the constructor of the interface, but its functionality is DEFINED in the module
        currently being 'loaded'; however, the code that is run in the sublclass implementation MUST follow the same
        rough syntax::

                return ProgramArgs(**vars(self._parser.parse_args(self._split_options(self._input_args))))

        Where ProgramArgs is the specific subclass that implements ProgramArgs for the current module.

        :return: A subclass of ProgramArgs containing parsed options
        """
        pass

    @abstractmethod
    def _ingest_data(self, parsed_options: ProgramArgs) -> AssociationPack:
        """Download, process, or store files and data required for module functionality

        This method is CALLED in the constructor of the interface, but its functionality is DEFINED in the module
        currently being 'loaded'.

        :param parsed_options: A subclass of ProgramArgs containing options parsed by argparse
        :return: A subclass of AssociationPack containing runtime parameters required by the module to function
        """
        pass


def conditional_import(package: str) -> Type[ModuleLoader]:
    """Import a requested module.

    This method takes a possible 'package' and tries to import it into the current projects namespace.
    This is to enable 'modules' that are not part of the standard distribution to be run without being part of the
    project at runtime. All functionality that runassociationtesting 'implements' has now been moved out of the primary
    runassociationtesting app/applet and placed in individual modules which use interfaces in the runassociationtesting
    package to interface with the DNANexus platform.

    This functionality is done via the 'execDepends' section of the DNAnexus dxapp.json that is used to build this app.
    Individual packages are included as json elements like::

        {
            "name": "burden",
            "package_manager": "git",
            "url":  "https://github.com/mrcepid-rap/mrcepid-runassociationtesting-burden.git",
            "build_commands": "pip3 install ."
        }

    When a DNANexus instance is requested, these packages are installed into the python3 namespace (by default using
    pip install .), and can then be imported ad-hoc with :func:`find_spec` + :func:`import_module` supplied by the
    importlib package. import_module returns a reference to the top-level instance of the module as created by
    `__init__` (e.g., `burden.__init__`). We then return a reference to the `LoadModule` class within this package (
    e.g., `burden.LoadModule` that **must** be implemented by a given module. This class can then be instantiated
    like normal as a reference in the :func:`main` of `mrcepid-runassociationtesting.py` (e.g. `LoadModule(
    output_prefix, input_args)`).

    For more information on this process and how to build a module of your own, see the developer README in this
    repository.

    :param package: The name of a python3 package installed via GitHub
    :return: An instance of a class that implements the interface ModuleLoader
    """

    loader = util.find_spec(package)
    if loader is not None:
        loader_package = import_module(package)
        return loader_package.LoadModule
    else:
        raise ModuleNotFoundError
