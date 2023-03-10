#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# Prior to using this script PLEASE make sure test_data either has the following files:
# base_covar.scrambled.tsv          covar.tsv               exclude.txt             include.txt
# pheno_covar.scrambled.tsv         pheno_covar.txt         transcripts.tsv.gz      base_covar.tsv
# covar.wrong_header.tsv            expected_results.tsv    pheno.tsv               pheno_covar.tsv
# pheno_covar.wrong_header.tsv      pheno.first.tsv         pheno.second.tsv
#
# OR that the script `generate_test_data.R` has been run to generate them!

import os
import sys
from pathlib import Path

import pytest
import pandas as pd

from typing import Any, List, Dict

# DO NOT move the sys.path.append() calls below - they are required to be able to import the runassociationtesting
# classes that follow on the DNANexus platform for proper testing.
sys.path.append('/')
sys.path.append('/runassociationtesting/')
from runassociationtesting.module_loader import conditional_import, ModuleLoader



test_folder = os.getenv('TEST_DIR')


@pytest.mark.parametrize(
    argnames=['input_str', 'expected_filename'],
    argvalues=zip(['file-Fx2x270Jx0j17zkb3kbBf6q2', f'{test_folder}/base_covar.tsv'],
                  ['hs38DH.fa.gz', 'base_covar.tsv'])
)
def test_dxfile_input(input_str: str, expected_filename: Any):
    """Tests that the ModuleLoader.dxfile_input can take (and find) a file on the DNoneexus platform that is either a
    file-ID or absolute path.

    :param input_str: Either a DNoneexus file-ID (e.g, file-12345...) or an absolute path.
    """
    assert ModuleLoader.dxfile_input(input_str).describe(fields={'name': True})['name'] == expected_filename


@pytest.mark.parametrize(
    argnames=['input_str', 'expected_exception'],
    argvalues=zip(['file-1234567890AbCdeFgHiJkLmN', 'file-12345', '/A random path/to/hell.txt', 'foo'],
                  [TypeError, TypeError, FileNotFoundError, TypeError])
)
def test_dxfile_exceptions(input_str: str, expected_exception: Exception):
    """Make sure that the correct exceptions are thrown when incorrect test_data is given.

    The inputs should test the following errors (in order):

    1. A DNoneexus file ID that doesn't actually point to a real file

    2. A string that looks like a DNoneexus file ID but actually isn't

    3. A path that doesn't point to anything

    4. A random standalone string

    :param input_str:
    :param expected_exception:
    """

    with pytest.raises(expected_exception):
        ModuleLoader.dxfile_input(input_str)


@pytest.mark.parametrize(
    argnames=['input_str', 'expected_list'],
    argvalues=zip(['1,2,3,4', 'param1,param2,param3'],
                  [['1', '2', '3', '4'], ['param1', 'param2', 'param3']])
)
def test_comma_str(input_str: str, expected_list: List[str]):
    """Test that the comma_str function in module_loader splits properly

    :param input_str: A comma-seperated string
    :param expected_list: A list representing the split string
    """

    assert ModuleLoader.comma_str(input_str) == expected_list


def test_split_arguments():
    """Test the custom argument parser that is built into RunAssociationTesting.

    I have tried to test as many possible argument parameters types (e.g., boolean flags, quote-surrounded,
    space-delimited lists etc.) as possible here, but may have missed a few. I will continue to check what possible
    arguments could be entered. Some should never occur, but I test here for future proofing (e.g., comma_list_opt is
    not used anywhere in the current codebase and crazy_arg is unlikely to ever happen)
    """
    input_args = '--boolean_opt --dx_opt file-1234567890AbCdeFgHiJkLmN --str_opt foo ' \
                 '--path_single_quote \'/Path to a file/file.txt\' --path_double_quote "/Path to a file/file.txt" ' \
                 '--path_no_quote /Path_to_a_file/file.txt ' \
                 '--list_opt arg1 arg2 arg3 arg4 --comma_list_opt arg1,arg2,arg3,arg4 ' \
                 '--dash_opt arg-with-dashes -a aShortArg --crazy_arg A$arg.with_#crÅzy_€har€t€rs --final_boolean_opt'
    output_arg_list = ['--boolean_opt', '--dx_opt', 'file-1234567890AbCdeFgHiJkLmN', '--str_opt', 'foo',
                       '--path_single_quote', '/Path to a file/file.txt', '--path_double_quote',
                       '/Path to a file/file.txt', '--path_no_quote', '/Path_to_a_file/file.txt',
                       '--list_opt', 'arg1', 'arg2', 'arg3', 'arg4', '--comma_list_opt', 'arg1,arg2,arg3,arg4',
                       '--dash_opt', 'arg-with-dashes', '-a', 'aShortArg', '--crazy_arg', 'A$arg.with_#crÅzy_€har€t€rs',
                       '--final_boolean_opt']

    assert ModuleLoader._split_options(input_args) == output_arg_list


def test_import():
    """Test the basics of the conditional_import function in ingest_data.py

    I have created a test repo (mrcepid-test_loader) that contains a minimal implementation of the classes required
    to start a module in runassociationtesting. This test will ask to find that module and return the class (
    ModuleLoader) that allows for startup happen. Individual modules should have a test to ensure the module-specific
    startup functionality works; module-specific startup is NOT tested here. However, this will test if the correct
    module structure is implemented through the parsing of imports.
    """

    module_loader = conditional_import(f'test_loader.loader')
    assert issubclass(module_loader, ModuleLoader)


# To regenerate this list, please use the script `generate_load_parameters.py`. This script will print a text-based
# format of this table that can be copy-pasted into this file. Do to how I get data onto DNANexus, this is the best I
# could do
load_parameters = [
    {'test_name': 'base mode', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 890, 'n_expected_columns': 47}},
    {'test_name': 'base mode scrambled base columns', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.scrambled.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 890, 'n_expected_columns': 47}},
    {'test_name': 'base mode scrambled pheno columns', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.scrambled.tsv', 'phenoname': 'quantPheno1', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 890, 'n_expected_columns': 47}},
    {'test_name': 'no sex specified', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': None, 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 890, 'n_expected_columns': 47}},
    {'test_name': 'wrong pheno name', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno4', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': 'phenoname was not found', 'n_expected_samples': 0, 'n_expected_columns': 0}},
    {'test_name': 'no FID pheno', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.wrong_header.tsv', 'phenoname': 'quantPheno1', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': 'Pheno file does not contain', 'n_expected_samples': 0, 'n_expected_columns': 0}},
    {'test_name': 'no named pheno', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno.tsv', 'phenoname': None, 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 900, 'n_expected_columns': 49}},
    {'test_name': 'binary pheno', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'catPheno3', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': True, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 890, 'n_expected_columns': 47}},
    {'test_name': 'named covars', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': 'pheno_covar.tsv', 'categorical_covariates': 'catCovar3', 'quantitative_covariates': 'quantCovar1 quantCovar2', 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 860, 'n_expected_columns': 50}},
    {'test_name': 'space-delimited pheno_covar', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.txt', 'phenoname': 'quantPheno1', 'covarfile': 'pheno_covar.txt', 'categorical_covariates': 'catCovar3', 'quantitative_covariates': 'quantCovar1 quantCovar2', 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 860, 'n_expected_columns': 50}},
    {'test_name': 'single covar', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': 'pheno_covar.tsv', 'categorical_covariates': None, 'quantitative_covariates': 'quantCovar1', 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 880, 'n_expected_columns': 48}},
    {'test_name': 'separate pheno and covars', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno.tsv', 'phenoname': 'quantPheno1', 'covarfile': 'covar.tsv', 'categorical_covariates': 'catCovar3', 'quantitative_covariates': 'quantCovar1 quantCovar2', 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 860, 'n_expected_columns': 50}},
    {'test_name': 'multiple phenos', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno.first.tsv pheno.second.tsv', 'phenoname': None, 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 900, 'n_expected_columns': 49}},
    {'test_name': 'no FID covar', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno.tsv', 'phenoname': 'quantPheno1', 'covarfile': 'covar.wrong_header.tsv', 'categorical_covariates': 'catCovar3', 'quantitative_covariates': 'quantCovar1 quantCovar2', 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': 'FID & IID column not found', 'n_expected_samples': 0, 'n_expected_columns': 0}},
    {'test_name': 'wrong covar name', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': 'pheno_covar.tsv', 'categorical_covariates': None, 'quantitative_covariates': 'quantCovar4', 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': 'Additional covariate file provided', 'n_expected_samples': 0, 'n_expected_columns': 0}},
    {'test_name': 'females', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '0', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 495, 'n_expected_columns': 47}},
    {'test_name': 'males', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '1', 'exclusion_list': None, 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 395, 'n_expected_columns': 47}},
    {'test_name': 'exclusion', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '2', 'exclusion_list': 'exclude.txt', 'inclusion_list': None}, 'outputs': {'error_type': None, 'n_expected_samples': 830, 'n_expected_columns': 47}},
    {'test_name': 'inclusion', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '2', 'exclusion_list': None, 'inclusion_list': 'include.txt'}, 'outputs': {'error_type': None, 'n_expected_samples': 700, 'n_expected_columns': 47}},
    {'test_name': 'exclusion inclusion', 'parameters': {'transcript_index': 'transcripts.tsv.gz', 'base_covariates': 'base_covar.tsv', 'phenofile': 'pheno_covar.tsv', 'phenoname': 'quantPheno1', 'covarfile': None, 'categorical_covariates': None, 'quantitative_covariates': None, 'is_binary': False, 'sex': '2', 'exclusion_list': 'exclude.txt', 'inclusion_list': 'include.txt'}, 'outputs': {'error_type': None, 'n_expected_samples': 690, 'n_expected_columns': 47}}
]


@pytest.mark.parametrize(
    argnames=['test'],
    argvalues=zip(load_parameters)
)
def test_load_general_options(test: Dict):
    """Test loading of options and test_data processing via the test_loader module

    Since all classes included in the primary mrcepid-runassociationtesting package are interfaces, I have created a
    'dummy' module (test_loader) that we can use to test. `test_loader` has no additional parameters (i.e.,
    argparse) nor does it have additional files / test_data that it ingests BEYOND that already specified in the
    interfaces within mrcepid-runassociationtesting.

    We are testing the following parameters – Need to test each one of these in turn (and some multiple times for
    different processing parameters).

    Required arguments:
    --phenofile, --transcript_index, --base_covariates

    Optional arguments:
    --phenoname, --covarfile, --categorical_covariates, --quantitative_covariates, --is_binary, --sex,
    --exclusion_list, --inclusion_list

    We actually test several parts of the module loading process here:

    1. The number of columns and samples in the combined phenotype/covariate file (always phenotypes_covariates.formatted.txt) is correct

    2. The number of samples in the SAMPLES_Include.txt file is correct

    3. The :func:`start_module` is able to run after startup. To test this, test_loader creates a blank test file with a name derived from the output prefix: '{output_prefix}.start_worked.txt'

    4. The correct information for this test is written to the AssociationPack class

    :param test: Test dictionary. Keys of 'test_name', 'parameters', and 'outputs'.
    """

    test_name = test['test_name']
    error_type = test['outputs']['error_type']
    n_expected_samples = test['outputs']['n_expected_samples']
    n_expected_columns = test['outputs']['n_expected_columns']

    # Parameters that point to a file in test_folder
    file_params = ['transcript_index', 'base_covariates', 'phenofile', 'covarfile', 'inclusion_list', 'exclusion_list']
    input_args = []
    for param, option in test['parameters'].items():
        if option is not None:
            if param == 'is_binary':
                if option is True:
                    input_args.append('--is_binary')
            elif param in file_params:
                if param == 'phenofile':
                    pheno_formatted = ' '.join([f'{test_folder}/{file}' for file in option.split()])
                    input_args.append(f'--{param} {pheno_formatted}')
                else:
                    input_args.append(f'--{param} {test_folder}/{option}')
            else:
                input_args.append(f'--{param} {option}')

    input_args = ' '.join(input_args)
    output_prefix = test_name.replace(' ', '_')

    if error_type is None:
        # 1. Run module loader
        module_loader = conditional_import(f'test_loader.loader')
        loaded_module = module_loader(output_prefix, input_args)
        association_pack = loaded_module.association_pack
        loaded_module.start_module()

        # 2. Check output files
        test_file = pd.read_csv('phenotypes_covariates.formatted.txt', sep='\\s+')
        n_col = len(test_file.columns)
        n_row = len(test_file)
        n_samples_file = 0
        with Path('SAMPLES_Include.txt').open('r') as sample_file:
            for _ in sample_file:
                n_samples_file += 1

        assert n_col == n_expected_columns, f'Actual column names: {test_file.columns}'
        assert n_row == n_expected_samples
        assert n_samples_file == n_expected_samples
        assert Path(f'{output_prefix}.start_worked.txt').exists()

        # 3. Check loaded association pack information:
        # Set lists of information:
        if test['parameters']['phenoname'] is None:
            test_pheno_names = ['quantPheno1', 'quantPheno2', 'catPheno3']
        else:
            test_pheno_names = [test['parameters']['phenoname']]

        if test['parameters']['categorical_covariates'] is None:
            test_cat_names = []
        else:
            test_cat_names = test['parameters']['categorical_covariates'].split()

        if test['parameters']['quantitative_covariates'] is None:
            test_quant_names = []
        else:
            test_quant_names = test['parameters']['quantitative_covariates'].split()

        assert association_pack.is_binary == test['parameters']['is_binary']
        assert association_pack.sex == (int(test['parameters']['sex']) if test['parameters']['sex'] is not None else 2)
        assert association_pack.threads == os.cpu_count()
        assert sorted(association_pack.pheno_names) == sorted(test_pheno_names)
        assert sorted(association_pack.found_categorical_covariates) == sorted(test_cat_names)
        assert sorted(association_pack.found_quantitative_covariates) == sorted(test_quant_names)

    else:
        with pytest.raises(Exception, match=test['outputs']['error_type']):
            module_loader = conditional_import(f'test_loader.loader')
            loaded_module = module_loader(output_prefix, input_args)
            loaded_module.start_module()



