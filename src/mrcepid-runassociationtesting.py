#!/usr/bin/env python
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/
import os
import dxpy
import tarfile
import pkg_resources

from pathlib import Path

from general_utilities.association_resources import generate_linked_dx_file
from general_utilities.job_management.subjob_utility import check_subjob_decorator
from general_utilities.mrc_logger import MRCLogger
from general_utilities.job_management.command_executor import CommandExecutor
from general_utilities.import_utils.module_loader.module_loader import conditional_import

MRC_LOGGER = MRCLogger()
LOGGER = MRC_LOGGER.get_logger()

# Before we do anything else, we MUST check if subjobs are being run through this script so that we can properly load
# the required dxpy.entry_point() decorators into the python classpath, so they can be found by DNANexus
loaded_module = check_subjob_decorator()
if loaded_module:
    LOGGER.info(f'Loaded dxpy.entrypoint module {loaded_module}')


@dxpy.entry_point('main')
def main(mode: str, output_prefix: str, input_args: str, testing_script: dict, testing_directory: str) -> dict:
    """This is the :func:`main()` method for all apps/applets required by DNANexus.

    This method is written as an entry point for all possible analyses using the 'runassociationtesting' suite and is
    an access point to individual modules that are stored independently of this app/applet using git (any online repo
    can be used: GitHub, GitLab, etc.). This method performs three primary tasks:

    1. A conditional import (see :func:`conditional_import()` in ModuleLoader for more information) of the module
    requested by ":param:`mode`".

    2. Launching the requested module via the :func:`start_module()` function that MUST be implemented by all modules
    that use the RunAssociationTesting framework.

    3. Processing the output returned by the module (must be file-based) into a single tarball that can be passed back
    to the DNANexus platform.

    :param mode: The module that the user wants to load for this job.
    :param output_prefix: A prefix to name the output tarball returned by this method.
    :param input_args: A string containing options that conform to optparse specification. These are processed
        downstream in the ModuleLoader class.
    :param testing_script: Script compatible with pytest. If not null, invoke the runassociationtesting testing suite
        via :func:`test`.
    :param testing_directory: Directory containing test files if in testing mode.
    :return: A dictionary with keys equal to all outputs in 'output_spec' from dxapp.json and values equal to files
        uploaded to the DNANexus platform by the generate_linked_dx_file() method in association_resources. For this
        app, will only be a single tarball of all expected outputs with name 'output_prefix.tar.gz'.
    """

    if testing_script is not None:
        LOGGER.info('Testing mode activated...')
        if testing_directory is None:
            raise ValueError(f'Testing mode invoked but -itesting_directory not provided!')
        output_tarball = test(output_prefix, testing_script, testing_directory)
        output = {'output_tarball': dxpy.dxlink(generate_linked_dx_file(output_tarball))}

    else:
        # Define the package to search for based on the 'mode' requested
        LOGGER.info(f'Attempting to load module {mode}...')
        module_loader = conditional_import(f'{mode}.loader')
        LOGGER.info(f'Loaded {mode} {pkg_resources.get_distribution(mode).version}')

        # All packages MUST have a 'start_module' class by definition
        LOGGER.info(f'Launching mode {mode}...')
        loaded_module = module_loader(output_prefix, input_args)
        loaded_module.start_module()

        # Create tar of all possible output files
        output_tarball = Path(f'{output_prefix}.assoc_results.tar.gz')
        LOGGER.info(f'Processing and writing outputs to {output_tarball.name}...')
        tar = tarfile.open(output_tarball, 'w:gz')

        for file in loaded_module.get_outputs():
            tar.add(file)
        tar.add(MRC_LOGGER.get_log_file_path())
        tar.close()

        # Have to do 'upload_local_file' to make sure the new file is registered with dna nexus
        output = {'output_tarball': dxpy.dxlink(generate_linked_dx_file(output_tarball))}

    return output


def test(output_prefix: str, testing_script: dict, testing_directory: str) -> Path:
    """Run the runassociationtesting testing suite.

    This method is invisible to the applet and can only be accessed by using API calls via dxpy.DXApplet() on
    a local machine. See the resources in the `./test/` folder for more information on running tests.

    :param output_prefix: A prefix to name the output tarball returned by this method.
    :param testing_script: The dxfile ID of the pytest-compatible script
    :param testing_directory: The name of the folder containing test resources on the DNANexus platform
    :return: Dict of containing the pytest log in a tar.gz to ensure compatibility with the main() method returns
    """

    LOGGER.info('Launching mrcepid-runassociationtesting with the testing suite')
    dxpy.download_dxfile(dxid=testing_script['$dnanexus_link'], filename='test.py')

    # I then set an environment variable that tells pytest where the testing directory is
    os.environ['TEST_DIR'] = testing_directory
    os.environ['CI'] = '500'  # Make sure logs aren't truncated
    LOGGER.info(f'TEST_DIR environment variable set: {os.getenv("TEST_DIR")}')

    # pytest always throws an error when a test fails, which causes the entire suite to fall apart (and,
    # problematically, not return the logfile...). So we catch a runtime error if thrown by run_cmd() and then return
    # the log that (hopefully) should already exist. This will fall apart if there is an issue with run_cmd that is
    # outside of running pytest.
    out_log = Path(f'pytest.{output_prefix}.log')
    try:
        CommandExecutor().run_cmd('pytest test.py', stdout_file=out_log)
    except RuntimeError:
        pass

    output_tarball = Path(f'{output_prefix}.assoc_results.tar.gz')
    tar = tarfile.open(output_tarball, 'w:gz')
    tar.add(out_log)
    tar.add(MRC_LOGGER.get_log_file_path())
    tar.close()

    return output_tarball


dxpy.run()
