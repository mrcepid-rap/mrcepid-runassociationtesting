#!/usr/bin/env python
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import sys
import dxpy
import tarfile

from general_utilities.association_resources import generate_linked_dx_file
from general_utilities.mrc_logger import MRCLogger

# Update system search paths to look for modules appropriately PRIOR to attempting actual import.
# We have to do this to get modules to run properly on DNANexus while still enabling easy programming
sys.path.append('/')
sys.path.append('/runassociationtesting/')

from runassociationtesting.module_loader import conditional_import

LOGGER = MRCLogger().get_logger()


@dxpy.entry_point('main')
def main(mode: str, output_prefix: str, input_args: str) -> dict:
    """This is the :func:`main()` method for all apps/applets required by DNANexus.

    This method is written as an entry point for all possible analyses using the 'runassociationtesting' suite and is
    simply an access point to individual modules that are stored independently of this app/applet on GitHub. This
    method performs three primary tasks:

    1. A conditional import (see :func:`conditional_import()` in ModuleLoader for more information) of the module
        requested by ":param:`mode`".
    2. Launching the requested module via the :func:`start_module()` function that MUST be implemented by all modules
        that use the RunAssociationTesting framework.
    3. Processing the output returned by the module into a single tarball that can be passed back to the DNANexus
        platform.

    :param mode: The module that the user wants to load for this job.
    :param output_prefix: A prefix to name the output tarball returned by this method.
    :param input_args: A string containing options that conform to optparse specification. These are processed
        downstream in the ModuleLoader class.
    :return: A dictionary with keys equal to all outputs in 'output_spec' from dxapp.json and values equal to files uploaded
        to the DNANexus platform by the generate_linked_dx_file() method in association_resources. For this app, will
        only be a single tarball of all expected outputs with name 'output_prefix.tar.gz'.
    """

    # Define the package to search for based on the 'mode' requested
    LOGGER.info(f'Attempting to load module {mode}...')
    module_loader = conditional_import(f'{mode}.loader')

    # All packages MUST have a 'start_module' class by definition
    LOGGER.info(f'Launching mode {mode}...')
    loaded_module = module_loader(output_prefix, input_args)
    loaded_module.start_module()

    # Create tar of all possible output files
    output_tarball = f'{output_prefix}.assoc_results.tar.gz'

    LOGGER.info(f'Processing and writing outputs to {output_tarball}...')
    tar = tarfile.open(output_tarball, "w:gz")
    for file in loaded_module.get_outputs():
        tar.add(file)
    tar.close()

    # Have to do 'upload_local_file' to make sure the new file is registered with dna nexus
    output = {"output_tarball": dxpy.dxlink(generate_linked_dx_file(output_tarball))}

    return output


dxpy.run()
