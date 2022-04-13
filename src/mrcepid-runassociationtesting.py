#!/usr/bin/env python
# mrcepid-runnassociationtesting 0.0.1
# Generated by dx-app-wizard.
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import tarfile
import sys

# We have to do this to get modules to run properly on DNANexus while still enabling easy editing in PyCharm
import dxpy

sys.path.append('/')
sys.path.append('/runassociationtesting/')

from runassociationtesting.association_resources import *
from runassociationtesting.covariate_processor import CovariateProcessor
from runassociationtesting.ingest_data import IngestData
from runassociationtesting.tool_runners.bolt_runner import BOLTRunner
from runassociationtesting.tool_runners.glm_runner import GLMRunner
from runassociationtesting.tool_runners.saige_runner import SAIGERunner
from runassociationtesting.tool_runners.staar_runner import STAARRunner
from runassociationtesting.tool_runners.regenie_runner import REGENIERunner
from runassociationtesting.extract_variants.extract_variants import ExtractVariants
from runassociationtesting.phewas.phewas import PheWAS


@dxpy.entry_point('main')
def main(association_tarballs, tool, mode, gene_ids, is_binary, sex,
         exclusion_list, inclusion_list, phenofile, phenoname, covarfile, categorical_covariates, quantitative_covariates, output_prefix, run_marker_tests,
         bgen_index, transcript_index, base_covariates, bed_file, fam_file, bim_file, low_MAC_list, sparse_grm, sparse_grm_sample):

    # Check required options fit before running anything:
    if mode == 'burden':
        if tool is None:
            dxpy.AppError("Must provide a tool if running 'burden' mode.")

    # Grab the data / docker resources necessary to run with the 'IngestData' class:
    ingested_data = IngestData(association_tarballs, phenofile, covarfile, inclusion_list, exclusion_list,
                               bgen_index, transcript_index, base_covariates,
                               bed_file, fam_file, bim_file, low_MAC_list, sparse_grm, sparse_grm_sample)

    # Process samples and build covariate / phenotype resources with the CovariateProcessor class
    # This does sample and covariate processing for all pipelines regardless of what we need to run
    # Also creates and returns an object of class AssocationPack that contains information necessary to run all tests
    covariate_processor = CovariateProcessor(ingested_data, phenoname, categorical_covariates, quantitative_covariates,
                                             gene_ids, sex, is_binary, run_marker_tests, output_prefix, mode)
    association_pack = covariate_processor.association_pack

    # Now do specific analysis depending on selected 'mode':
    if mode == 'burden':
        if association_pack.is_snp_tar:
            dxpy.AppError("Burden tests currently do not allow for a SNP-based tar file. Please use the 'extract' or 'phewas' modes.")
        # Run the selected tool
        if tool == 'bolt':
            print("Running BOLT")
            tool_run = BOLTRunner(association_pack)
        elif tool == 'saige':
            print("Running SAIGE")
            tool_run = SAIGERunner(association_pack)
        if tool == 'staar':
            print("Running STAAR")
            tool_run = STAARRunner(association_pack)
        elif tool == 'glm':
            print("Running Linear Models")
            tool_run = GLMRunner(association_pack)
        elif tool == 'regenie':
            print("Running REGENIE")
            tool_run = REGENIERunner(association_pack)

    elif mode == 'extract':
        if association_pack.gene_ids is None and association_pack.is_snp_tar is False:
            dxpy.AppError("Must provide gene_ids if running 'extract' mode without a SNP-list tar.")
        tool_run = ExtractVariants(association_pack)

    elif mode == 'phewas':
        if association_pack.gene_ids is None and association_pack.is_snp_tar is False:
            dxpy.AppError("Must provide gene_ids if running 'extract' mode without a SNP-list tar.")
        tool_run = PheWAS(association_pack)

    # Create tar of all possible output files
    if output_prefix is None:
        output_tarball = "assoc_results.tar.gz"
    else:
        output_tarball = output_prefix + ".assoc_results.tar.gz"

    tar = tarfile.open(output_tarball, "w:gz")
    for file in tool_run.outputs:
        tar.add(file)
    tar.close()

    ## Have to do 'upload_local_file' to make sure the new file is registered with dna nexus
    output = {"output_tarball": dxpy.dxlink(dxpy.upload_local_file(output_tarball))}

    return output


dxpy.run()
