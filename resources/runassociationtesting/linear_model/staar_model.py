from typing import List
from runassociationtesting.association_resources import run_cmd


# Generate the NULL model for STAAR
def staar_null(phenoname: str, is_binary: bool,
               found_quantitative_covariates: List[str], found_categorical_covariates: List[str]) -> None:

    # I have made a custom script in order to generate the STAAR Null model:
    # located in /usr/bin/runSTAAR_Null.R
    # This generates an RDS output file containing the NULL model
    # See the README.md for more information on these parameters
    cmd = "Rscript /prog/runSTAAR_Null.R " \
          "/test/phenotypes_covariates.formatted.txt " + \
          phenoname + " " + \
          str(is_binary)
    if len(found_quantitative_covariates) > 0:
        cmd += " " + ','.join(found_quantitative_covariates)
    else:
        cmd += " NULL"
    if len(found_categorical_covariates) > 0:
        cmd += " " + ','.join(found_categorical_covariates)
    else:
        cmd += " NULL"
    run_cmd(cmd, True)


# Run rare variant association testing using STAAR
# Returns the finished chromosome to aid in output file creation
def staar_genes(tarball_prefix: str, chromosome: str, phenoname: str, has_gene_info: bool) -> tuple:

    # I have made a custom script in order to run the actual per-gene model in STAAR:
    # located in /usr/bin/runSTAAR_Genes.R
    # This generates a text output file of p.values
    # See the README.md for more information on these parameters
    cmd = f'Rscript /prog/runSTAAR_Genes.R ' \
            f'/test/{tarball_prefix}.{chromosome}.STAAR.matrix.rds ' \
            f'/test/{tarball_prefix}.{chromosome}.variants_table.STAAR.tsv ' \
            f'/test/{phenoname}.STAAR_null.rds ' + \
            f'{phenoname} ' \
            f'{tarball_prefix} ' \
            f'{chromosome} '

    # If a subset of genes has been requested, do it here.
    if has_gene_info:
        cmd += f'/test/staar.gene_list'
    else:
        cmd += f'none'  # This is always none when doing a genome-wide study.

    run_cmd(cmd, True)

    return tarball_prefix, chromosome, phenoname
