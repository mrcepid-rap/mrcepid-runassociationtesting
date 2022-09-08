from os.path import exists

import dxpy
import numpy as np
import pandas as pd
import statsmodels.api as sm

from runassociationtesting.association_resources import get_chromosomes, run_cmd
from runassociationtesting.ingest_data import AssociationPack


class LinearModelPack:

    # This is a helper class to store information for running linear models
    #
    # phenotypes: Phenotypes/covariates for every individual in a Pandas df
    # phenoname: name of the phenotype of interest
    # model_family: family from statsmodels families
    # model_formula: Formatted formula for all linear models
    # null_model: a null model containing just IID and residuals
    # n_model: Number of individuals with values in the pheno/covariate file (do this here to save compute)

    def __init__(self, phenotypes: pd.DataFrame, phenoname: str, model_family: sm.families, model_formula: str,
                 null_model: pd.DataFrame):

        self.phenotypes = phenotypes
        self.phenoname = phenoname
        self.model_family = model_family
        self.model_formula = model_formula
        self.null_model = null_model
        self.n_model = len(self.phenotypes)


class LinearModelResult:

    def __init__(self, p_val_init: float, n_car: int, cMAC: int, n_model: int, ENST: str, maskname: str,
                 pheno_name: str, p_val_full: float, effect: float, std_err: float):

        self.p_val_init = p_val_init
        self.n_car = n_car
        self.cMAC = cMAC
        self.n_model = n_model
        self.ENST = ENST
        self.maskname = maskname
        self.pheno_name = pheno_name
        self.p_val_full = p_val_full
        self.effect = effect
        self.std_err = std_err
        self.n_noncar_affected = 0
        self.n_noncar_unaffected = 0
        self.n_car_affected = 0
        self.n_car_unaffected = 0

    def set_carrier_stats(self, n_noncar_affected: int, n_noncar_unaffected: int,
                          n_car_affected: int, n_car_unaffected: int):

        self.n_noncar_affected = n_noncar_affected
        self.n_noncar_unaffected = n_noncar_unaffected
        self.n_car_affected = n_car_affected
        self.n_car_unaffected = n_car_unaffected

    def todict(self):
        return {'p_val_init': self.p_val_init,
                'n_car': self.n_car,
                'cMAC': self.cMAC,
                'n_model': self.n_model,
                'ENST': self.ENST,
                'maskname': self.maskname,
                'pheno_name': self.pheno_name,
                'p_val_full': self.p_val_full,
                'effect': self.effect,
                'std_err': self.std_err,
                'n_noncar_affected': self.n_noncar_affected,
                'n_noncar_unaffected': self.n_noncar_unaffected,
                'n_car_affected': self.n_car_affected,
                'n_car_unaffected': self.n_car_unaffected}


# Setup linear models:
def linear_model_null(phenotype: str, association_pack: AssociationPack) -> LinearModelPack:

    # load covariates and phenotypes
    pheno_covars = pd.read_csv("phenotypes_covariates.formatted.txt",
                               sep=" ",
                               index_col="FID",
                               dtype={'IID': str})
    pheno_covars.index = pheno_covars.index.astype(str)

    # Check if a binary trait and make sure there is more than one level after filtering
    if (association_pack.is_binary is True and len(pheno_covars[phenotype].value_counts()) > 1) or \
            association_pack.is_binary is False:

        # Decide what model family to use:
        if association_pack.is_binary:
            family = sm.families.Binomial()
        else:
            family = sm.families.Gaussian()

        # And finally define the formula to be used by all models:
        # Make sure to define additional covariates as requested by the user...
        form_null = f'{phenotype} ~ sex + age + age_squared + C(wes_batch) ' \
                    f'{" + ".join(["PC%s" % x for x in range(1, 11)])}'
        if len(association_pack.found_quantitative_covariates) > 0:
            for covar in association_pack.found_quantitative_covariates:
                form_null += ' + ' + covar
        if len(association_pack.found_categorical_covariates) > 0:
            for covar in association_pack.found_categorical_covariates:
                form_null += ' + C(' + covar + ')'

        form_full = form_null + ' + has_var'

        print('Using the following formula for GLMs: ')
        print(form_full)

        # Build the null model and extract residuals:
        sm_results_null = sm.GLM.from_formula(form_null, data=pheno_covars, family=sm.families.Gaussian()).fit()
        null_table = sm_results_null.resid_response.to_frame()
        null_table = null_table.rename(columns={0: 'resid'})

        # Start building a set of models we want to test.
        return LinearModelPack(phenotypes=pheno_covars,
                               phenoname=phenotype,
                               model_family=family,
                               model_formula=form_full,
                               null_model=null_table)
    else:
        raise dxpy.AppError(f'Phenotype {phenotype} has no individuals after filtering, exiting...')


# load genes/genetic data we want to test/use:
# For each tarball prefix, we want to make ONE dict for efficient querying of variants
def load_tarball_linear_model(tarball_prefix: str, is_snp_tar: bool, is_gene_tar: bool) -> tuple:

    print("Loading tarball prefix: " + tarball_prefix)
    geno_tables = []
    for chromosome in get_chromosomes(is_snp_tar, is_gene_tar):
        # This handles the genes that we need to test:
        if exists(tarball_prefix + "." + chromosome + ".STAAR.matrix.rds"):
            # This handles the actual genetic data:
            # This just makes a sparse matrix with columns: sample_id, gene name, genotype, ENST
            # All information is derived from the sparse STAAR matrix files
            cmd = "Rscript /prog/sparseMatrixProcessor.R " + \
                  "/test/" + tarball_prefix + "." + chromosome + ".STAAR.matrix.rds " + \
                  "/test/" + tarball_prefix + "." + chromosome + ".variants_table.STAAR.tsv " + \
                  tarball_prefix + " " + \
                  chromosome
            run_cmd(cmd, True)

            # And read in the resulting table
            geno_table = pd.read_csv(tarball_prefix + "." + chromosome + ".lm_sparse_matrix.tsv",
                                     sep="\t",
                                     dtype={'FID': str})

            # What I understand here is that 'sum()' will only sum on numeric columns, so I don't have to worry
            # about the varID column being drawn in
            geno_table = geno_table.groupby(['ENST', 'FID']).sum()
            geno_tables.append(geno_table)

    # And concatenate the final data_frame together:
    # The structure is a multi-indexed pandas DataFrame, where:
    # Index 0 = ENST
    # Index 1 = FID
    genetic_data = pd.concat(geno_tables)
    print("Finished loading tarball prefix: " + tarball_prefix)

    return tarball_prefix, genetic_data


# Helper function for run_linear_model() that makes a copy of a model and adds individuals w or w/o a variant
def add_individuals_with_variant(model_frame: pd.DataFrame, indv_w_var: pd.DataFrame):

    internal_frame = pd.DataFrame.copy(model_frame)

    # And then we merge in the genotype information
    internal_frame = pd.merge(internal_frame, indv_w_var, how='left', on='FID')
    internal_frame['has_var'] = internal_frame['gt'].transform(lambda x: 0 if np.isnan(x) else x)
    internal_frame = internal_frame.drop(columns=['gt'])

    return internal_frame


def get_null_result(n_car: int, cMAC: int, n_indvs: int, gene_name: str, mask_name: str,
                    pheno_name: str) -> LinearModelResult:
    result = LinearModelResult(p_val_init=float('nan'),
                               n_car=n_car,
                               cMAC=cMAC,
                               n_model=n_indvs,
                               ENST=gene_name,
                               maskname=mask_name,
                               pheno_name=pheno_name,
                               p_val_full=float('nan'),
                               effect=float('nan'),
                               std_err=float('nan'))

    return result


# Run association testing using GLMs
def run_linear_model(linear_model_pack: LinearModelPack, genotype_table: pd.DataFrame, gene: str,
                     mask_name: str, is_binary: bool, always_run_corrected: bool = False) -> LinearModelResult:

    # Now successively iterate through each gene and run our model:
    # I think this is straight-forward?
    # This just extracts the pandas dataframe that is relevant to this particular gene:
    # First if is to ensure that the gene is actually in the index.
    # Note that the typing error is due to pandas having two instances of DataFrame, one with multiple indices, and
    # one with a single index. The former allows 'levels' calls, the latter does not.
    if gene in genotype_table.index.levels[0]:
        indv_w_var = genotype_table.loc[gene]

        # We have to make an internal copy as this pandas.DataFrame is NOT threadsafe...
        # (psssttt... I still am unsure if it is...)
        internal_frame = add_individuals_with_variant(linear_model_pack.null_model, indv_w_var)
        n_car = len(internal_frame.loc[internal_frame['has_var'] >= 1])
        cMAC = internal_frame['has_var'].sum()

        if n_car <= 2:
            gene_dict = get_null_result(n_car, cMAC, linear_model_pack.n_model, gene, mask_name,
                                        linear_model_pack.phenoname)
        else:
            sm_results = sm.GLM.from_formula('resid ~ has_var',
                                             data=internal_frame,
                                             family=sm.families.Gaussian()).fit()

            internal_frame = add_individuals_with_variant(linear_model_pack.phenotypes, indv_w_var)

            # If we get a significant result here, re-test with the full model to get accurate
            # beta/p. value/std. err. OR if we are running a phewas, always calculate the full model
            if sm_results.pvalues['has_var'] < 1e-4 or always_run_corrected:

                sm_results_full = sm.GLM.from_formula(linear_model_pack.model_formula,
                                                      data=internal_frame,
                                                      family=linear_model_pack.model_family).fit()

                gene_dict = LinearModelResult(p_val_init=sm_results.pvalues['has_var'],
                                              n_car=n_car,
                                              cMAC=cMAC,
                                              n_model=sm_results.nobs,
                                              ENST=gene,
                                              maskname=mask_name,
                                              pheno_name=linear_model_pack.phenoname,
                                              p_val_full=sm_results_full.pvalues['has_var'],
                                              effect=sm_results_full.params['has_var'],
                                              std_err=sm_results_full.bse['has_var'])
                # If we are dealing with a binary phenotype we also want to provide the "Fisher's" table
                if is_binary:
                    phenoname = linear_model_pack.phenoname
                    gene_dict.set_carrier_stats(
                        n_noncar_affected=len(internal_frame.query(f'has_var == 0 & {phenoname} == 1')),
                        n_noncar_unaffected=len(internal_frame.query(f'has_var == 0 & {phenoname} == 0')),
                        n_car_affected=len(internal_frame.query(f'has_var >= 1 & {phenoname} == 1')),
                        n_car_unaffected=len(internal_frame.query(f'has_var >= 1 & {phenoname} == 0')))

            else:
                gene_dict = LinearModelResult(p_val_init=sm_results.pvalues['has_var'],
                                              n_car=n_car,
                                              cMAC=cMAC,
                                              n_model=sm_results.nobs,
                                              ENST=gene,
                                              maskname=mask_name,
                                              pheno_name=linear_model_pack.phenoname,
                                              p_val_full=float('nan'),
                                              effect=float('nan'),
                                              std_err=float('nan'))
                # If we are dealing with a binary phenotype we also want to provide the "Fisher's" table
                if is_binary:
                    phenoname = linear_model_pack.phenoname
                    gene_dict.set_carrier_stats(
                        n_noncar_affected=len(internal_frame.query(f'has_var == 0 & {phenoname} == 1')),
                        n_noncar_unaffected=len(internal_frame.query(f'has_var == 0 & {phenoname} == 0')),
                        n_car_affected=len(internal_frame.query(f'has_var >= 1 & {phenoname} == 1')),
                        n_car_unaffected=len(internal_frame.query(f'has_var >= 1 & {phenoname} == 0')))
    else:
        gene_dict = get_null_result(0, 0, linear_model_pack.n_model, gene, mask_name,
                                    linear_model_pack.phenoname)

    return gene_dict
