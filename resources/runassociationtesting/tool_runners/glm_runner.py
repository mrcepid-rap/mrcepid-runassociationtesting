import csv

import statsmodels.api as sm
import numpy as np

from os.path import exists
from ..association_resources import *
from ..thread_utility import ThreadUtility


class LinearModelPack:

    # This is a helper class to store information for running linear models:
    # 'phenotypes': Phenotypes/covariates for every individual in a Pandas df
    # 'model_family': family from statsmodels families
    # 'model_formula': Formatted formula for all linear models
    # 'null_model': a null model containing just IID and residuals
    # 'n_model': Number of individuals with values in the pheno/covariate file (do this here to save compute)
    # 'genotypes': A dictionary storing "genotypes" for all masks provided to this tool (empty at instantiation, is filled by self._load_tarball_linear_model)
    def __init__(self, phenotypes: pd.DataFrame, phenoname: str, model_family: sm.families, model_formula: str, null_model: pd.DataFrame):

        self.phenotypes = phenotypes
        self.phenoname = phenoname
        self.model_family = model_family
        self.model_formula = model_formula
        self.null_model = null_model
        self.n_model = len(self.phenotypes)


class GLMRunner:

    def __init__(self, association_pack: AssociationPack, gene_infos: list = None):

        # Dereference constructor variables
        self._association_pack = association_pack
        self._gene_infos = gene_infos

        # 1. Do setup for the linear models.
        # This will load all variants, genes, and phenotypes into memory to allow for parallelization
        # This function returns a class of type LinearModelPack containing info for running GLMs
        print("Loading data and running null Linear Model")
        returned_pack = self._linear_model_null(self._association_pack.pheno_names[0])
        if returned_pack is None:
            raise dxpy.AppError("Phenotype " + returned_pack.phenoname + " has no individuals after filtering, exiting...")

        # 2. Load the tarballs INTO separate genotypes dictionaries
        print("Loading Linear Model genotypes")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A GLM thread failed',
                                       incrementor=10,
                                       thread_factor=2)

        for tarball_prefix in association_pack.tarball_prefixes:
            thread_utility.launch_job(self.load_tarball_linear_model,
                                      tarball_prefix = tarball_prefix,
                                      is_snp_tar = self._association_pack.is_snp_tar,
                                      is_gene_tar = self._association_pack.is_gene_tar)
        future_results = thread_utility.collect_futures()
        genotype_packs = {}
        for result in future_results:
            tarball_prefix, genotype_dict = result
            genotype_packs[tarball_prefix] = genotype_dict

        # 3. Iterate through every model / gene (in linear_model_pack['genes']) pair and run a GLM
        print("Submitting Linear Models to threads")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A GLM thread failed',
                                       incrementor=500,
                                       thread_factor=1)

        # This 'if/else' allows us to provide a subset of models that we want to run, if selected
        for model in genotype_packs:
            if self._gene_infos is None:
                for gene in genotype_packs[model].index.levels[0]: # level[0] in this DataFrame is ENST
                    thread_utility.launch_job(self._linear_model_genes,
                                              linear_model_pack = returned_pack,
                                              genotype_table = genotype_packs[model],
                                              gene = gene,
                                              mask_name = model)
            else:
                for gene_info in self._gene_infos:
                    thread_utility.launch_job(self._linear_model_genes,
                                              linear_model_pack = returned_pack,
                                              genotype_table = genotype_packs[model],
                                              gene = gene_info.name,
                                              mask_name = model)

        # As futures finish, write unformatted results:
        fieldnames = ['ENST','maskname','pheno','p_val_init','n_car','cMAC','n_model',
                      'p_val_full','effect','std_err']
        # Binary traits get an additional set of fields to describe the confusion matrix.
        if self._association_pack.is_binary:
            fieldnames.extend(['n_noncar_affected', 'n_noncar_unaffected', 'n_car_affected', 'n_car_unaffected'])

        lm_stats_file = open(association_pack.output_prefix + '.lm_stats.tmp', 'w')
        lm_stats_writer = csv.DictWriter(lm_stats_file,
                                         delimiter = "\t",
                                         fieldnames=fieldnames)

        lm_stats_writer.writeheader()
        future_results = thread_utility.collect_futures()
        for result in future_results:
            finished_gene = result
            lm_stats_writer.writerow(finished_gene)
        lm_stats_file.close()

        # 5. Annotate unformatted results and print final outputs
        print("Annotating Linear Model results")
        self.outputs = process_linear_model_outputs(self._association_pack, self._gene_infos)

    # Setup linear models:
    def _linear_model_null(self, phenotype: str) -> LinearModelPack:

        # load covariates and phenotypes
        pheno_covars = pd.read_csv("phenotypes_covariates.formatted.txt",
                                   sep=" ",
                                   index_col="FID",
                                   dtype={'IID': str})
        pheno_covars.index = pheno_covars.index.astype(str)

        # Check if a binary trait and make sure there is more than one level after filtering
        if (self._association_pack.is_binary is True and len(pheno_covars[phenotype].value_counts()) > 1) or \
                self._association_pack.is_binary is False:

            # Decide what model family to use:
            if self._association_pack.is_binary:
                family = sm.families.Binomial()
            else:
                family = sm.families.Gaussian()

            # And finally define the formula to be used by all models:
            # Make sure to define additional covariates as requested by the user...
            form_null = phenotype + ' ~ sex + age + age_squared + C(wes_batch) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10'
            if len(self._association_pack.found_quantitative_covariates) > 0:
                for covar in self._association_pack.found_quantitative_covariates:
                    form_null += ' + ' + covar
            if len(self._association_pack.found_categorical_covariates) > 0:
                for covar in self._association_pack.found_categorical_covariates:
                    form_null += ' + C(' + covar + ')'

            form_full = form_null + ' + has_var'

            print('Using the following formula for GLMs: ')
            print(form_full)

            # Build the null model and extract residuals:
            sm_results_null = sm.GLM.from_formula(form_null, data=pheno_covars, family=sm.families.Gaussian()).fit()
            null_table = sm_results_null.resid_response.to_frame()
            null_table = null_table.rename(columns={0:'resid'})

            # Start building a set of models we want to test.
            return LinearModelPack(phenotypes=pheno_covars,
                                   phenoname=phenotype,
                                   model_family=family,
                                   model_formula=form_full,
                                   null_model=null_table)
        else:
            return None

    # load genes/genetic data we want to test/use:
    # For each tarball prefix, we want to make ONE dict for efficient querying of variants
    # Note that this **MUST BE A STATIC METHOD** as it is accessed statically by the phewas class
    @staticmethod
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

    # Run rare variant association testing using GLMs
    def _linear_model_genes(self, linear_model_pack: LinearModelPack, genotype_table: pd.DataFrame, gene: str, mask_name: str) -> dict:

        # Now successively iterate through each gene and run our model:
        # I think this is straight-forward?
        # This just extracts the pandas dataframe that is relevant to this particular gene:
        # First if is to ensure that the gene is actually in the index.
        if gene in genotype_table.index.levels[0]:
            indv_w_var = genotype_table.loc[gene]

            # We have to make an internal copy as this pandas.DataFrame is NOT threadsafe...
            # (psssttt... I still am unsure if it is...)
            internal_frame = pd.DataFrame.copy(linear_model_pack.null_model)

            # And then we merge in the genotype information
            internal_frame = pd.merge(internal_frame, indv_w_var, how='left', on='FID')
            internal_frame['has_var'] = internal_frame['gt'].transform(lambda x: 0 if np.isnan(x) else x)
            internal_frame = internal_frame.drop(columns=['gt'])
            n_car = len(internal_frame.loc[internal_frame['has_var'] >= 1])
            cMAC = internal_frame['has_var'].sum()

            if n_car <= 2:
                gene_dict = {'p_val_init': 'NA',
                             'n_car': n_car,
                             'cMAC': cMAC,
                             'n_model': linear_model_pack.n_model,
                             'ENST': gene,
                             'maskname': mask_name,
                             'pheno': linear_model_pack.phenoname,
                             'p_val_full': 'NA',
                             'effect': 'NA',
                             'std_err': 'NA'}
                if self._association_pack.is_binary:
                    gene_dict['n_noncar_affected'] = 'NA'
                    gene_dict['n_noncar_unaffected'] = 'NA'
                    gene_dict['n_car_affected'] = 'NA'
                    gene_dict['n_car_unaffected'] = 'NA'
            else:
                sm_results = sm.GLM.from_formula('resid ~ has_var',
                                                 data=internal_frame,
                                                 family=sm.families.Gaussian()).fit()

                internal_frame = pd.DataFrame.copy(linear_model_pack.phenotypes)

                # And then we merge in the genotype information
                internal_frame = pd.merge(internal_frame, indv_w_var, how='left', on='FID')
                internal_frame['has_var'] = internal_frame['gt'].transform(lambda x: 0 if np.isnan(x) else x)
                internal_frame = internal_frame.drop(columns=['gt'])

                # If we get a significant result here, re-test with the full model to get accurate beta/p. value/std. err.
                # OR if we are running a phewas, always calculate the full model
                if sm_results.pvalues['has_var'] < 1e-4 or self._association_pack.mode == "extract":

                    sm_results_full = sm.GLM.from_formula(linear_model_pack.model_formula,
                                                          data=internal_frame,
                                                          family=linear_model_pack.model_family).fit()
                    gene_dict = {'p_val_init': sm_results.pvalues['has_var'],
                                 'n_car': n_car,
                                 'cMAC': cMAC,
                                 'n_model': sm_results.nobs,
                                 'ENST': gene,
                                 'maskname': mask_name,
                                 'pheno': linear_model_pack.phenoname,
                                 'p_val_full': sm_results_full.pvalues['has_var'],
                                 'effect': sm_results_full.params['has_var'],
                                 'std_err': sm_results_full.bse['has_var']}
                    # If we are dealing with a binary phenotype we also want to provide the "Fisher's" table
                    if self._association_pack.is_binary:
                        gene_dict['n_noncar_affected'] = len(internal_frame.query(str.format('has_var == 0 & {phenoname} == 1', phenoname=linear_model_pack.phenoname)))
                        gene_dict['n_noncar_unaffected'] = len(internal_frame.query(str.format('has_var == 0 & {phenoname} == 0', phenoname=linear_model_pack.phenoname)))
                        gene_dict['n_car_affected'] = len(internal_frame.query(str.format('has_var >= 1 & {phenoname} == 1', phenoname=linear_model_pack.phenoname)))
                        gene_dict['n_car_unaffected'] = len(internal_frame.query(str.format('has_var >= 1 & {phenoname} == 0', phenoname=linear_model_pack.phenoname)))

                else:
                    gene_dict = {'p_val_init': sm_results.pvalues['has_var'],
                                 'n_car': n_car,
                                 'cMAC': cMAC,
                                 'n_model': sm_results.nobs,
                                 'ENST': gene,
                                 'maskname': mask_name,
                                 'pheno': linear_model_pack.phenoname,
                                 'p_val_full': 'NA',
                                 'effect': 'NA',
                                 'std_err': 'NA'}
                    # If we are dealing with a binary phenotype we also want to provide the "Fisher's" table
                    if self._association_pack.is_binary:
                        gene_dict['n_noncar_affected'] = len(internal_frame.query(str.format('has_var == 0 & {phenoname}=="1"', phenoname=linear_model_pack.phenoname)))
                        gene_dict['n_noncar_unaffected'] = len(internal_frame.query(str.format('has_var == 0 & {phenoname}=="0"', phenoname=linear_model_pack.phenoname)))
                        gene_dict['n_car_affected'] = len(internal_frame.query(str.format('has_var >= 1 & {phenoname}=="1"', phenoname=linear_model_pack.phenoname)))
                        gene_dict['n_car_unaffected'] = len(internal_frame.query(str.format('has_var >= 1 & {phenoname}=="0"', phenoname=linear_model_pack.phenoname)))
        else:
            gene_dict = {'p_val_init': 'NA',
                         'n_car': 0,
                         'cMAC': 0,
                         'n_model': linear_model_pack.n_model,
                         'ENST': gene,
                         'maskname': mask_name,
                         'pheno': linear_model_pack.phenoname,
                         'p_val_full': 'NA',
                         'effect': 'NA',
                         'std_err': 'NA'}
            if self._association_pack.is_binary:
                gene_dict['n_noncar_affected'] = 'NA'
                gene_dict['n_noncar_unaffected'] = 'NA'
                gene_dict['n_car_affected'] = 'NA'
                gene_dict['n_car_unaffected'] = 'NA'

        return gene_dict