import csv
import numpy as np
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError

from ..association_pack import AssociationPack
from ..association_resources import *
from ..tool_runners.glm_runner import GLMRunner
from ..thread_utility import ThreadUtility


class PheWAS:

    def __init__(self, association_pack: AssociationPack):

        self._association_pack = association_pack

        # 1. Figure out genes/SNPlist to run...
        gene_ENST_to_run = []
        if self._association_pack.is_snp_tar:
            gene_ENST_to_run.append('ENST00000000000')
        else:
            genes_to_run = self._association_pack.gene_ids
            transcripts_table = build_transcript_table()
            for gene in genes_to_run:
                gene_info = get_gene_id(gene, transcripts_table)
                gene_ENST_to_run.append(gene_info.name)

        # 2. Load the tarballs INTO separate genotypes dictionaries
        print("Loading Linear Model genotypes")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A GLM thread failed',
                                       incrementor=10,
                                       thread_factor=2)

        for tarball_prefix in association_pack.tarball_prefixes:
            thread_utility.launch_job(GLMRunner.load_tarball_linear_model,
                                      tarball_prefix = tarball_prefix,
                                      is_snp_tar = self._association_pack.is_snp_tar)
        future_results = thread_utility.collect_futures()
        genotype_packs = {}
        for result in future_results:
            tarball_prefix, genotype_dict = result
            genotype_packs[tarball_prefix] = genotype_dict

        # 3. Add a variable for all gene/mask combinations to the pandas data frame.
        # I think this is the most efficient way to do this on a large memory machine where I can store everything in
        # one massive data frame. Everything else I have tried takes a very long time!

        # First load the model dictionary
        self._model_dictionary = pd.read_csv("phenotypes_covariates.formatted.txt",
                                       sep=" ",
                                       index_col="FID",
                                       dtype={'IID': str})
        self._model_dictionary.index = self._model_dictionary.index.astype(str)

        # And then iterate through all possible combinations and add a column for mask+gene pairs
        for model in genotype_packs:
            for gene in gene_ENST_to_run:
                indv_w_var = genotype_packs[model][gene]
                var_name = '%s_%s' % (model, gene)
                # Need to reserved characters for formula writing from the string name:
                var_name = var_name.translate(
                    str.maketrans({'-': '_', '+': '_', '(': '_', ')': '_', '~': '_', '*': '_'}))
                self._model_dictionary[var_name] = np.where(self._model_dictionary.index.isin(indv_w_var), 1, 0)

        # 4. Now run the actual models in parallel
        # I think this is straight-forward?
        print("Submitting linear models...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A GLM thread failed',
                                       incrementor=10,
                                       thread_factor=2)

        # First set the base covariates
        covars = ['sex', 'age', 'age_squared', 'wes_batch'] + ['PC' + str(PC) for PC in range(1, 11)]

        if len(self._association_pack.found_quantitative_covariates) > 0:
            for covar in self._association_pack.found_quantitative_covariates:
                covars.append(covar)
        if len(self._association_pack.found_categorical_covariates) > 0:
            for covar in self._association_pack.found_categorical_covariates:
                covars.append('C(' + covar + ')')

        for pheno in self._association_pack.pheno_names:
            for model in genotype_packs:
                for gene in gene_ENST_to_run:
                    var_name = '%s_%s' % (model, gene)
                    # Need to reserved characters for formula writing from the string name:
                    var_name = var_name.translate(
                        str.maketrans({'-': '_', '+': '_', '(': '_', ')': '_', '~': '_', '*': '_'}))
                    form_formated = pheno + ' ~ ' + ' + '.join(covars) + ' + ' + var_name
                    thread_utility.launch_job(self._linear_model_phewas,
                                              id_column=var_name,
                                              formula=form_formated,
                                              phenoname=pheno,
                                              gene=gene,
                                              mask_name=model,
                                              is_binary=self._association_pack.is_binary)

        print("All GLM threads submitted")
        # Here we're collecting futures and writing the unformatted results at the same time
        fieldnames = ['ENST','maskname','pheno_name','p_val_init','n_car','n_model',
                      'p_val_full','effect','std_err','model_run']
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

        # 4. Annotate unformatted results and print final outputs
        print("Annotating Linear Model results")
        if self._association_pack.is_snp_tar:
            os.rename(association_pack.output_prefix + '.lm_stats.tmp',
                      association_pack.output_prefix + '.SNP.glm.stats.tsv')
            self.outputs = [association_pack.output_prefix + '.SNP.glm.stats.tsv']
        else:
            GLMRunner.process_linear_model_outputs(self._association_pack.output_prefix, gene_ENST_to_run)
            self.outputs = [association_pack.output_prefix + '.genes.glm.stats.tsv.gz',
                            association_pack.output_prefix + '.genes.glm.stats.tsv.gz.tbi']

        # This next line doesn't matter for standard burden tests, but is required for phewas when we don't pre-remove
        # individuals with NA phenotypes:
        # pheno_covars = pheno_covars[pheno_covars[phenotype].isna() == False]

    def _linear_model_phewas(self, id_column: str, formula: str, phenoname: str, gene: str, mask_name: str, is_binary: bool) -> dict:

        # Calculate Number of Carriers
        n_car = len(self._model_dictionary.loc[self._model_dictionary[id_column] == 1])

        # Default return:
        gene_dict = {'p_val_init': 'NA',
                     'n_car': n_car,
                     'n_model': len(self._model_dictionary['IID']),
                     'ENST': gene,
                     'maskname': mask_name,
                     'pheno_name': phenoname,
                     'p_val_full': 'NA',
                     'effect': 'NA',
                     'std_err': 'NA',
                     'model_run': 'false'}
        if is_binary:
            gene_dict['n_noncar_affected'] = len(self._model_dictionary.query(
                str.format('{has_var} == 0 & {phenoname} == 1', has_var=id_column, phenoname=phenoname)))
            gene_dict['n_noncar_unaffected'] = len(self._model_dictionary.query(
                str.format('{has_var} == 0 & {phenoname} == 0', has_var=id_column, phenoname=phenoname)))
            gene_dict['n_car_affected'] = len(self._model_dictionary.query(
                str.format('{has_var} == 1 & {phenoname} == 1', has_var=id_column, phenoname=phenoname)))
            gene_dict['n_car_unaffected'] = len(self._model_dictionary.query(
                str.format('{has_var} == 1 & {phenoname} == 0', has_var=id_column, phenoname=phenoname)))

        if n_car <= 2 or (is_binary and (gene_dict['n_car_affected'] <= 2)):
            return gene_dict
        else:
            try:
                sm_results = sm.GLM.from_formula(formula,
                                                 data=self._model_dictionary,
                                                 family=sm.families.Binomial() if is_binary else sm.families.Gaussian(),
                                                 missing='drop').fit()

                gene_dict['p_val_init'] = sm_results.pvalues[id_column]
                gene_dict['n_car'] = n_car
                gene_dict['n_model'] = sm_results.nobs
                gene_dict['ENST'] = gene
                gene_dict['maskname'] = mask_name
                gene_dict['pheno_name'] = phenoname
                gene_dict['p_val_full'] = sm_results.pvalues[id_column]
                gene_dict['effect'] = sm_results.params[id_column]
                gene_dict['std_err'] = sm_results.bse[id_column]
                gene_dict['model_run'] = 'true'

                # If we are dealing with a binary phenotype we also want to provide the "Fisher's" table
                return gene_dict

            except PerfectSeparationError:

                return gene_dict