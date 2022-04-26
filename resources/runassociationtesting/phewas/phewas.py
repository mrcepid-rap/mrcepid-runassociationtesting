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

        # 3. Now build individual phenotype dictionaries:
        # Successively iterate through each gene and run our model:
        # I think this is straight-forward?
        print("Submitting linear models...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A GLM thread failed',
                                       incrementor=10,
                                       thread_factor=2)

        for pheno in self._association_pack.pheno_names:
            covars = ['IID', pheno, 'sex', 'age', 'age_squared', 'wes_batch'] + ['PC' + str(PC) for PC in range(1,11)]
            form_formated = pheno + ' ~ has_var + sex + age + age_squared + C(wes_batch) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10'
            if len(self._association_pack.found_quantitative_covariates) > 0:
                for covar in self._association_pack.found_quantitative_covariates:
                    form_formated += ' + ' + covar
                    covars.append(covar)
            if len(self._association_pack.found_categorical_covariates) > 0:
                for covar in self._association_pack.found_categorical_covariates:
                    form_formated += ' + C(' + covar + ')'
                    covars.append(covar)

            for model in genotype_packs:
                for gene in gene_ENST_to_run:
                    thread_utility.launch_job(self._linear_model_phewas,
                                              indv_w_var=genotype_packs[model][gene],
                                              formula=form_formated,
                                              phenoname=pheno,
                                              gene=gene,
                                              mask_name=model,
                                              is_binary=self._association_pack.is_binary)

        print("All GLM threads submitted")
        # Here we're collecting futures and writing the unformatted results at the same time
        fieldnames = ['ENST','maskname','pheno_name','p_val_init','n_car','n_model',
                      'p_val_full','effect','std_err']
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
            self._process_linear_model_outputs(genes_to_run)
            self.outputs = [association_pack.output_prefix + '.genes.glm.stats.tsv.gz',
                            association_pack.output_prefix + '.genes.glm.stats.tsv.gz.tbi']

        # This next line doesn't matter for standard burden tests, but is required for phewas when we don't pre-remove
        # individuals with NA phenotypes:
        # pheno_covars = pheno_covars[pheno_covars[phenotype].isna() == False]

    @staticmethod
    def _linear_model_phewas(indv_w_var: list, formula: str, phenoname: str, gene: str, mask_name: str, is_binary: bool) -> dict:

        model_dictionary = pd.read_csv("phenotypes_covariates.formatted.txt",
                                   sep=" ",
                                   index_col="FID",
                                   dtype={'IID': str})
        model_dictionary.index = model_dictionary.index.astype(str)

        # Remove phenotype NAs
        model_dictionary = model_dictionary[model_dictionary[phenoname].isna() == False]

        # Add a column for having a variable
        model_dictionary['has_var'] = np.where(model_dictionary.index.isin(indv_w_var), 1, 0)

        # Calculate Number of Carriers
        n_car = len(model_dictionary.loc[model_dictionary['has_var'] == 1])

        # Default return:
        gene_dict = {'p_val_init': 'NA',
                     'n_car': n_car,
                     'n_model': len(model_dictionary['IID']),
                     'ENST': gene,
                     'maskname': mask_name,
                     'pheno_name': phenoname,
                     'p_val_full': 'NA',
                     'effect': 'NA',
                     'std_err': 'NA'}
        if is_binary:
            gene_dict['n_noncar_affected'] = 'NA'
            gene_dict['n_noncar_unaffected'] = 'NA'
            gene_dict['n_car_affected'] = 'NA'
            gene_dict['n_car_unaffected'] = 'NA'

        if n_car <= 2:
            return gene_dict
        else:
            try:
                sm_results = sm.GLM.from_formula(formula,
                                                 data=model_dictionary,
                                                 family=sm.families.Binomial() if is_binary else sm.families.Gaussian()).fit()

                gene_dict = {'p_val_init': sm_results.pvalues['has_var'],
                             'n_car': n_car,
                             'n_model': sm_results.nobs,
                             'ENST': gene,
                             'maskname': mask_name,
                             'pheno_name': phenoname,
                             'p_val_full': sm_results.pvalues['has_var'],
                             'effect': sm_results.params['has_var'],
                             'std_err': sm_results.bse['has_var']}

                # If we are dealing with a binary phenotype we also want to provide the "Fisher's" table
                if is_binary:
                    gene_dict['n_noncar_affected'] = len(model_dictionary.query(str.format('has_var == 0 & {phenoname} == 1', phenoname=phenoname)))
                    gene_dict['n_noncar_unaffected'] = len(model_dictionary.query(str.format('has_var == 0 & {phenoname} == 0', phenoname=phenoname)))
                    gene_dict['n_car_affected'] = len(model_dictionary.query(str.format('has_var == 1 & {phenoname} == 1', phenoname=phenoname)))
                    gene_dict['n_car_unaffected'] = len(model_dictionary.query(str.format('has_var == 1 & {phenoname} == 0', phenoname=phenoname)))

                return gene_dict

            except PerfectSeparationError:

                return gene_dict