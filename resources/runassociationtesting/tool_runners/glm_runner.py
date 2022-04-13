import csv
import statsmodels.api as sm
import numpy as np

from os.path import exists
from pandas.core.series import Series
from ..association_pack import AssociationPack
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

    def __init__(self, association_pack: AssociationPack, genes_to_run: list = None):

        # 1. Do setup for the linear models.
        # This will load all variants, genes, and phenotypes into memory to allow for parallelization
        # This function returns a class of type LinearModelPack containing info for running GLMs
        # If we are doing a PheWAS mode, we set a flag to not run a null model
        print("Loading data and running null Linear Model")
        self._association_pack = association_pack
        model_packs = {}
        for pheno in self._association_pack.pheno_names:
            returned_pack = self._linear_model_null(pheno)
            if returned_pack is None:
                print("Phenotype " + pheno + " has no individuals after filtering, skipping...")
            else:
                model_packs[pheno] = returned_pack

        # 2. Load the tarballs INTO separate genotypes dictionaries
        print("Loading Linear Model genotypes")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A GLM thread failed',
                                       incrementor=10,
                                       thread_factor=2)

        for tarball_prefix in association_pack.tarball_prefixes:
            thread_utility.launch_job(self._load_tarball_linear_model,
                                      tarball_prefix = tarball_prefix,
                                      is_snp_tar = self._association_pack.is_snp_tar)
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
        for pheno in model_packs:
            for model in genotype_packs:
                if genes_to_run is None:
                    for gene in genotype_packs[model]:
                        thread_utility.launch_job(self._linear_model_genes,
                                                  linear_model_pack = model_packs[pheno],
                                                  genotype_table = genotype_packs[model],
                                                  gene = gene,
                                                  mask_name = model,
                                                  mode = self._association_pack.mode)
                else:
                    for gene in genes_to_run:
                        thread_utility.launch_job(self._linear_model_genes,
                                                  linear_model_pack = model_packs[pheno],
                                                  genotype_table = genotype_packs[model],
                                                  gene = gene,
                                                  mask_name = model,
                                                  mode = self._association_pack.mode,
                                                  is_binary = self._association_pack.is_binary)

        # 4. Write unformatted results:
        print("Writing initial Linear Model results")

        # Binary traits get an additional set of fields to describe the confusion matrix.
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

        # 5. Annotate unformatted results and print final outputs
        print("Annotating Linear Model results")
        if self._association_pack.is_snp_tar:
            os.rename(association_pack.output_prefix + '.lm_stats.tmp',
                      association_pack.output_prefix + '.SNP.glm.stats.tsv')
            self.outputs = [association_pack.output_prefix + '.SNP.glm.stats.tsv']
        else:
            self._process_linear_model_outputs(genes_to_run)
            self.outputs = [association_pack.output_prefix + '.genes.glm.stats.tsv.gz',
                            association_pack.output_prefix + '.genes.glm.stats.tsv.gz.tbi']

    # Setup linear models:
    def _linear_model_null(self, phenotype: str) -> LinearModelPack:

        # load covariates and phenotypes
        pheno_covars = pd.read_csv("phenotypes_covariates.formatted.txt",
                                   sep=" ",
                                   index_col="FID",
                                   dtype={'IID': str})
        pheno_covars.index = pheno_covars.index.astype(str)
        # This next line doesn't matter for standard burden tests, but is required for phewas when we don't pre-remove
        # individuals with NA phenotypes:
        pheno_covars = pheno_covars[pheno_covars[phenotype].isna() == False]

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
    @staticmethod
    def _load_tarball_linear_model(tarball_prefix: str, is_snp_tar: bool) -> tuple:

        print("Loading tarball prefix: " + tarball_prefix)
        geno_tables = []
        for chromosome in get_chromosomes(is_snp_tar):
            # This handles the genes that we need to test:
            if exists(tarball_prefix + "." + chromosome + ".BOLT.bgen"):
                # This handles the actual genetic data:
                # load genetic data
                # first convert into format we can use:
                cmd = "plink2 --threads 2 --bgen /test/" + tarball_prefix + "." + chromosome + ".BOLT.bgen 'ref-last' " \
                        "--sample /test/" + tarball_prefix + "." + chromosome + ".BOLT.sample " \
                        "--export bcf --out /test/lm." + tarball_prefix + "." + chromosome + " " + \
                        "--keep-fam /test/SAMPLES_Include.txt"
                run_cmd(cmd, True)

                # This just makes a sparse matrix with columns: sample_id, gene name, genotype
                cmd = "bcftools query -i \"GT='alt'\" -f \"[%SAMPLE\\t%ID\\t%GT\\n]\" " \
                        "/test/lm." + tarball_prefix + "." + chromosome + ".bcf > " \
                        "lm." + tarball_prefix + "." + chromosome + ".tsv"
                run_cmd(cmd, True)

                geno_table = pd.read_csv("lm." + tarball_prefix + "." + chromosome + ".tsv",
                                         sep = "\t",
                                         names = ['eid', 'gene', 'gt'])
                # bgen stores samples in eid_eid format, so get just the first eid
                geno_table[['eid','eid2']] = geno_table['eid'].str.split('_', 1, expand=True)
                geno_table = geno_table.drop('eid2', axis=1)

                # Get all possible genes found and aggregate across gene IDs to get a list of individuals with a given gene
                geno_table = geno_table.groupby('gene').agg({'eid': Series.to_list})
                geno_tables.append(geno_table)

        # And concatenate the final data_frame together:
        genetic_data = pd.concat(geno_tables).to_dict()
        genetic_data = genetic_data['eid'] # no idea why there is a top level 'eid' key
        print("Finished loading tarball prefix: " + tarball_prefix)

        return tarball_prefix, genetic_data

    # Run rare variant association testing using GLMs
    @staticmethod
    def _linear_model_genes(linear_model_pack: LinearModelPack, genotype_table: dict, gene: str, mask_name: str, mode: str, is_binary: bool) -> dict:

        # Now successively iterate through each gene and run our model:
        # I think this is straight-forward?
        indv_w_var = genotype_table[gene]

        # We have to make an internal copy as this pandas.DataFrame is NOT threadsafe...
        # (psssttt... I still am unsure if it is...)
        internal_frame = pd.DataFrame.copy(linear_model_pack.null_model)
        internal_frame['has_var'] = np.where(internal_frame.index.isin(indv_w_var), 1, 0)
        n_car = len(internal_frame.loc[internal_frame['has_var'] == 1])
        if n_car <= 2:
            gene_dict = {'p_val_init': 'NA',
                         'n_car': n_car,
                         'n_model': linear_model_pack.n_model,
                         'ENST': gene,
                         'maskname': mask_name,
                         'pheno_name': linear_model_pack.phenoname,
                         'p_val_full': 'NA',
                         'effect': 'NA',
                         'std_err': 'NA'}
            if is_binary:
                gene_dict['n_noncar_affected'] = 'NA'
                gene_dict['n_noncar_unaffected'] = 'NA'
                gene_dict['n_car_affected'] = 'NA'
                gene_dict['n_car_unaffected'] = 'NA'
        else:
            sm_results = sm.GLM.from_formula('resid ~ has_var',
                                             data=internal_frame,
                                             family=sm.families.Gaussian()).fit()

            internal_frame = pd.DataFrame.copy(linear_model_pack.phenotypes)
            internal_frame['has_var'] = np.where(internal_frame.index.isin(indv_w_var), 1, 0)

            # If we get a significant result here, re-test with the full model to get accurate beta/p. value/std. err.
            # OR if we are running a phewas, always calculate the full model
            if sm_results.pvalues['has_var'] < 1e-4 or mode == "phewas":

                sm_results_full = sm.GLM.from_formula(linear_model_pack.model_formula,
                                                      data=internal_frame,
                                                      family=linear_model_pack.model_family).fit()
                gene_dict = {'p_val_init': sm_results.pvalues['has_var'],
                             'n_car': n_car,
                             'n_model': sm_results.nobs,
                             'ENST': gene,
                             'maskname': mask_name,
                             'pheno_name': linear_model_pack.phenoname,
                             'p_val_full': sm_results_full.pvalues['has_var'],
                             'effect': sm_results_full.params['has_var'],
                             'std_err': sm_results_full.bse['has_var']}
                # If we are dealing with a binary phenotype we also want to provide the "Fisher's" table
                if is_binary:
                    gene_dict['n_noncar_affected'] = len(internal_frame.query(str.format('has_var == 0 & {phenoname} == 1', phenoname=linear_model_pack.phenoname)))
                    gene_dict['n_noncar_unaffected'] = len(internal_frame.query(str.format('has_var == 0 & {phenoname} == 0', phenoname=linear_model_pack.phenoname)))
                    gene_dict['n_car_affected'] = len(internal_frame.query(str.format('has_var == 1 & {phenoname} == 1', phenoname=linear_model_pack.phenoname)))
                    gene_dict['n_car_unaffected'] = len(internal_frame.query(str.format('has_var == 1 & {phenoname} == 0', phenoname=linear_model_pack.phenoname)))

            else:
                gene_dict = {'p_val_init': sm_results.pvalues['has_var'],
                             'n_car': n_car,
                             'n_model': sm_results.nobs,
                             'ENST': gene,
                             'maskname': mask_name,
                             'pheno_name': linear_model_pack.phenoname,
                             'p_val_full': 'NA',
                             'effect': 'NA',
                             'std_err': 'NA'}
                # If we are dealing with a binary phenotype we also want to provide the "Fisher's" table
                if is_binary:
                    gene_dict['n_noncar_affected'] = len(internal_frame.query(str.format('has_var == 0 & {phenoname}=="1"', phenoname=linear_model_pack.phenoname)))
                    gene_dict['n_noncar_unaffected'] = len(internal_frame.query(str.format('has_var == 0 & {phenoname}=="0"', phenoname=linear_model_pack.phenoname)))
                    gene_dict['n_car_affected'] = len(internal_frame.query(str.format('has_var == 1 & {phenoname}=="1"', phenoname=linear_model_pack.phenoname)))
                    gene_dict['n_car_unaffected'] = len(internal_frame.query(str.format('has_var == 1 & {phenoname}=="0"', phenoname=linear_model_pack.phenoname)))

        return gene_dict

    # Process/annotate Linear Model output file
    def _process_linear_model_outputs(self, genes_to_run: list):

        # read in the GLM stats file:
        glm_table = pd.read_csv(open(self._association_pack.output_prefix + ".lm_stats.tmp", 'r'), sep = "\t")

        # Now process the gene table into a useable format:
        # First read in the transcripts file
        transcripts_table = pd.read_csv(gzip.open('transcripts.tsv.gz', 'rt'), sep = "\t")
        transcripts_table = transcripts_table.rename(columns={'#chrom':'chrom'})
        transcripts_table = transcripts_table.set_index('ENST')
        transcripts_table = transcripts_table[transcripts_table['fail'] == False]
        transcripts_table = transcripts_table.drop(columns=['syn.count','fail.cat','fail'])

        # Limit to genes we care about if running only a subset:
        if genes_to_run is not None:
            transcripts_table = transcripts_table.loc[genes_to_run]

        # Test what columns we have in the 'SNP' field so we can name them...
        field_one = glm_table.iloc[1]
        field_one = field_one['maskname'].split("-")
        field_names = []
        if len(field_one) == 2: # This could be the standard naming format... check that column [1] is MAF/AC
            if 'MAF' in field_one[1] or 'AC' in field_one[1]:
                field_names.extend(['MASK','MAF'])
            else: # This means we didn't hit on MAF/AC in column [2] and a different naming convention is used...
                field_names.extend(['var1','var2'])
        else:
            for i in range(2,len(field_one) + 1):
                field_names.append('var%i' % i)

        # Process the 'SNP' column into separate fields and remove
        glm_table[field_names] = glm_table['maskname'].str.split("-",expand = True)
        glm_table = glm_table.drop(columns=['maskname'])

        # Now merge the transcripts table into the gene table to add annotation and the write
        glm_table = pd.merge(transcripts_table, glm_table, on='ENST', how="left")
        with open(self._association_pack.output_prefix + '.genes.glm.stats.tsv', 'w') as gene_out:

            # Sort just in case
            glm_table = glm_table.sort_values(by=['chrom','start','end'])

            glm_table.to_csv(path_or_buf=gene_out, index = False, sep="\t", na_rep='NA')
            gene_out.close()

            # And bgzip and tabix...
            cmd = "bgzip /test/" + self._association_pack.output_prefix + '.genes.glm.stats.tsv'
            run_cmd(cmd, True)
            cmd = "tabix -S 1 -s 2 -b 3 -e 4 /test/" + self._association_pack.output_prefix + '.genes.glm.stats.tsv.gz'
            run_cmd(cmd, True)