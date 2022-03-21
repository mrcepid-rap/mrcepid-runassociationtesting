import pandas as pd

from ..association_pack import AssociationPack

class REGENIERunner:

    def __init__(self, association_pack: AssociationPack):

        self._association_pack = association_pack

        # 1. Run step 1 of regenie
        self._run_regenie_step_one()

        # 2. Mask masks together to generate required REGENIE inputs for step 2

    def _run_regenie_step_one(self):

        cmd = 'regenie ' \
              '--step 1 ' \
              '--bed /test/genetics/UKBB_450K_Autosomes_QCd_WBA ' \
              '--exclude /test/genetics/UKBB_450K_Autosomes_QCd.low_MAC.snplist ' \
              '--covarFile /test/phenotypes_covariates.formatted.txt ' \
              '--phenoFile /test/phenotypes_covariates.formatted.txt ' \
              '--keep /test/SAMPLES_Include.txt ' \
              '--bsize 100 ' \
              '--lowmem ' \
              '--lowmem-prefix /test/ ' \
              '--out /test/fit_out ' \
              '--threads ' + self._association_pack.threads + ' ' \
              '--phenoCol ' + self._association_pack.pheno_name + ' '

        if len(self._association_pack.found_quantitative_covariates) > 0:
            quant_covars_join = ','.join(self._association_pack.found_quantitative_covariates)
            cmd = cmd + '--covarColList PC{1:10},age,age_squared,sex,' + quant_covars_join + ' '
        else:
            cmd = cmd + '--covarColList PC{1:10},age,age_squared,sex '

        if len(self._association_pack.found_categorical_covariates) > 0:
            cat_covars_join = ','.join(self._association_pack.found_categorical_covariates)
            cmd = cmd + '--catCovarList wes_batch,' + cat_covars_join + ' '
        else:
            cmd = cmd + '--catCovarList wes_batch '

        if self._association_pack.is_binary:
            cmd = cmd + '--bt'

    def _make_regenie_mask(self, chromosome: str):

        for tarball_prefix in self._association_pack.tarball_prefixes:
            curr_index = pd.read_csv(tarball_prefix + '.' + chromosome + '.variants_table.STAAR.tsv')

