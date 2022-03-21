class AssociationPack:

    def __init__(self, tarball_prefixes, bgen_dict, is_binary, sex, threads, run_marker_tests, output_prefix,
                 pheno_name, found_quantitative_covariates, found_categorical_covariates):

        self.tarball_prefixes = tarball_prefixes
        self.bgen_dict = bgen_dict
        self.is_binary = is_binary
        self.sex = sex
        self.threads = threads
        self.run_marker_tests = run_marker_tests
        self.output_prefix = output_prefix
        self.pheno_name = pheno_name
        self.found_quantitative_covariates = found_quantitative_covariates
        self.found_categorical_covariates = found_categorical_covariates
