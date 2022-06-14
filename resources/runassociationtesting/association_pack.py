class AssociationPack:

    def __init__(self, tarball_prefixes: list, bgen_dict: dict, is_binary: bool, sex: int, threads: int,
                 run_marker_tests: bool, output_prefix: str, pheno_names: list, mode: str, is_snp_tar: bool, is_gene_tar: bool,
                 found_quantitative_covariates: list, found_categorical_covariates: list, gene_ids: list):

        self.tarball_prefixes = tarball_prefixes
        self.bgen_dict = bgen_dict
        self.is_binary = is_binary
        self.sex = sex
        self.threads = threads
        self.run_marker_tests = run_marker_tests
        self.output_prefix = output_prefix
        self.pheno_names = pheno_names
        self.mode = mode
        self.is_snp_tar = is_snp_tar
        self.is_gene_tar = is_gene_tar
        self.is_non_standard_tar = is_snp_tar or is_gene_tar
        self.found_quantitative_covariates = found_quantitative_covariates
        self.found_categorical_covariates = found_categorical_covariates
        self.gene_ids = gene_ids
