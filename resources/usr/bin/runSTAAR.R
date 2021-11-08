#!/usr/bin/Rscript --vanilla

library(Matrix)
library(STAAR)
library(data.table)

args = commandArgs(trailingOnly = T)

matrix_file = args[1]
variants_file = args[2]
annotations_file = args[3]
covariates_file = args[4]
pheno_name = args[5]
filename_prefix = args[6]
is_binary = as.logical(args[7])

# Load RDS genotype matrix
genotypes <- readRDS(matrix_file)

# Load list of variants
variants <- fread(variants_file)
variants[,rownum:=.I]

# Load and use REGENIE annotations to get gene IDs for testing:
annotations <- fread(annotations_file, header = F)
setnames(annotations,names(annotations), c("chrID", "GENEID", "VARTYPE"))

# Merge with the variants file so we can use later:
variants[,chrID:=paste0("chr", ID)]
variants <- merge(variants, annotations, by = "chrID")

# Load covariates:
data_for_STAAR <- fread(covariates_file)

# Build data table for STAAR
# Remember! Not necessary to exclude individuals as the covariate and phenotype files take care of that
data.cols <- c("FID", "age", "batch", "sex", paste0("PC", seq(1,10)), pheno_name)
data_for_STAAR <- data_for_STAAR[,..data.cols]
data_for_STAAR[,FID:=as.character(FID)]

# Trim the genotypes down to those in the covar file:
poss <- rownames(genotypes)
poss <- poss[poss %in% data_for_STAAR[,FID]]
genotypes <- genotypes[poss, 1:ncol(genotypes)]
data_for_STAAR <- data_for_STAAR[FID %in% poss]

# Sort the matrix so it lines up with the data frame:
genotypes <- genotypes[data_for_STAAR[,FID],1:ncol(genotypes)]

# Fit the null model:
covariates <- c("age","sex","batch",paste0("PC",seq(1,10)))
cov.string <- paste(covariates, collapse=" + ")
formated.formula <- as.formula(paste(pheno_name, cov.string,sep=" ~ "))
if (is_binary) {
  obj_nullmodel <- fit_null_glm(formated.formula, data=data_for_STAAR, family="binomial")
} else {
  obj_nullmodel <- fit_null_glm(formated.formula, data=data_for_STAAR, family="gaussian")
}
# Get list of all possible genes and make a data.table to iterate through:
poss.genes <- unique(variants[,GENEID])
gene.results <- data.table(geneID=poss.genes, n.samps=nrow(genotypes))

staar.gene <- function(gene) {
    
  current_GENE <- Matrix(genotypes[,variants[GENEID == gene, rownum]])
  if (sum(current_GENE) > 1) {
    staar_result <- STAAR(genotype = current_GENE, obj_nullmodel = obj_nullmodel, rare_maf_cutoff = 1)
    return(list(staar_result$results_STAAR_O, staar_result$num_variant, staar_result$cMAC))
  } else {
    return(list(NaN, 0, 0))
  }
  
}

gene.results[,c("staar.O.p","n.var","cMAC"):=staar.gene(geneID),by=1:nrow(gene.results)]

fwrite(gene.results, paste0("/test/", filename_prefix, ".STAAR_results.tsv"), sep = "\t", quote = F, row.names = F, col.names = T, na = "NaN")





