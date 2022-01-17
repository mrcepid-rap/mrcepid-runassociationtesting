#!/usr/bin/Rscript --vanilla

# Load required libraries
library(Matrix)
library(STAAR)
library(data.table)

# Read in inputs. They are:
# 1. [1] A sparse matrix[i][j] where i = samples and j = variants
# 2. [2] Variants file that has the same length as j with variant IDs, positions, and MAF
# 3. [3] Annotations file that has the same length as j with variant gene information
# 4. [4] Filtered covariates + phenotypes file
# 5. [5] Phenotype name. Must be identical to some column in [4]
# 6. [6] File output prefix
# 7. [7] Current chromosome being run for output purposes
# 8. [8] Boolean flag for if our phenotype is binary
args = commandArgs(trailingOnly = T)
matrix_file = args[1]
variants_file = args[2]
annotations_file = args[3]
covariates_file = args[4]
pheno_name = args[5]
filename_prefix = args[6]
chromosome = args[7]
is_binary = as.logical(args[8])

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

# Load GRM:
sparse_kinship <- readMM("/test/genetics/sparseGRM_450K_Autosomes_QCd_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx")
sparse_kinship_samples <- fread("/test/genetics/sparseGRM_450K_Autosomes_QCd_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
rownames(sparse_kinship) <- as.character(sparse_kinship_samples[,V1])
colnames(sparse_kinship) <- as.character(sparse_kinship_samples[,V1])

# Build data table for STAAR
# Remember! Not necessary to exclude individuals from the covariate / phenotype file as we take care of that in the main applet code
data.cols <- c("FID", "age", "batch", "wes_batch", "sex", paste0("PC", seq(1,10)), pheno_name)
data_for_STAAR <- data_for_STAAR[,..data.cols]
data_for_STAAR[,FID:=as.character(FID)]

# Trim the genotypes/sparse kinship mtx down to individuals included in the covariate file:
poss <- rownames(genotypes) # this gets possible individuals in the genotype matrix
poss <- poss[poss %in% data_for_STAAR[,FID]] # Subset to those samples that exist in the covariate / pheno file
genotypes <- genotypes[poss, 1:ncol(genotypes)] # And then use that list to pare down the genotype matrix
data_for_STAAR <- data_for_STAAR[FID %in% poss] # Just make sure everything is the same between the genotype matrix and covariate / pheno file
genotypes <- genotypes[data_for_STAAR[,FID],1:ncol(genotypes)] # This sorts the genotype matrix so it lines up with the covariate / pheno data frame
sparse_kinship <- sparse_kinship[data_for_STAAR[,FID],data_for_STAAR[,FID]] # Pares down the GRM to the same individuals in the covariate / pheno data frame

# Fit the null model for STAAR:
# These lines just autoformat our formula for association testing
if (length(unique(data_for_STAAR[,sex])) == 1) {
  covariates <- c("age","batch","wes_batch",paste0("PC",seq(1,10)))
} else {
  covariates <- c("age","sex","batch","wes_batch",paste0("PC",seq(1,10)))
}
cov.string <- paste(covariates, collapse=" + ")
formated.formula <- as.formula(paste(pheno_name, cov.string,sep=" ~ "))
# And run either a linear or logistic model according to is_binary
if (is_binary) {
  obj_nullmodel <- fit_null_glmmkin(formated.formula, data=data_for_STAAR, id="FID", family=binomial(link="logit"), kins = sparse_kinship)
} else {
  obj_nullmodel <- fit_null_glmmkin(formated.formula, data=data_for_STAAR, id="FID", family=gaussian(link="identity"), kins = sparse_kinship)
}
# Get list of all possible genes and make a data.table to iterate through:
poss.genes <- unique(variants[,GENEID])
gene.results <- data.table(geneID=poss.genes, n.samps=nrow(genotypes))

# This function just takes one gene and runs it through STAAR
staar.gene <- function(gene) {
  # Grabs all columns in the genotype matrix that are for the given gene
  current_GENE <- Matrix(genotypes[,variants[GENEID == gene, rownum]])

  # I don't exclude variants that don't exist in the subset of individuals with a given phenotype
  # So we have to check here how many variants we actually have for current_GENE
  tot_vars <- 0
  for (i in 1:ncol(current_GENE)) {
    if (sum(current_GENE[,i]) > 0) {
      tot_vars <- tot_vars + 1
    }
  }
  cMAC <- sum(current_GENE)
  
  # Only run STAAR if there is greater than one non-ref variant
  if (tot_vars > 1) {
    staar_result <- STAAR(genotype = current_GENE, obj_nullmodel = obj_nullmodel, rare_maf_cutoff = 1)
    return(list(staar_result$results_STAAR_O, tot_vars, cMAC))
  } else {
    # Else just return NaN to indicate the test was not run
    return(list(NaN, tot_vars, cMAC))
  }
  
}

# This just uses data.frame functionality to run the function staar.gene on each row of the table (i.e. each gene)
gene.results[,c("staar.O.p","n.var","cMAC","n.var.calc","cMAC.calc"):=staar.gene(geneID),by=1:nrow(gene.results)]

# And write the final output table
fwrite(gene.results, paste0("/test/", paste(filename_prefix, chromosome, "STAAR_results.tsv", sep = ".")), sep = "\t", quote = F, row.names = F, col.names = T, na = "NaN")

