#!/usr/bin/Rscript --vanilla

# Load required libraries
library(Matrix)
library(STAAR)
library(data.table)

# Read in inputs. They are:
# 1. [1] A sparse matrix[i][j] where i = samples and j = variants
# 2. [2] Variants file that has the same length as j with variant IDs, positions, and MAF
# 3. [3] The null model file from runSTAAR_Null.R
# 4. [4] Pheno name for including in output table
# 5. [5] File output prefix
# 6. [6] Current chromosome being run for output purposes
# 7. [7] If only running a subset of genes, define here. Default is 'none'
args = commandArgs(trailingOnly = T)
matrix_file = args[1]
variants_file = args[2]
null_model_file = args[3]
pheno_name = args[4]
filename_prefix = args[5]
chromosome = args[6]
genes_file = args[7]

# Load RDS genotype matrix
genotypes <- readRDS(matrix_file)

# Load list of variants
variants <- fread(variants_file)

# Load the null model file:
obj_nullmodel <- readRDS(null_model_file)

# Load the gene list if necessary:
if (genes_file != "none") {
  valid.genes <- unique(variants[,ENST]) # Need a list of valid genes to ensure that we only analyse genes on the current chromosome...
  poss.genes <- fread(genes_file, header = F)
  poss.genes <- poss.genes[V1 %in% valid.genes,V1]
} else {
  poss.genes <- unique(variants[,ENST])
}

# Trim the genotypes/sparse kinship mtx down to individuals included in the null model file:
poss <- obj_nullmodel$id_include # this gets possible individuals from the null model
genotypes <- genotypes[poss, 1:ncol(genotypes)] # And then use that list to pare down the genotype matrix

# Make a data.table to iterate through that contains all genes defined in poss.genes above:
gene.results <- data.table(geneID=poss.genes, n.samps=nrow(genotypes), pheno_name=pheno_name, mask=filename_prefix, relatedness.correction=obj_nullmodel$corrected_for_relateds)
gene.results[,SNP:=paste(geneID,mask,sep="-")] # This sets an ID similar to BOLT to allow identical processing

# This function just takes one gene and runs it through STAAR
staar.gene <- function(gene) {
  # Grabs all columns in the genotype matrix that are for the given gene
  current_GENE <- Matrix(genotypes[,variants[ENST == gene, column]])

  # I don't exclude variants that don't exist in the subset of individuals with a given phenotype
  # So we have to check here how many variants we actually have for current_GENE
  tot_vars <- 0
  for (i in 1:ncol(current_GENE)) {
    if (sum(current_GENE[,i]) > 0) {
      tot_vars <- tot_vars + 1
    }
  }
  cMAC <- sum(current_GENE)
  
  # Only run STAAR if there is greater than one non-ref variant and there are cases
  if (obj_nullmodel$zero_cases == TRUE | tot_vars <= 1) {
    # Else just return NaN to indicate the test was not run
    return(list(NaN, NaN, NaN, NaN, tot_vars, cMAC))
  } else {
    staar_result <- STAAR(genotype = current_GENE, obj_nullmodel = obj_nullmodel, rare_maf_cutoff = 1)
    return(list(staar_result$results_STAAR_O,
                staar_result$results_STAAR_S_1_25[[1]],
                staar_result$results_STAAR_B_1_1[[1]],
                staar_result$results_STAAR_A_1_25[[1]],
                tot_vars,
                cMAC))
  }
  
}

# This just uses data.frame functionality to run the function staar.gene on each row of the table (i.e. each gene)
gene.results[,c("staar.O.p","staar.SKAT.p","staar.burden.p","staar.ACAT.p","n_var","cMAC"):=staar.gene(geneID),by=1:nrow(gene.results)]

# And write the final output table
# Remove the geneID and mask values as they are redundant at this point:
gene.results[,geneID:=NULL]
gene.results[,mask:=NULL]
fwrite(gene.results, paste0("/test/", paste(filename_prefix, pheno_name, chromosome, "STAAR_results.tsv", sep = ".")), sep = "\t", quote = F, row.names = F, col.names = T, na = "NaN")