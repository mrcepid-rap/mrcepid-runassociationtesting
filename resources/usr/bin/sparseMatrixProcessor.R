#!/usr/bin/Rscript --vanilla

# Load required libraries
library(Matrix)
library(data.table)

# Read in inputs. They are:
# 1. [1] A sparse matrix[i][j] where i = samples and j = variants
# 2. [2] Variants file that has the same length as j with variant IDs, positions, and MAF
# 3. [3] File output prefix
# 4. [4] Current chromosome being run for output purposes

args = commandArgs(trailingOnly = T)
matrix_file = args[1]
variants_file = args[2]
filename_prefix = args[3]
chromosome = args[4]

sparse_matrix <- readRDS(matrix_file)
variant_index = fread(variants_file)

sparse_matrix <- data.table(FID=rownames(sparse_matrix)[sparse_matrix@i+1],
                   varID=colnames(sparse_matrix)[sparse_matrix@j+1],
                   gt=sparse_matrix@x)

sparse_matrix <- merge(sparse_matrix, variant_index[,c("varID","ENST")], by = "varID")

fwrite(sparse_matrix, paste0("/test/", paste(filename_prefix, chromosome, "lm_sparse_matrix.tsv", sep = ".")), sep = "\t", quote = F, row.names = F, col.names = T, na = "NaN")
