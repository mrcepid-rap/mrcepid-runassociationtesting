#!/usr/bin/env Rscript

# This script is intended for the generation of test data required for the mrcepid-runassociationtesting test suite.

# Load necessary modules
library(data.table)
library(tidyverse)


# Make sure we get reproducible results
set.seed(1234)

# Now we generate a set of files:
# 1. A base coviarate file
# 2. A phenotype and covariate file
# 3. Permutations of these files to do testing (see below)

# Set fake IDs
# We are first going to create 1000 'participants'. The first 995 will be full participants,
eid <- sample(x=c(1000000:6999999), size = 1000, replace = FALSE)
base_covar_test <- data.table(eid=eid)

# Set PCs (Value of the PC does not matter for this test, just need to ensure it's a float
for (PC in c(1:40)) {
  PC_data <- c(rep(NA, 5), rnorm(995, mean=sample(c(-5:5), 1), sd=sample(seq(0.1,2, by=0.1), 1)))
  col_name <- paste0("22009-0.",PC)
  base_covar_test[,eval(col_name):=PC_data]
}

# Set sex, genetics QC pass, and wes batch (need to strictly define individuals for testing)
# The 'pyramids' below indicate the number of samples after we simulate sex, genetics pass, and WES batch
# 5 NA;           547 females;                               448 males
# 5 NA;       500 pass, 47 fail;                         400 pass, 48 fail
# 5 NA; 200 200k, 230 450k, 50 470k, 20 NA;       180 200k, 200 450k, 30 470k, 10 NA;
#       20 200k, 10 450k, 5 470k, 2 NA            20 200k, 10 450k, 5 470k, 3 NA
base_covar_test[,`22001-0.0`:=c(rep(NA, 5), rep(0,547),rep(1,448))]
base_covar_test[,genetics_qc_pass:=c(rep(0, 5), rep(1, 500), rep(0, 47), rep(1, 400), rep(0, 48))]
base_covar_test[,wes.batch:=c(rep(NA, 5),
                              rep("200k", 200),rep("450k", 230),rep("470k", 50), rep(NA, 20), # Female pass
                              rep("200k", 20),rep("450k", 10),rep("470k", 5), rep(NA, 2), # Female fail
                              rep("200k", 180),rep("450k", 200),rep("470k", 30), rep(NA, 10), # Male pass
                              rep("200k", 20),rep("450k", 10),rep("470k", 5), rep(NA, 3))] # Male fail

# Set age and array (Can be random as no filtering is done at this stage)
base_covar_test[,`21003-0.0`:=sample(c(40:70), size=1000, replace = T)]
base_covar_test[,`22000-0.0`:=c(rep('believenan', 5), sample(c(paste0('axiom',c(1:95)),paste0('bileve', c(1:11))), size=995, replace = T))]

# Write base covariates (in tab-delim and space-delim formats, and with scrambled columns)
base_covar_test <- sample_n(base_covar_test, 1000)
cols <- c("eid",paste0("22009-0.",1:40),"22000-0.0","22001-0.0","21003-0.0","22001-0.0","21003-0.0","22001-0.0","wes.batch","genetics_qc_pass")
cols_scrambled <- sample(cols)
fwrite(base_covar_test[,..cols], "test_data/base_covar.tsv", sep="\t", row.names=F, col.names=T, quote=F)
fwrite(base_covar_test[,..cols_scrambled], "test_data/base_covar.scrambled.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Simulate phenotypes for the whole population
base_covar_test[,quantPheno1:=rnorm(1000, mean=0, sd=2)]
base_covar_test[,quantPheno2:=rpois(1000, lambda = 5)]
base_covar_test[,catPheno3:=sample(c(0,1), size = 1000, replace = T, prob = c(0.9,0.1))]

# Simulate covariates for the whole population
base_covar_test[,quantCovar1:=rnorm(1000)]
base_covar_test[,quantCovar2:=rpois(1000, 10)]
base_covar_test[,catCovar3:=sample(x=c('red','blue','green'), size=1000, replace=TRUE)]

# Build 3 pheno/covar files and generate missing test_data
# 10 (5 female/5 male) individuals INDEPENDENTLY set to NA for each phenotype/covariate to simulate missing test_data\
# But only for individuals not already NA to ensure easy testing of numbers
pheno_covar <- base_covar_test[genetics_qc_pass == 1 & !is.na(`22001-0.0`),c("eid","22001-0.0","quantPheno1","quantPheno2","catPheno3","quantCovar1","quantCovar2","catCovar3")]

exclude_pheno_covar <- function(col) {
  na.samples.female <- pheno_covar[!is.na(quantPheno1) & !is.na(quantPheno2) & !is.na(catPheno3) & !is.na(quantCovar1) & !is.na(quantCovar2) & !is.na(catCovar3) & `22001-0.0` == 0, sample(eid, 5, replace = F)]
  na.samples.male <- pheno_covar[!is.na(quantPheno1) & !is.na(quantPheno2) & !is.na(catPheno3) & !is.na(quantCovar1) & !is.na(quantCovar2) & !is.na(catCovar3) & `22001-0.0` == 1, sample(eid, 5, replace = F)]
  na.samples <- c(na.samples.female, na.samples.male)
  pheno_covar[,eval(col):=ifelse(eid %in% na.samples, NA, get(col))]
}

for (col in c("quantPheno1","quantPheno2","catPheno3","quantCovar1","quantCovar2","catCovar3")) {
  exclude_pheno_covar(col)
}

# Set column names correctly for this file format
pheno_covar[,FID:=eid]
pheno_covar[,IID:=eid]
pheno_covar[,eid:=NULL]

# 3 files:
# 1. Combined pheno/covar
# 2. Combined pheno/covar space-delimited
# 3. Combined pheno/covar with scrambled columns to ensure that the header is read properly
cols <- c("FID","IID","quantPheno1","quantPheno2","catPheno3","quantCovar1","quantCovar2","catCovar3")
fwrite(pheno_covar[,..cols], "test_data/pheno_covar.tsv", sep="\t", row.names=F, col.names=T, quote=F)
cols <- c("FID","IID","quantPheno1","quantPheno2","catPheno3","quantCovar1","quantCovar2","catCovar3")
fwrite(pheno_covar[,..cols], "test_data/pheno_covar.txt", sep=" ", row.names=F, col.names=T, quote=F)
cols_wrong <- c("IID","quantPheno1","quantPheno2","catPheno3","quantCovar1","quantCovar2","catCovar3")
fwrite(pheno_covar[,..cols_wrong], "test_data/pheno_covar.wrong_header.tsv", sep="\t", row.names=F, col.names=T, quote=F)
cols <- sample(cols)
fwrite(pheno_covar[,..cols], "test_data/pheno_covar.scrambled.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# 3 files:
# 1. Seperate phenotype file
# 2/3. Two phenofiles
# 4. Seperate covariate file
# 5. A covariate file missing the "FID" column to ensure the correct error is thrown
cols <- c("FID","IID","quantPheno1","quantPheno2","catPheno3")
fwrite(pheno_covar[,..cols], "test_data/pheno.tsv", sep="\t", row.names=F, col.names=T, quote=F)
cols <- c("FID","IID","quantPheno1","quantPheno2")
fwrite(pheno_covar[,..cols], "test_data/pheno.first.tsv", sep="\t", row.names=F, col.names=T, quote=F)
cols <- c("FID","IID","catPheno3")
fwrite(pheno_covar[,..cols], "test_data/pheno.second.tsv", sep="\t", row.names=F, col.names=T, quote=F)
cols <- c("FID","IID","quantCovar1","quantCovar2","catCovar3")
fwrite(pheno_covar[,..cols], "test_data/covar.tsv", sep="\t", row.names=F, col.names=T, quote=F)
cols_wrong <- c("IID","quantCovar1","quantCovar2","catCovar3")
fwrite(pheno_covar[,..cols_wrong], "test_data/covar.wrong_header.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Sample Inclusion / Exclusion
# Generate a list of individuals that have no missing data (to make for easy counting)
# And get a list of 700 individuals we want to include
samples <- pheno_covar[!is.na(quantPheno1) & !is.na(quantPheno2) & !is.na(catPheno3) & !is.na(quantCovar1) & !is.na(quantCovar2) & !is.na(catCovar3), IID]
include <- sample(samples, 700)
fwrite(data.table(include), "test_data/include.txt", row.names=F, col.names=F, quote=F)

# Sample Exclusion (we want some samples in Inclusion to make sure this works correctly when we provide both an
# inclusion and exclusion file)
exclude <- c(sample(samples[!samples %in% include], 50), sample(samples[samples %in% include], 10))
fwrite(data.table(exclude), "test_data/exclude.txt", row.names=F, col.names=F, quote=F)
