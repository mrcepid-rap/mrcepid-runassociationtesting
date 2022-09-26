#!/usr/bin/Rscript --vanilla

# Load required libraries
library(Matrix)
library(STAAR)
library(data.table)

# Read in inputs. They are:
# 1. [1] Filtered covariates + phenotypes file
# 2. [2] Phenotype name. Must be identical to some column in [1]
# 3. [3] Boolean flag for if our phenotype is binary
# 4. [4] Additional quantitative covariates
# 5. [5] Additional categorical covariates

args = commandArgs(trailingOnly = T)
covariates_file = args[1]
pheno_name = args[2]
is_binary = as.logical(args[3])
quant_covars = args[4]
cat_covars = args[5]

# Load covariates:
data_for_STAAR <- fread(covariates_file)

# Set covariates
data.cols <- c("FID", "age", "age_squared", "wes_batch", "sex", paste0("PC", seq(1,10)), pheno_name)

# Need to exclude NA individuals when running PheWAS data.
data_for_STAAR <- data_for_STAAR[!is.na(get(pheno_name))]

# Determine if there are additional covariates we need to consider:
if (quant_covars != "NULL") {
  quant_covars <- strsplit(quant_covars,",")[[1]]
  data.cols <- c(data.cols, quant_covars)
} else {
  quant_covars <- c()
}

if (cat_covars != "NULL") {
  cat_covars <- strsplit(cat_covars, ",")[[1]]
  data.cols <- c(data.cols, cat_covars)
  for (covar in cat_covars) {
    data_for_STAAR[,eval(covar):=as.character(get(covar))]
  }
} else {
  cat_covars <- c()
}

# Load GRM:
sparse_kinship <- readMM("/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx")
sparse_kinship_samples <- fread("/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt")
rownames(sparse_kinship) <- as.character(sparse_kinship_samples[,V1])
colnames(sparse_kinship) <- as.character(sparse_kinship_samples[,V1])

# Build data table for STAAR
# Remember! Not necessary to exclude individuals from the covariate / phenotype file as we take care of that in the main applet code
data_for_STAAR <- data_for_STAAR[,..data.cols]
data_for_STAAR[,FID:=as.character(FID)]

# Trim the genotypes/sparse kinship mtx down to individuals included in the covariate file:
sparse_kinship <- sparse_kinship[data_for_STAAR[,FID],data_for_STAAR[,FID]] # Pares down the GRM to the same individuals in the covariate / pheno data frame

# This is a very strange bit - STAAR won't run if the minimum value is < 0.05, so I threshold individuals to fix this:
sparse_kinship@x <- sapply(sparse_kinship@x, function(x) ifelse(x < 0.05, 0.05, x))

# Fit the null model for STAAR:
# These lines just autoformat our formula for association testing
if (length(unique(data_for_STAAR[,sex])) == 1) {
  covariates <- c("age", "age_squared","wes_batch",paste0("PC",seq(1,10)),quant_covars,cat_covars)
} else {
  covariates <- c("age", "age_squared","sex","wes_batch",paste0("PC",seq(1,10)),quant_covars,cat_covars)
}
cov.string <- paste(covariates, collapse=" + ")
formated.formula <- as.formula(paste(pheno_name, cov.string,sep=" ~ "))

# And run either a linear or logistic model according to is_binary, with the NULL set based on relatedness of samples
if (length(sparse_kinship@x[sparse_kinship@x < 0.5]) == 0) {
  cat("No related samples found, using GLM to fit STAAR null")
  if (is_binary) {
    if (nrow(data_for_STAAR[get(pheno_name) == 1]) == 0) { # For PheWAS mode need to check if we have any cases. If not, build a 'fake' null model so we can signal to step 2...
      obj_nullmodel <- list('zero_cases' = TRUE)  
    } else {
      obj_nullmodel <- fit_null_glm(formated.formula, data=data_for_STAAR, family="binomial")
    }
  } else {
    obj_nullmodel <- fit_null_glm(formated.formula, data=data_for_STAAR, family="gaussian")
  }
  obj_nullmodel$id_include = sparse_kinship@Dimnames[[1]] # This adds a variable to this S3 object make it easier to keep the same samples regardless of model type
  obj_nullmodel$corrected_for_relateds = FALSE
} else {
  cat("Related samples found, using LMM to fit STAAR null")
  if (is_binary) {
    obj_nullmodel <- tryCatch(
      {
        if (nrow(data_for_STAAR[get(pheno_name) == 1]) == 0) { # For PheWAS mode need to check if we have any cases. If not, build a 'fake' null model so we can signal to step 2...
          obj_nullmodel <- list('zero_cases' = TRUE,
                                'corrected_for_relateds' = FALSE,
                                'id_include' = sparse_kinship@Dimnames[[1]])
          obj_nullmodel
        } else {
          obj_nullmodel <- fit_null_glmmkin(formated.formula, data=data_for_STAAR, id="FID", family=binomial(link="logit"), kins = sparse_kinship)
          obj_nullmodel$corrected_for_relateds = TRUE
          obj_nullmodel$zero_cases = FALSE
          obj_nullmodel
        }
        },
      error=function(cond) {
        cat(paste0("STAAR model failed for phenotype ", pheno_name, ". Trying model uncorrected for relatedness..."))
        obj_nullmodel <- fit_null_glm(formated.formula, data=data_for_STAAR, family="binomial")
        obj_nullmodel$id_include = sparse_kinship@Dimnames[[1]] # This adds a variable to this S3 object make it easier to keep the same samples regardless of model type
        obj_nullmodel$corrected_for_relateds = FALSE
        obj_nullmodel$zero_cases = FALSE
        obj_nullmodel
        }
      )
  } else {
    obj_nullmodel <- tryCatch(
      {
        obj_nullmodel <- fit_null_glmmkin(formated.formula, data=data_for_STAAR, id="FID", family=gaussian(link="identity"), kins = sparse_kinship)
        obj_nullmodel$corrected_for_relateds = TRUE
        obj_nullmodel$zero_cases = FALSE
        obj_nullmodel
      },
      error=function(cond) {
        cat(paste0("STAAR model failed for phenotype ", pheno_name, ". Trying model uncorrected for relatedness..."))
        obj_nullmodel <- fit_null_glm(formated.formula, data=data_for_STAAR, family="gaussian")
        obj_nullmodel$id_include = sparse_kinship@Dimnames[[1]] # This adds a variable to this S3 object make it easier to keep the same samples regardless of model type
        obj_nullmodel$corrected_for_relateds = FALSE
        obj_nullmodel$zero_cases = FALSE
        obj_nullmodel
        }
    )
  }
}

# Save the null model:
saveRDS(object = obj_nullmodel, file = paste0("/test/", paste(pheno_name, "STAAR_null.rds", sep = ".")))
