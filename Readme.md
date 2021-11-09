# RunAssociationTesting (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
   * [Background](#background)
      + [1. BOLT-LMM](#1-bolt-lmm)
      + [2. SAIGE-GENE](#2-saige-gene)
      + [3. STAAR](#3-staar)
      + [4. Generalised Linear Models (GLMs)](#4-generalised-linear-models-glms)
   * [Dependencies](#dependencies)
      + [Docker](#docker)
      + [Resource Files](#resource-files)
- [Methodology](#methodology)
   * [Covariate processing](#covariate-processing)
   * [BOLT](#bolt)
      + [Inputs](#inputs)
      + [Command Line Example](#command-line-example)
      + [Outputs](#outputs)
   * [SAIGE-GENE](#saige-gene)
      + [Inputs](#inputs-1)
      + [Command Line Example](#command-line-example-1)
      + [Outputs](#outputs-1)
   * [STAAR](#staar)
      + [Inputs](#inputs-2)
      + [Command line example](#command-line-example)
      + [Outputs](#outputs-2)
   * [GLMs](#glms)
      + [Inputs](#inputs-3)
      + [Command line example](#command-line-example-1)
      + [Outputs](#outputs-3)
- [Running on DNANexus](#running-on-dnanexus)
   * [Inputs](#inputs-4)
   * [Outputs](#outputs-4)
   * [Command line example](#command-line-example-2)
      + [Runtime Examples, System Requirements, and Output Expectations](#runtime-examples-system-requirements-and-output-expectations)
      + [Batch Running](#batch-running)

## Introduction

**PLEASE NOTE - ** this applet is very much a work in progress. Please understand that some standard functionalities have
not been built into this applet. I have tried to document various "TO DOs" within this README where applicable, but it is
by no means exhaustive!

This applet performs rare variant burden testing using one of four different methods:

* [BOLT](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)
* [SAIGE-GENE](https://github.com/weizhouUMICH/SAIGE/)
* [STAAR](https://github.com/xihaoli/STAAR)
* GLMs – vanilla linear/logistic models implemented with python's [statsmodels module](https://www.statsmodels.org/stable/index.html)

This README makes use of DNANexus file and project naming conventions. Where applicable, an object available on the DNANexus
platform has a hash ID like:

* file – `file-1234567890ABCDEFGHIJKLMN`
* project – `project-1234567890ABCDEFGHIJKLMN`

Information about files and projects can be queried using the `dx describe` tool native to the DNANexus SDK:

```commandline
dx describe file-1234567890ABCDEFGHIJKLMN
```

**Note:** This README pertains to data included as part of the DNANexus project "MRC - Variant Filtering" (project-G2XK5zjJXk83yZ598Z7BpGPk)

### Background

One of the primary goals of genetics is to determine if genetic variant(s) have a demonstrable effect on a given phenotype.
For common variants (MAF > ~.1%), a GWAS approach is typically used with imputed genotypes from genotyping chips. This 
approach is not as effective as for rare variants (MAF < 0.1%) as the number of individuals that carry a given variant
may be very small and thus we may not have the required power to identify an association. Thus, aggregation approaches
are commonly employed whereby variants of a likely similar consequence or effect are "merged" within or across genes to
perform testing. 

Aggregation can be done in multiple different ways. In this applet, we have implemented four different approaches. We 
also attempted to implement a fifth, REGENIE, but ran into issues with run-time and compute required. Each of the four 
implemented approaches has advantages and disadvantages. I have summarised some primary examples in the sections for 
each tool:

#### 1. [BOLT-LMM](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)

Original Publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4342297/ 

BOLT was originally designed to perform association testing for common variants. BOLT does not, by default, allow for 
rare variant burden testing. Nonetheless, here we have implemented an aggregation approach whereby we create "dummy" variables
to trick BOLT into performing rare variant burden tests. For more information, please see the [section in the methodology](#bolt) 
that covers BOLT in more detail.

Advantages:

* Relatively Fast
* Handles cryptic relatedness

Disadvantages:

* Cannot handle quantification of number of qualifying variants – only has a binary "yes/no" that an individual has a qualifying variant

#### 2. [SAIGE-GENE](https://github.com/weizhouUMICH/SAIGE/)

Original Publication: https://www.medrxiv.org/content/10.1101/2021.07.12.21260400v1

Like BOLT, SAIGE was originally designed to perform association testing for common variants. Unlike BOLT, the authors of
SAIGE have adapted it for running rare variant association tests. Crucially, SAIGE implements the [SKAT-O](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3415556/)
methodology to handle known issues with p. value approximation with extreme case-control imbalance of binary traits which
reduces false-positive rates when performing association testing. SAIGE-GENE performs two tests: 1) it first performs
a test where the total number of qualifying variants (e.g. PTV, missense, etc.) per individual per gene are summed and 2)
a per-variant test to look for an association at all individual qualifying variants. For more information, please see the
[section in the methodology](#saige-gene) that covers SAIGE-GENE in more detail.

Advantages:

* Handles cryptic relatedness
* Adjusts for variance of rare variants
* Natively handles rare variant aggregation

Disadvantages

* Relatively slow
* High memory usage

#### 3. [STAAR](https://github.com/xihaoli/STAAR)

Original Publication: https://www.nature.com/articles/s41588-020-0676-4

STAAR is a rare variant burden testing method that we have implemented to provide an alternate approach to the two methods
listed above. In brief, it uses a similar approach to SAIGE-GENE, whereby qualifying variants are aggregated per 
individual per gene. The primary difference in this approach is that it uses a pre-built set of familial relationships
to handle cryptic relatedness between individuals being tested rather than calculating this relationship during execution.
STAAR also ostensibly contains additional functionality to assign different weights to variants based on a variety of 
different annotations (i.e. methylation, open chromatin, etc.). We have not implemented this as we are 1) primarily dealing
with coding sequence and 2) aggregating variants ourselves. We may implement this model as time or necessity allows. For
more information, please see the [section in the methodology](#staar) that covers STAAR in more detail.

Advantages:

* Very fast

Disadvantages:

* Implemented in R as a package rather than as a command-line application (requires creation of wrapper script)
* Does not natively handle cryptic relatedness

#### 4. Generalised Linear Models (GLMs)

Original Publication: N/A

GLMs are among the most basic approach that can be implemented for rare variant burden testing. In the case of GLMs, we 
have implemented a simple aggregation approach whereby number of qualifying variants per individual per gene are counted
and tested for association with a given phenotype using either a linear or logistic model. The current approach has been
implemented using the [statsmodels module](https://www.statsmodels.org/stable/index.html) for python3. For more information, 
please see the [section in the methodology](#glms) that covers GLMs in more detail.

Advantages:

* Very fast
* The least complex model (less of a black box)

Disadvantages:

* Does not control for case-control imbalance when calculating p. value
* Does not control for cryptic relatedness

### Dependencies

#### Docker

This applet uses [Docker](https://www.docker.com/) to supply dependencies to the underlying AWS instance
launched by DNANexus. The Dockerfile used to build dependencies is available as part of the MRCEpid organisation at:

https://github.com/mrcepid-rap/dockerimages/blob/main/associationtesting.Dockerfile

This Docker image is built off of the primary 20.04 Ubuntu distribution available via [dockerhub](https://hub.docker.com/layers/ubuntu/library/ubuntu/20.04/images/sha256-644e9b64bee38964c4d39b8f9f241b894c00d71a932b5a20e1e8ee8e06ca0fbd?context=explore).
This image is very light-weight and only provides basic OS installation. Other basic software (e.g. wget, make, and gcc) need
to be installed manually. For more details on how to build a Docker image for use on the UKBiobank RAP, please see:

https://github.com/mrcepid-rap#docker-images

In brief, the primary **bioinformatics software** dependencies required by this Applet (and provided in the associated Docker image)
are:

* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [plink1.9](https://www.cog-genomics.org/plink2)
* [plink2](https://www.cog-genomics.org/plink/2.0/)
* [R](https://www.r-project.org/) - v4.1.1 – and specifically the modules:
    * [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
    * [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
    
This Docker image also contains installs for three approaches used in this applet, BOLT, SAIGE, and STAAR. Each of these
approaches have detailed install instructions on their respective websites implemented in the Dockerfile referenced above.

This applet also uses the "exec depends" functionality provided as part of the DNANexus applet building. Please see
[DNANexus documentation on RunSpec](https://documentation.dnanexus.com/developer/apps/execution-environment#external-utilities)
for more information. This functionality allows one to specify python packages the applet needs via pip/apt. Here we install:

* [pandas](https://pandas.pydata.org/)
* [statsmodels](https://www.statsmodels.org/stable/index.html)

See `dxapp.json` for how this is implemented for this applet.

This list is not exhaustive and does not include dependencies of dependencies and software needed
to acquire other resources (e.g. wget). See the referenced Dockerfile for more information.

I have written a custom script (`runSTAAR.R`) for generating a file that we need to run the tool [STAAR](https://github.com/xihaoli/STAAR).
This custom script is placed into the directory:

`resources/usr/bin`

and in accordance with dependency [instructions](https://documentation.dnanexus.com/developer/apps/dependency-management/asset-build-process)
from DNANexus, all resources stored in this folder are then included with the built app at `/usr/bin/` on the launched AWS instance.

#### Resource Files

This applet makes use of the filtered genotype data and genetic relatedness matrices generated using the [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms)
applet included in this repository. Please see that repository for more information.

## Methodology

This applet is step 5 (mrc-runassociationtesting) of the rare variant testing pipeline developed by Eugene Gardner for the UKBiobank
RAP at the MRC Epidemiology Unit:

![](https://github.com/mrcepid-rap/.github/blob/main/images/RAPPipeline.png)

This applet uses the final output of [mrcepid-mergecollapsevariants](https://github.com/mrcepid-rap/mrcepid-mergecollapsevariants.git)
to perform rare variant burden testing. As such, the filtering for variants that the user provided beginning with [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants.git)
will determine what variants are processed by this applet.

This applet proceeds in two basic steps:

1. Covariate processing
2. Rare variant burden testing using selected approaches. For each tool/method, I have documented:
    1. Select inputs generated by [mrcepid-mergecollapsevariants](https://github.com/mrcepid-rap/mrcepid-mergecollapsevariants.git)
    2. An example command-line
    3. Example outputs

The user selects the variant test that they want to run at runtime. In theory, one can run all four approaches using one 
command-line execution, but this may not be feasible due to runtime and compute resources required. Please see the [inputs](#inputs)
and [batch running](#batch-running) sections for more information.

In theory, each tool could have been implemented using individual applets. This was decided against in order to simplify
the inputs provided to each tool. By including covariate and sample filtering within the applet as the first step, we can ensure
that covariates, samples, and phenotypes are run **identically** for all four tools/methods.

Here we have documented in more detail covariate processing and then the implementation of each of the four tools documented
in the [background](#background) section of this README. Note that when we use the term "<file_prefix>", that refers to the
prefix that was provided to mrcepid-collapsevariants and mrcepid-mergecollapsevariants.

### Covariate processing

This applet uses two **TAB-DELIMITED** files for providing phenotypes/covariates:

1. A phenotype to test for an association to rare variants

Please see [inputs](#inputs) for how this file is structured.   

2. Covariates to control for during association testing:

The covariates file is currently hard-coded into the applet. **This will be changed in a future version!** This file has the 
DNANexus hash-id of: `file-G529v8jJxQq1jgxx4gbPpPQj`

The details of this file are not important at the moment but, in brief, contain a raw version of the database output from
the UK Biobank data [implemented on the UKBiobank RAP](https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-phenotypic-data-as-a-file).
Currently, these are (with [showcase](https://biobank.ndph.ox.ac.uk/showcase/) ID):

* Age at assessment - [21003](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21003)
* Genetic sex - [22001](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22001)
* Genetic principal components - [22009](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22009)
* Genotype measurement batch - [22000](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22000)

These covariates are automatically included when performing rare variant burden testing. Sex is only included as a covariate
when performing an association without specifying a sex on the [commandline](#inputs). 

Covariates and samples are processed in the following steps:

1. Generate an exclusion list of individuals NOT to be tested. This is **currently hardcoded** to exclude all related, non-European 
   genetic ancestry individuals. How this file was generated is described in detail as part of [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms#1-selecting-individuals).

2. Subset to individuals retained AFTER filtering genotying data and selecting for individuals who have whole exome sequencing.
   How this file was generated is described in detail as part of [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms#methodology).
   
3. Read the phenotypes file and exclude individuals who have missing (e.g. NaN/NA) data.

4. Ingest the covariates file and restrict to individuals of the requested sex (see [inputs](#inputs)) that were not excluded 
   by steps (1-3). Format covariates into a file suitable for burden testing.
   
5. Create a plink file of the genetic data that only includes individuals that will be used during rare variant burden 
   testing.
   
Individual rare variant burden tests are then run according to [user input](#inputs). 

### BOLT

BOLT roughly proceedes in two steps. First, BOLT computes a genetic relatedness matrix (GRM) from genetic data curated during
the previous step of this applet. Next, it tests for the association a variant provided in a bgen file for an association
with the phenotype of interest. Normally, this bgen file is a collection of phased SNP data from genotying chips. Here,
we have created a "dummy" bgen file that instead encodes presence/absence of a qualifying variant per-individual. This dummy
variable is then used by BOLT to test for an association with a given phenotype.

#### Inputs

A bgen file of genes to test for association with our phenotype of interest. This file is initially formated as a plink 
.ped file like the following:

```text
1000000 1000000 0 0 0 -9 A A A A
1000001 1000001 0 0 0 -9 A C A A
1000002 1000002 0 0 0 -9 A C A A
1000003 1000003 0 0 0 -9 A A A C
1000004 1000004 0 0 0 -9 A A A A
```

The first 6 columns are standard [.ped format](https://www.cog-genomics.org/plink/1.9/formats#ped). The next four columns
are the representation of the dummy variable. Each pair of columns (i.e. 7 & 8, 9 & 10) represent one gene that we are testing. 
Individuals WITH a qualifying variant (e.g. a PTV, missense, etc) are labelled as heterozygous (A C); individuals without a 
qualifying variant are labelled as homozygous (A A). To be clear, this means that individuals 1000001 and 1000002 have a 
qualifying variant in Gene1, and 1000003 has a qualifying variant in Gene2. An accompanying [.map file](https://www.cog-genomics.org/plink/1.9/formats#map)
is provided with this of the format:

```text
chr1 ENSG000000001 0 1000000
chr1 ENSG000000002 0 2000000
```

This file is then converted to .bed using plink2 and then .bgen 1.2 8-bit format using plink2:

```commandline
plink --make-bed --file <file_prefix>.BOLT --out <file_prefix>.BOLT
plink2 --export bgen-1.2 'bits='8 --bfile <file_prefix>.BOLT --out <file_prefix>.BOLT"
```

#### Command Line Example

```commandline
bolt --bfile=UKBB_200K_Autosomes_QCd_WBA /                                       # Filtered genetic data
        --exclude=UKBB_200K_Autosomes_QCd.low_MAC.snplist /                      # List of SNP IDs 
        --phenoFile=phenotypes_covariates.formatted.txt /                        # formated phenotype + covariate file generated during step 1
        --phenoCol=<pheno_name> /                                                # phenotype name extracted from the provided phenotype file
        --covarFile=phenotypes_covariates.formatted.txt /                        # BOLT requires you to provide this twice...
        --covarCol=batch /                                                       # batch is field 22000, we provide here as a categorical variable (covarCol)
        --covarCol=sex /                                                         # sex as a categorical variable
        --qCovarCol=age /                                                        # age as a continuous variable
        --qCovarCol=PC{1:10} /                                                   # First 10 principal components as a continuous variable
        --covarMaxLevels=110 /                                                   # Since there are 110 batches of arrays (batch), we have to tell BOLT to allow for more categorical bins
        --LDscoresFile=BOLT-LMM_v2.3.5/tables/LDSCORE.1000G_EUR.tab.gz /         # LD scores pre-computed by BOLT
        --geneticMapFile=BOLT-LMM_v2.3.5/tables/genetic_map_hg19_withX.txt.gz /  # Genetic map pre-computed by BOLT
        --lmmInfOnly /                                                           # Tell bolt to compute only infinitesimal mixed model association statistics
        --numThreads=32 /                                                        # Number of threads
        --statsFile=<file_prefix>.tarball_prefix.stats.gz /                      # This sets the output file of basic stats
        --verboseStats /                                                         # Give verbose stats
        --bgenFile=bolt_input.bgen /                                             # This is the bgen dummy file we create above
        --sampleFile=bolt_input.sample /                                         # and the associated .sample file
        --statsFileBgenSnps=<file_prefix>.bgen.stats.gz                          # output statistics for testing dummy variables
```

**Note:** BOLT does not differentiate between binary and continuous traits.

#### Outputs

Two files that use the standard BOLT [output](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html#x1-470008):

1. For all SNPs listed in the file provided to --bfile - not particularly useful and provided for posterity (`<file_prefix>.stats.gz`)
2. For all genes listed in the dummy bgen file provided to --bgenFile - the primary output of this applet for BOLT (`<file_prefix>.bgen.stats.gz`)

### SAIGE-GENE

#### Inputs

SAIGE-GENE requires 2 input files that are created during mrcepid-mergecollapsevariants:

1. A VCF format file of rare variants we want to test.

This file was generated during mrcepid-collapsevariants for each VCF file and merged across VCFs during mrcepid-collapsevariants
into a plink .ped format file. This .ped file was then converted into a VCF that contained ONLY variants we want to test using
a command like:

```commandline
plink2 --pfile <file_prefix> --export vcf --out <file_prefix>.SAIGE
```

2. A "groupFile" that tells SAIGE which variants to combine for a given gene.

Each gene we are testing is represented by one tab-delimited row, where the first column is the gene name and subsequent columns
are variants we want to test for that gene. Variant ID is represented as the VCF fields CHROM:POS_REF/ALT and **NOT** by the ID
column.

```text
ENSG000000001 1:1000000_A/T 1:1000010_T/G   1:1000020_G/A
ENSG000000002 1:2000000_G/C 1:2000030_A/C   1:2000050_ATC/A 1:2000000_G/GATC 
```

#### Command Line Example

SAIGE-GENE proceedes in two steps:

1. Fitting the null GLMM (slow):

```commandline
step1_fitNULLGLMM.R
          --plinkFile=UKBB_200K_Autosomes_QCd_WBA /                          # Filtered genetic data
          --phenoFile=phenotypes_covariates.formatted.txt /                  # formated phenotype + covariate file generated during step 1
          --phenoCol=<pheno_name> /                                          # phenotype name extracted from the provided phenotype file
          --isCovariateTransform=FALSE /                                     # Should SAIGE inverse normalise covariates (NO!)
          --sampleIDColinphenoFile=IID /                                     # Sample ID name in our covariates file
          --outputPrefix=SAIGE_OUT /                                         # Output name for step 1
          --outputPrefix_varRatio=SAIGE_OUT_cate /                           # Output name for step 1 of CATE ratios
          --sparseGRMFile=sparseGRM_200K.sparseGRM.mtx /                     # Sparse GRM pre-computed during mrcepid-buildgrms
          --sparseGRMSampleIDFile=sparseGRM_200K.sampleIDs.txt /             # Associated sample file for the GRM
          --nThreads=32 /                                                    # Number of threads
          --LOCO=FALSE /                                                     # Should we do leave-one-chrom-out (NO!)
          --skipModelFitting=FALSE /                                         # Should we skip model fitting and go straight to step2 (NO!)
          --IsSparseKin=TRUE /                                               # Do we provide a sparse GRM (YES!)
          --isCateVarianceRatio=TRUE /                                       # Should be compute CATE variance ratios of rare variants?
          --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age,sex /  # Covariates to control for. Note – Sex is excluded when testing one sex (see source code)
          --traitType=binary/quantitative                                    # set according to type of trait as provided by user
```

**Note:** I have shortened the name of the sparseGRMFile for readability (see source code).

2. Performing rare variant burden tests (fast):

```commandline
step2_SPAtests.R
          --vcfFile=saige_input.vcf.gz \                                # Input vcf file I document above
          --vcfField=GT \                                               # Hardcoded INFO field to check for presence absence in the vcf file
          --GMMATmodelFile=SAIGE_OUT.rda \                              # File generated by step1 above
          --varianceRatioFile=SAIGE_OUT_cate.varianceRatio.txt \        # File generated by step1 above
          --LOCO=FALSE \                                                # Should we do leave-one-chrom-out (NO!)
          --SAIGEOutputFile=<file_prefix>.SAIGE_OUT.SAIGE.gene.txt \    # Output file from this step
          --groupFile=<file_prefix>.SAIGE.groupFile.txt \               # Input groupFile I document above
          --sparseSigmaFile=SAIGE_OUT_cate.sparseSigma.mtx \            # File generated by step1 above
          --IsSingleVarinGroupTest=TRUE \                               # Should we also test individual variants (YES!)
          --MACCutoff_to_CollapseUltraRare=0.5 \                        # Minimum allele count to include a variant. 0.5 means we include all variants (even singletons)
          --IsOutputBETASEinBurdenTest=TRUE                             # Should we output betas and std. errors (YES!)?
```

**Note:** I have shortened the name of the sparseSigmaFile for readability (see source code).

#### Outputs

Two tab-delimited files:

1. Per-gene burden test summary statistics (`<file_prefix>.SAIGE_OUT.SAIGE.gene.txt`).
2. Per-variant burden test summary statistics (`<file_prefix>.SAIGE_OUT.SAIGE.gene.txt_single`)

The headers should be self-explanatory, but please see the [SAIGE documentation](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE#output-files-3)
if you need more information.

### STAAR

STAAR is unique here in that it does not have a command-line tool associated with it and is instead provided as an R package.
To run STAAR, we created an R wrapper script that ingests the data required to run burden tests, and formats the phenotypic
and genetic data into tables / matrices that are compatible with the STAAR package. This commented script is provided in 
this repository at:

`resources/usr/bin/runSTAAR.R`

#### Inputs

STAAR requires the following inputs:

1. A sparse matrix created by the R package 'Matrix'. This file is created in the mrcepid-mergecollapsevariants step of 
   this workflow and saved in the binary .rds format. Specifically, it a lgCMatrix object that stores only information for
   individuals that have given variant as a value of "1":
   
```text
        chr1:1000000:A:T chr1:1000010:T:G chr1:1000020:G:A chr1:2000000:G:C chr1:2000030:A:C
1000000                .                .                .                .                1
1000001                .                1                .                .                .
1000002                .                .                1                .                .
1000003                .                .                .                .                1
1000004                1                .                .                .                . 
```

Column identifiers are _identical_ to that found in the bgen file created for BOLT. Note the differences from the input
[provided to SAIGE-GENE](#saige-gene).

2. A tab-delimited file of variant information:

```text
1:1000000:A:T   1   1000000 1.2e-6
1:1000010:T:G   1   1000010 2.4e-6
1:1000020:G:A   1   1000020 1.2e-6
1:2000000:G:C   1   2000000 1.4e-5
1:2000030:A:C   1   2000030 4.8e-6
```

Column 4 is minor allele frequency.

3. A tab-delimited file of gene annotations:

```text
chr1:1000000:A:T    ENSG000000001   PTV
chr1:1000010:T:G    ENSG000000001   PTV
chr1:1000020:G:A    ENSG000000001   PTV
chr1:2000000:G:C    ENSG000000002   PTV
chr1:2000030:A:C    ENSG000000002   PTV
```

Column 3 is the class of variant being tested.

#### Command line example

This example command-line is to run the script that we have created for this applet. All commands are just space delimited:

```commandline
Rscript runSTAAR.R /                                     
          <file_prefix>.STAAR.matrix.rds \              # File (1) from above
          <file_prefix>.variants_table.STAAR.tsv \      # File (2) from above
          <file_prefix>.REGENIE.annotation \            # File (3) from above – REGENIE is a legacy name and needs to be fixed
          phenotypes_covariates.formatted.txt \         # formated phenotype + covariate file generated during step 1
          <pheno_name> \                                # phenotype name extracted from the provided phenotype file
          <file_prefix> \                               # prefix from mrcepid-mergecollapsevariants
          is_binary                                     # Is this a binary trait? [TRUE/FALSE]
```

Please see the source code for the R script cited above, but in brief, we first fit a null model:

```r
obj_nullmodel <- fit_null_glm(formated.formula, data=data_for_STAAR, family="binomial")
```

and then run a loop that tests the effect of having a rare variant on the given phenotype using STAAR:

```r
for (gene in genes) {
  # We subset the genotype matrix (genotypes; file (1) above) using the rownumbers in our variant file (file (2) above)
  current_GENE <- Matrix(genotypes[,variants[GENEID == gene, rownum]])
  # And then run STAAR using this subset matrix
  staar_result <- STAAR(genotype = current_GENE, obj_nullmodel = obj_nullmodel, rare_maf_cutoff = 1)
}
```

#### Outputs

We output a single tab-delimited file with the following columns (`PTV.STAAR_results.tsv`):

1. geneID - ENSG ID of the gene being tested (e.g. ENSG00000001)
2. n.samps – number of samples run through STAAR()
3. staar.O.p – p. value based on SKAT-O
4. n.var – Number of variants considered for this gene
5. cMAC – Cumulative minor allele count of all variants considered for this gene

### GLMs

The method to perform Generalised Linear Models (GLMs) as part of this applet is implemented using the [statsmodels](https://www.statsmodels.org/stable/index.html)
python package.

#### Inputs

GLMs require 2 files:

1. Variant data is identical to that for [BOLT](#bolt). To ingest this data into Python, we convert it into a sparse matrix
using `bcftools query`:

```commandline
# Convert to bcf file with plink2
plink2 --bgen <file_prefix>.BOLT.bgen 'ref-last' --export bcf --out lm
# Convert to sparse genotype matrix
bcftools query -i "GT='alt'" -f "[%SAMPLE\t%ID\t%GT\n]" lm.bcf > lm.tsv
```

This command creates a tab-delimited file like:

```text
1000000 ENSG000000001 0/1
1000001 ENSG000000001 0/1
1000003 ENSG000000002 0/1
```

2. A list of genes to process:

```text
ENSG000000001
ENSG000000002
```

#### Command line example

This method has no command line as it is run within the applet's source code. For specifics on the method used, please see 
the applet source code. Briefly, this runs in two basic steps:

```python
import pandas as pd
import statsmodels.api as sm

# 1. read file (1) from above as a pandas dataframe:
geno_table = pd.read_csv("lm.tsv",
                         sep = "\t",
                         names = ['eid', 'gene', 'gt'])

# 2. Loop through all genes from file (2):
for gene in genes:
    if is_binary == True:
       family = sm.familes.Binomial()
    else:
       family = sm.families.Gaussian()
       sm.GLM.from_formula(pheno_name + ' ~ has_var + sex + age + batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10', 
                           data=pheno_covars, 
                           family=family).fit()
```

#### Outputs

We output a single tab-delimited file with the following columns (`<file_prefix>.LM_results.tsv`):

1. p_val – raw p. value
2. effect – beta / log OR
3. std_err – standard error
4. n_car – Number of individuals with a qualifying variant in the given gene
1. gene - ENSG ID of the gene being tested (e.g. ENSG00000001)

## Running on DNANexus

### Inputs

|input|description             |
|---- |------------------------|
|association_tarball  | Path/hash ID of the output from [mrcepid-mergecollapsevariants](https://github.com/mrcepid-rap/mrcepid-mergecollapsevariants) that you wish to use for rare variant burden testing. |
|phenofile | Phenotypes file – see below for more information on the format of this file |
|run_bolt    | run BOLT? Default: FALSE |
|run_regenie | run REGENIE? Default: FALSE. **WARNING** Using this option currently does nothing as REGENIE is not properly implemented |
|run_saige   | run SAIGE? Default: FALSE |
|run_staar   | run STAAR? Default: FALSE |
|run_linear_model   | run GLMs? Default: FALSE |
|is_binary   | Is the given trait in the phenofile binary? |
|sex   | Run only one sex or both sexes be run (0 = female, 1 = male, 2 = both) [2]? |
|output_prefix   | Prefix to use for naming output tar file of association statistics. Default is to use the file name 'assoc_stats.tar.gz' |

**Phenotypes File**

A tab-delimited file with three columns and a header. Column 1 and 2 **MUST** be named FID and IID in that order. Column 3 can
have any non-whitespace name (`pheno_name` below) that represents the phenotype that you want to test. The applet automatically 
uses this name when performing variant testing and creating output. Values for binary traits MUST be either 0/1/NA/NaN while
values for continuous traits can be any float or NA/NaN. Individuals with NA/NaN values are automatically excluded during 
testing.

```text
FID IID pheno_name
1000000 1000000 0
1000001 1000001 1
1000002 1000002 1
1000003 1000003 0
1000003 1000003 NA
```

### Outputs

|output                 | description       |
|-----------------------|-------------------|
|output_tarball         |  Output tarball containing   |

output_tarball is either named `assoc_stats.tar.gz` by default. If the parameter `output_prefix` is provided a file like (set `-ioutput_prefix="PTV"`):

`PTV.assoc_stats.tar.gz`

Would be created. This tar.gz file will contain files that are specific to the tool that was requested, but at most can contain 6 files:

1. `<file_prefix>.SAIGE_OUT.SAIGE.gene.txt` (SAIGE-GENE output)
2. `<file_prefix>.SAIGE_OUT.SAIGE.gene.txt` (SAIGE-GENE output)
3. `<file_prefix>.STAAR_results.tsv` (STAAR output)
4. `<file_prefix>.bgen.stats.gz` (BOLT output)
5. `<file_prefix>.stats.gz` (BOLT output)
6. `<file_prefix>.lm_stats.tsv` (GLM output)

Please see each tool's respective output section for more information on these outputs.

### Command line example

If this is your first time running this applet within a project other than "MRC - Variant Filtering", please see our
organisational documentation on how to download and build this app on the DNANexus Research Access Platform:

https://github.com/mrcepid-rap

This applet can be run in a few different ways:

1. Run one tool/method at a time:

```commandline
dx run mrcepid-runassociationtesting --priority low --destination results/ \ 
        -iassociation_tarball=file-G58jv3QJb67q641b4z3kpQpb \
        -iphenofile=file-G5PFv00JXk80Qfx04p9X56X5 \
        -iis_binary=true \
        -ioutput_prefix="PTV" \
        -isex=2 \
        -irun_bolt=true
```

2. Run all tools/methods:

```commandline
dx run mrcepid-runassociationtesting --priority low --destination results/ \ 
        -iassociation_tarball=file-G58jv3QJb67q641b4z3kpQpb \
        -iphenofile=file-G5PFv00JXk80Qfx04p9X56X5 \
        -iis_binary=true \
        -ioutput_prefix="PTV" \
        -isex=2 \
        -irun_bolt=true \
        -irun_saige=true \
        -irun_staar=true \
        -irun_linear_model=true
```

**Note:** boolean parameters must be provided in json-compatible case ('true' or 'false').

Brief I/O information can also be retrieved on the command line:

```commandline
dx run mrcepid-runassociationtesting --help
```

I have set a sensible (and tested) default for compute resources on DNANexus for running BOLT and SAIGE. This is baked into the json used for building
the app (at `dxapp.json`) so setting an instance type when running either of those tools or all tools at once is unnecessary.
This current default is for a mem1_ssd1_v2_x36 instance (36 CPUs, 72 Gb RAM, 900Gb storage). If running either STAAR or GLMs
it is recommended to change instance type to save money. Do this by providing the flag `--instance-type mem1_ssd1_v2_x8`.

#### Runtime Examples, System Requirements, and Output Expectations

| Tool | Per-gene p. value | Per-variant p. value | Accurate β / OR | Runtime§ | Cost§ |
| ---- | ----------------- | -------------------- | --------------- | ------- | ---- |
| BOLT | <span style="color:green">**TRUE**</span> | <span style="color:red">**FALSE**</span> | <span style="color:green">**TRUE**</span> / <span style="color:red">**FALSE**</span> | 2.5hrs | £0.58 |
| SAIGE | <span style="color:green">**TRUE**</span> | <span style="color:green">**TRUE**</span> | <span style="color:green">**TRUE**</span> / <span style="color:red">**FALSE**</span> | 20hrs | £4.67 |
| STAAR | <span style="color:green">**TRUE**</span> | <span style="color:red">**FALSE**</span> | <span style="color:red">**FALSE**</span> / <span style="color:red">**FALSE**</span> | 0.1hrs | £0.01‡ |
| GLM | <span style="color:green">**TRUE**</span> | <span style="color:red">**FALSE**</span> | <span style="color:green">**TRUE**</span> / <span style="color:green">**TRUE**</span> | 0.1hrs | £0.01‡ |

§Numbers are representative runtimes and cost for a real analysis using two VCF slices on the effect of rare variant burden on having type II diabetes. 

‡Cost for STAAR and GLM is using smaller instances (`mem1_ssd1_v2_x8`).

#### Batch Running

This applet is not compatible with batch running.

