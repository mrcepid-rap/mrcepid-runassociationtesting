# RunAssociationTesting (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
  * [Changelog](#changelog)
  * [Background](#background)
    + [Burden Tools Implemented](#burden-tools-implemented)
      - [1. BOLT-LMM](#1-bolt-lmmhttpsalkesgroupbroadinstituteorgbolt-lmmbolt-lmm_manualhtml)
      - [2. SAIGE-GENE+](#2-saige-genehttpsgithubcomsaigegitsaige)
      - [3. STAAR](#3-staarhttpsgithubcomxihaolistaar)
      - [4. REGENIE](#4-regeniehttpsrgcgithubgithubioregenie)
      - [5. Generalised Linear Models (GLMs)](#5-generalised-linear-models-glms)
  * [Dependencies](#dependencies)
    + [Docker](#docker)
    + [Resource Files](#resource-files)
- [Methodology](#methodology)
  * [Covariate processing](#covariate-processing)
  * [Burden Tests](#burden-tests)
    + [BOLT](#bolt)
      - [Inputs](#inputs)
      - [Command Line Example](#command-line-example)
    + [SAIGE-GENE+](#saige-gene-)
      - [Inputs](#inputs-1)
      - [Command Line Example](#command-line-example-1)
    + [STAAR](#staar)
      - [Inputs](#inputs-2)
      - [Command line example](#command-line-example)
    + [REGENIE](#regenie)
      - [Inputs](#inputs-3)
      - [Command line example](#command-line-example-1)
    + [GLMs](#glms)
      - [Inputs](#inputs-4)
      - [Command line example](#command-line-example-2)
  * [Variant Extraction](#variant-extraction)
  * [PheWAS](#phewas)
- [Running on DNANexus](#running-on-dnanexus)
  * [Inputs](#inputs-5)
    + [Mode](#mode)
    + [Association Tarballs](#association-tarballs)
    + [Phenotypes File](#phenotypes-file)
    + [Inclusion / Exclusion lists](#inclusion--exclusion-lists)
    + [Additional Covariate (Quantitative / Categorical) File](#additional-covariate--quantitative---categorical--file)
    + [Gene IDs](#gene-ids)
  * [Outputs](#outputs)
    + [Per-gene output](#per-gene-output)
    + [Extract and PheWAS mode outputs](#extract-and-phewas-mode-outputs)
    + [Per-marker output](#per-marker-output)
  * [Command line example](#command-line-example-3)
    + [Selecting an Instance Type](#selecting-an-instance-type)
    + [Batch Running](#batch-running)

## Introduction

**PLEASE NOTE - ** this applet is very much a work in progress. Understand that some standard functionalities within 
individual tools used by this applet have not yet been implemented. I have tried to document various "TO DOs" within this
README where applicable, but it is by no means exhaustive!

This applet performs rare variant burden testing and associated functions (e.g. extracting variants/samples and various 
ad hoc tests). Rare variant burden tests are implemented with one of four different methods:

* [BOLT](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)
* [SAIGE-GENE+](https://github.com/saigegit/SAIGE)
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

### Changelog

* v.1.2.4
  * Made filename updates to the codebase to reflect release of the 470k dataset
  * Code which loads VEP annotations was identical in multiple locations. Switched to a single function in association_resources.py

* v1.2.3
  * Fixed a bug related to parsing of .tar.gz files provided from collapsevariants 

* v1.2.2
  * Added language to the README to support running in gene-list mode (see mrcepid-collapsevariants for more details)
  * Additional bug-fix updates for v1.2.0
    * Added better error message for when a gene symbol is not found
  * Fixed a bug in GLMs where iteration over gene lists was not happening properly

* v1.2.1
  * Bug-fix update for v1.2.0
  * Adds improved output formatting in PheWAS mode

* v1.2.0
  * Implemented REGENIE
    * Please see the documentation below for how REGENIE is implemented in this app
  * The extract and phewas modes now run a STAAR model in addition to a standard GLM. Per-gene p. values from the various STAAR models are appended to the end of the table provided as output.
  * We have implemented a bug fix for STAAR which will exclude/include the GRM depending on issues with case/control imbalance for binary traits
    * Traits which fail during fitting of the initial null model will be run WITHOUT a GRM. This is now indicated in the STAAR output via the 'relatedness.correction' column as a boolean TRUE/FALSE value.
  * SAIGE-GENE per-marker tests now use multiple threads (4) to try and speed-up runtime.
  * Behind-the-scenes modification of the code-base to prepare for gene-set testing
  * Several behind-the-scenes modifications to improve readability of the code-base

* v1.1.1
  * Modified how STAAR handles the GRM.
    * There were issues with individuals that had pair-wise relatedness < 0.05. Individuals with relatedness < 0.05 are now thresholded to 0.05.
    * Refer to [this issue](https://github.com/xihaoli/STAAR/issues/11) on GitHub for our discussions with the authors of STAAR.
    * This should not have a major effect as the min relatedness coefficient was 0.0442 

* v1.1.0
  * Updated SAIGE to v1.0.4
    * Just adds support for default SAMPLE files
    * Also fixed minor issues with the SAIGE runtime
    * Have now implemented the GRM provided by Davies as the primary GRM for SAIGE and STAAR
  * Modified how GLMs are run
    * GLMs now run 'variant counts' as a continuous variable. i.e. hets are coded as '1' and homs are coded as '2'
    * This was done to prepare for multi-gene testing in the future
  * When run in 'extract' or 'phewas' mode, GLMs will now always report std errors and p. values regardless of initial p. value  
  * The PheWAS mode has been overhauled and should be much faster now
  * Multiple phenofiles can be provided at once, so long as they are run using identical model-familes (i.e. binary/or continuous). See this readme for how to specify multiple phenofiles.
  * Changed the default GRM to be identical to that generated by [Bycroft et al.](https://www.nature.com/articles/s41586-018-0579-z)

* v1.0.0 – Initial numbered release.
  * Made code object oriented for easier maintainability
  * Added new ‘tool’ option, replaces individual run_xyz flags
  * Added new ‘mode’ option, current modes are ‘burden’ and ‘extract’. Will add more later.
  * ExtractVariants is now deprecated. All functionality has been added to ‘runassociationtesting’
  * Added modified ‘gene_ids’ option. In addition to taking a single-gene, this will also take a comma-separated list of SYMBOLS/ENSTs to process multiple genes at once
  * Changed how the sparse GRM is built – now using default UKBB
  * -iassociation_tarballs will now take either a LIST of tarballs, or a single tar file
  * default output is now assoc_results.tar.gz
  * added ‘phenoname’ field. Not required, but can now use a phenotype file with more than 3 columns by providing the name of the requested column
  * Added new ‘phewas’ mode.
      * providing `phenoname` throws an error
  * Updated SAIGE to v1.0.1
      * Implemented faster per-marker tests for SAIGE (similar to BOLT now)
      * No longer use variance Ratio Estimation
      * SAIGE works well with updated GRM (i.e. good p. values)
  * Added a per-SNP test mode

### Background

One of the primary goals of genetics is to determine if genetic variant(s) have a demonstrable effect on a given phenotype.
For common variants (MAF > ~.1%), a GWAS approach is typically used with imputed genotypes from genotyping chips. This 
approach is not as effective for rare variants (MAF < 0.1%) as the number of individuals that carry a given variant
may be very small; thus we may not have the required power to identify an association. Aggregation approaches
are commonly employed whereby variants of a likely similar consequence or effect are "merged" within or across genes to
perform testing. This code performs several functions via a 'mode' flag. Please see the [command-line examples](#command-line-example)
below for more information on how this works. In brief, this applet has three different modes:

1. Burden Testing via various approaches (`burden`)
    * See [Burden Tools Implemented](#burden-tools-implemented) for more information
2. Variant Extraction (`extract`)
    * Extracts variants/samples/phenotypes identically to how set-up is performed for burden testing in point (1.) above.
3. PheWAS (`phewas`)
    * Uses a GLM as [described below](#4-generalised-linear-models-glms) to run multiple gene-phenotype combinations at once

#### Burden Tools Implemented

Aggregation can be done in multiple different ways. In this applet, we have implemented four different approaches. We 
also attempted to implement a fifth, REGENIE, but ran into issues with run-time and compute required. Implementation of 
REGENIE is in active development. Each of the four implemented approaches has advantages and disadvantages. I have 
summarised some primary examples in the sections below for each tool.

##### 1. [BOLT-LMM](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)

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

##### 2. [SAIGE-GENE+](https://github.com/saigegit/SAIGE)

Original Publication: https://www.medrxiv.org/content/10.1101/2021.07.12.21260400v2

Like BOLT, SAIGE was originally designed to perform association testing for common variants. Unlike BOLT, the authors of
SAIGE have adapted it for running rare variant association tests in the form of SAIGE-GENE+. Crucially, SAIGE-GENE+ implements 
the [SKAT-O](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3415556/) methodology to handle known issues with p. value
approximation with extreme case-control imbalance of binary traits which reduces false-positive rates when performing
association testing. SAIGE-GENE+ performs two tests: 1) it first performs a test where the total number of qualifying
variants (e.g. PTV, missense, etc.) per individual per gene are summed and 2) a per-variant test to look for an
association at all individual qualifying variants. For more information, please see the [section in the methodology](#saige-gene) 
that covers SAIGE-GENE+ in more detail.

Advantages:

* Handles cryptic relatedness
* Adjusts for variance of rare variants
* Natively handles rare variant aggregation

Disadvantages

* Relatively slow
* High memory usage

##### 3. [STAAR](https://github.com/xihaoli/STAAR)

Original Publication: https://www.nature.com/articles/s41588-020-0676-4

STAAR is a rare variant burden testing method that we have implemented to provide an alternate approach to the two methods
listed above. In brief, it uses a similar approach to SAIGE-GENE+, whereby qualifying variants are aggregated per 
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

##### 4. [REGENIE](https://rgcgithub.github.io/regenie/)

Original Publication: https://www.nature.com/articles/s41588-021-00870-7

REGENIE is very similar to, and has a similar history to, BOLT and SAIGE. Originally developed for GWAS, it now has 
methodology for performing rare variant burden tests.

Advantages:

* Controls for relatedness
* Runs gene-burden tests by default
* Implements several different statistical models (e.g. SKAT, SKAT-O, ACAT, ACAT-O)

Disadvantages:

* Memory, time, and CPU intensive

##### 5. Generalised Linear Models (GLMs)

Original Publication: N/A

GLMs are among the most basic approach that can be implemented for rare variant burden testing. In the case of GLMs, we 
have implemented a simple aggregation approach whereby number of qualifying variants per individual per gene are counted
and tested for association with a given phenotype using either a linear or logistic model. The current approach has been
implemented using the [statsmodels module](https://www.statsmodels.org/stable/index.html) for python3. For more information, 
please see the [section in the methodology](#glms) that covers GLMs in more detail.

Advantages:

* Relatively fast
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
* [qctool](https://www.well.ox.ac.uk/~gav/qctool_v2/index.html)
* [bgenix](https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md)
* [R](https://www.r-project.org/) - v4.1.1 – and specifically the modules:
    * [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
    * [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
    
This Docker image also contains executables for four approaches used in this applet, BOLT, SAIGE, STAAR, and REGENIE. Each of these
approaches have detailed install instructions on their respective websites implemented in the Dockerfile referenced above.

This applet also uses the "exec depends" functionality provided as part of the DNANexus applet building. Please see
[DNANexus documentation on RunSpec](https://documentation.dnanexus.com/developer/apps/execution-environment#external-utilities)
for more information. This functionality allows one to specify python packages the applet needs via pip/apt. Here we install:

* [pandas](https://pandas.pydata.org/)
* [statsmodels](https://www.statsmodels.org/stable/index.html)

See `dxapp.json` for how this is implemented for this applet.

This list is not exhaustive and does not include dependencies of dependencies and software needed
to acquire other resources (e.g. wget). See the referenced Dockerfile for more information.

I have written two custom scripts (`runSTAAR_Null.R` and `runSTAAR_Genes.R`) for generating a file that we need to run 
the tool [STAAR](https://github.com/xihaoli/STAAR). These custom scripts are placed into the directory:

`resources/usr/bin`

and in accordance with dependency [instructions](https://documentation.dnanexus.com/developer/apps/dependency-management/asset-build-process)
from DNANexus, all resources stored in this folder are then included with the built app at `/usr/bin/` on the launched AWS instance.

#### Resource Files

This applet makes use of the filtered genotype data and genetic relatedness matrices generated using the [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms)
applet included in this repository. Please see that repository for more information.

## Methodology

This applet is step 5 (mrc-runassociationtesting) of the rare variant testing pipeline developed by Eugene Gardner for the UKBiobank
RAP at the MRC Epidemiology Unit:

![](https://github.com/mrcepid-rap/.github/blob/main/images/RAPPipeline.v3.png)

This applet uses the final output of [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants.git)
to perform rare variant burden testing. As such, the filtering for variants that the user provided beginning with [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants.git)
will determine what variants are processed by this applet.

This applet proceeds in two basic steps:

1. Covariate/phenotype processing and individual-level filtering
2. Running various modules 
   1. Rare variant burden testing using selected approaches (`burden`). For each tool/method, I have documented:
       1. Select inputs generated by [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants.git)
       2. An example command-line
       3. Example outputs
   2. Variant Extraction (`extract`)
   3. PheWAS (`phewas`)

The user selects the variant test/mode that they want to run at runtime. Please see the [inputs](#inputs)
and [batch running](#batch-running) sections for more information.

### Covariate processing

This applet uses at least two **TAB-DELIMITED** files for providing phenotypes/covariates:

1. A phenotype to test for an association to rare variants. This file can have multiple phenotypes. The specific phenotype 
   to test can be specified via the command-line or all phenotypes can be tested if running in 'PheWAS' mode. 

Please see [inputs](#inputs) for how this file is structured.

2. Standard covariates to control for during association testing:

A default set of standard covariates is current used by this applet applet. This file has the DNANexus hash-id of: `file-G8ZvZFQJJv8qZQ9P12KVjYjy`

In brief, this file contain a raw version of the database output from the UK Biobank data 
[implemented on the UKBiobank RAP](https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-phenotypic-data-as-a-file).
Currently, these are (with [showcase](https://biobank.ndph.ox.ac.uk/showcase/) ID if applicable):

* Age at assessment - [21003](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21003)
* Genetic sex - [22001](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22001)
* The First 10 genetic principal components - [22009.1-10](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22009)
* WES Batch - Generated by pulling .fam files for the respective UKBB WES batches and extracting individuals unique to each fam file

These covariates are automatically included when performing rare variant burden testing. Sex is only included as a covariate
when performing an association without specifying a sex on the [commandline](#inputs). 

3. Additional quantitative and categorical covariates. Please see the [inputs](#inputs) section for more information on 
   how these files should be formatted.

Covariates and samples are processed in the following steps:

1. Generate an exclusion list of individuals NOT to be tested by limiting to/keeping samples in the files provided to 
   input parameters `inclusion_list`/`exclusion_list`. Three such files are provided as part of this prject. How these files
   were generated is described in detail as part of [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms#1-selecting-individuals).

2. Subset to individuals retained AFTER filtering genotying data and selecting for individuals who have whole exome sequencing.
   How this file was generated is described in detail as part of [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms#methodology).
   
3. Read the phenotypes file and exclude individuals who have missing (e.g. NaN/NA) data.

4. Ingest the covariates file(s) and restrict to individuals of the requested sex (see [inputs](#inputs)) that were not excluded 
   by steps (1-3). Format covariates into a file suitable for burden testing.
   
5. Create a plink file of the genetic data that only includes individuals that will be used during rare variant burden 
   testing.
   
One of three individual tools are then run according to [user input](#inputs):

1. Burden tests as implemented in BOLT, SAIGE, STAAR, or a GLM
2. Variant/sample extraction
3. PheWAS on a set of gene/phenotype combinations

Methodology for each is given in the following sections.

### Burden Tests

In theory, each tool could have been implemented using individual applets. This was decided against in order to simplify
the inputs provided to each tool. By including covariate and sample filtering within the applet as the first step, we can ensure
that covariates, samples, and phenotypes are run **identically** for all four tools/methods.

Here we have documented in more detail covariate processing and then the implementation of each of the four tools documented
in the [background](#background) section of this README. Note that when we use the term "<file_prefix>", that refers to the
prefix that was provided to mrcepid-collapsevariants and mrcepid-mergecollapsevariants.

#### BOLT

BOLT roughly proceeds in two steps. First, BOLT computes a genetic relatedness matrix (GRM) from genetic data curated during
the previous step of this applet. Next, it tests for the association a variant provided in a bgen file for an association
with the phenotype of interest. Normally, this bgen file is a collection of phased SNP data from genotying chips. Here,
we have created a "dummy" bgen file that instead encodes presence/absence of a qualifying variant per-individual. This dummy
variable is then used by BOLT to test for an association with a given phenotype. We test for both per-gene and per-marker 
(i.e. SNV/InDels) association with our phenotype of interest.

##### Inputs

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
chr1 ENST000000001 0 1000000
chr1 ENST000000002 0 2000000
```

This file is then converted to .bed using plink2 and then .bgen 1.2 8-bit, ref-last format using plink2:

```commandline
plink --make-bed --file <file_prefix>.BOLT --out <file_prefix>.BOLT
plink2 --export bgen-1.2 'bits='8 --bfile <file_prefix>.BOLT --out <file_prefix>.BOLT"
```

The above steps are done per-chromosome and provided via the `--bgenSamplesFileList` argument as described below.

To perform per-marker tests, we simply convert the VCF file used for [SAIGE-GENE+](#saige-gene) into bgen format using plink2, 
and run exactly the same command as described below. Inputs and outputs are essentially identical. 

##### Command Line Example

```commandline
bolt --bfile=UKBB_200K_Autosomes_QCd_WBA /                                       # Filtered genetic data
        --exclude=UKBB_200K_Autosomes_QCd.low_MAC.snplist /                      # List of SNP IDs 
        --phenoFile=phenotypes_covariates.formatted.txt /                        # formated phenotype + covariate file generated during step 1
        --phenoCol=<pheno_name> /                                                # phenotype name extracted from the provided phenotype file
        --covarFile=phenotypes_covariates.formatted.txt /                        # BOLT requires you to provide this twice...
        --covarCol=batch /                                                       # batch is field 22000, we provide here as a categorical variable (covarCol)
        --covarCol=sex /                                                         # sex as a categorical variable
        --covarCol=wes_batch /                                                   # Batch a given individual's WES is derived from (50k / 200k / 450k) 
        --qCovarCol=age /                                                        # age as a continuous variable
        --qCovarCol=PC{1:10} /                                                   # First 10 principal components as a continuous variable
        --covarMaxLevels=110 /                                                   # Since there are 110 batches of arrays (batch), we have to tell BOLT to allow for more categorical bins
        --LDscoresFile=BOLT-LMM_v2.3.5/tables/LDSCORE.1000G_EUR.tab.gz /         # LD scores pre-computed by BOLT
        --geneticMapFile=BOLT-LMM_v2.3.5/tables/genetic_map_hg19_withX.txt.gz /  # Genetic map pre-computed by BOLT
        --lmmInfOnly /                                                           # Tell bolt to compute only infinitesimal mixed model association statistics
        --numThreads=32 /                                                        # Number of threads
        --statsFile=<file_prefix>.tarball_prefix.stats.gz /                      # This sets the output file of basic stats
        --verboseStats /                                                         # Give verbose stats
        --bgenSampleFileList=bolt_input.bgen /                                   # This is a list (one file for each chromosome) of dummy bgen files created above
        --statsFileBgenSnps=<file_prefix>.bgen.stats.gz                          # output statistics for testing dummy variables
```

**Note:** BOLT does not differentiate between binary and continuous traits.

#### SAIGE-GENE+

##### Inputs

SAIGE-GENE+ requires 2 input files that are created during mrcepid-mergecollapsevariants:

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
ENST000000001   1:1000000_A/T   1:1000010_T/G   1:1000020_G/A
ENST000000001   foo foo foo
ENST000000002   1:2000000_G/C   1:2000030_A/C   1:2000050_ATC/A 1:2000000_G/GATC
ENST000000002   foo foo foo foo 
```

Note the interleaved row for each gene. This is used by SAIGE to define a mask for each variant. We create a dummy value
so we can define our own masks.

Files are created per-chromosome to enable faster parallelization on DNA Nexus.

##### Command Line Example

SAIGE-GENE+ proceedes in two steps:

1. Fitting the null GLMM (done for the whole genome):

```commandline
step1_fitNULLGLMM.R
          --phenoFile=phenotypes_covariates.formatted.txt /                             # formated phenotype + covariate file generated during step 1
          --phenoCol=<pheno_name> /                                                     # phenotype name extracted from the provided phenotype file
          --isCovariateTransform=FALSE /                                                # Should SAIGE inverse normalise covariates (NO!)
          --sampleIDColinphenoFile=IID /                                                # Sample ID name in our covariates file
          --outputPrefix=<phenoname>.SAIGE_OUT /                                        # Output name for step 1
          --sparseGRMFile=sparseGRM_200K.sparseGRM.mtx /                                # Sparse GRM pre-computed during mrcepid-buildgrms
          --sparseGRMSampleIDFile=sparseGRM_200K.sampleIDs.txt /                        # Associated sample file for the GRM
          --nThreads=64 /                                                               # Number of threads (default 64)
          --LOCO=FALSE /                                                                # Should we do leave-one-chrom-out (NO!)
          --skipModelFitting=FALSE /                                                    # Should we skip model fitting and go straight to step2 (NO!)
          --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age,sex,wes_batch /   # Covariates to control for. Note – Sex is excluded when testing one sex (see source code)
          --qCovarColList=wes_batch                                                     # Categorical covariates to control for from the overal covarColList
          --traitType=binary/quantitative \                                             # set according to type of trait as provided by user
          --useSparseGRMtoFitNULL=TRUE                                                  # Use the sparse GRM to fit the null model rather than estimating the variance ratio from genetic data (massively improves runtime)
```     

**Note:** I have shortened the name of the sparseGRMFile for readability (see source code).

2. Performing rare variant burden tests (done per chromosome):

```commandline
step2_SPAtests.R
          --vcfFile=saige_input.vcf.gz \                                          # Input vcf file I document above
          --vcfField=GT \                                                         # Hardcoded INFO field to check for presence absence in the vcf file
          --GMMATmodelFile=SAIGE_OUT.rda \                                        # File generated by step1 above
          --sparseGRMFile=sparseGRM_200K.sparseGRM.mtx /                          # Sparse GRM pre-computed during mrcepid-buildgrms
          --sparseGRMSampleIDFile=sparseGRM_200K.sampleIDs.txt /                  # Associated sample file for the GRM
          --LOCO=FALSE \                                                          # Should we do leave-one-chrom-out (NO!)
          --SAIGEOutputFile=<file_prefix>.<chromosome>.SAIGE_OUT.SAIGE.gene.txt \ # Output file from this step
          --groupFile=<file_prefix>.<chromosome>.SAIGE.groupFile.txt \            # Input groupFile I document above
          --is_output_moreDetails=TRUE                                            # Output additional information in the output file (het/hom counts, carrier status, etc.)
          --maxMAF_in_groupTest=0.5 \                                             # Minimum allele count to include a variant. 0.5 means we include all variants (even singletons)
          --maxMissing=1                                                          # We define our own missingness, so set a value to ignore
          --chrom=<chromosome>                                                    # Dummy chromosome value
          --annotation_in_groupTest=foo                                           # Dummy value from the groupFile
```

**Note:** I have shortened the name of the sparseSigmaFile for readability (see source code).

#### STAAR

STAAR is unique here in that it does not have a command-line tool associated with it and is instead provided as an R package.
To run STAAR, we created two R wrapper scripts that 

1. Generates a null model (`runSTAAR_Null.R`).
2. Ingests the data required to run burden tests, and formats the phenotypic and genetic data into tables / matrices that are compatible with
the STAAR package (`runSTAAR_Genes.R`).
   
These commented scripts are provided in this repository at:

`resources/usr/bin/`

##### Inputs

STAAR requires the following inputs:

1. A sparse matrix created by the R package 'Matrix'. This file is created in the mrcepid-collapsevariants step of 
   this workflow and saved in the binary .rds format. Specifically, it a lgCMatrix object that stores only information for
   individuals that have given variant as a value of "1" or "2" depending on het/hom status.
   
```text
        chr1:1000000:A:T chr1:1000010:T:G chr1:1000020:G:A chr1:2000000:G:C chr1:2000030:A:C
1000000                .                .                .                .                1
1000001                .                2                .                .                .
1000002                .                .                1                .                .
1000003                .                .                .                .                1
1000004                1                .                .                2                . 
```

Column identifiers are _identical_ to that found in the bgen file created for BOLT. Note the differences from the input
[provided to SAIGE-GENE+](#saige-gene).

2. A tab-delimited file of variant information with a header. varID is identical to that stored for SAIGE:

```text
varID   chrom   pos ENST    column
1:1000000:A:T   1   1000000 ENST00000000001 1
1:1000010:T:G   1   1000010 ENST00000000001 2
1:1000020:G:A   1   1000020 ENST00000000001 3
1:2000000:G:C   1   2000000 ENST00000000002 4
1:2000030:A:C   1   2000030 ENST00000000002 5
```

3. A sparse GRM. 

Currently, we use the same GRM created for SAIGE as part of the [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms). 
This input is hard-coded into the `runSTAAR.R` script.

##### Command line example

This example command-line is to run both scripts that we have created for this applet for one chromosome. All commands are
space delimited:

```commandline
Rscript runSTAAR_Null.R \
          phenotypes_covariates.formatted.txt \               # Formated phenotypes/covariates to derive the NULL STAAR model
          <pheno_name> \                                      # phenotype name extracted from the provided phenotype file
          <is_binary> \                                       # Is the trait binary (false for quantitative)?
          <quant_covars> \                                    # Names of additional quantitative covariates to include in model (NULL for none)
          <catagorical_covars>                               # Names of additional catagorical covariates to include in model (NULL for none)

Rscript runSTAAR_Genes.R \                                     
          <file_prefix>.<chr>.STAAR.matrix.rds \              # File (1) from above
          <file_prefix>.<chr>.variants_table.STAAR.tsv \      # File (2) from above
          <pheno_name>.STAAR_null.rds                         # Output null model from runSTAAR_Null.R
          <pheno_name> \                                      # phenotype name extracted from the provided phenotype file
          <file_prefix> \                                     # prefix from mrcepid-collapsevariants
          <chr>                                               # Current chromosome we are running to name output files
```

Please see the source code for the R scripts cited above, but in brief, we first fit a null model (in `runSTAAR_Null.R)`:

```r
obj_nullmodel <- fit_null_glmmkin(formated.formula, data=data_for_STAAR, id="FID", family=binomial(link="logit"), kins = sparse_kinship)
```

and then run a loop that tests the effect of having a rare variant on the given phenotype using STAAR (in `runSTAAR_Genes.R)`:

```r
for (gene in genes) {
  # We subset the genotype matrix (genotypes; file (1) above) using the rownumbers in our variant file (file (2) above)
  current_GENE <- Matrix(genotypes[,variants[GENEID == gene, rownum]])
  # And then run STAAR using this subset matrix
  staar_result <- STAAR(genotype = current_GENE, obj_nullmodel = obj_nullmodel, rare_maf_cutoff = 1)
}
```

#### REGENIE

##### Inputs

Like BOLT, REGENIE requires the genetic data to first estimate a null model and control for relatedness between study 
participants. It then uses the derived null model to perform rare variant burden tests individually on masks/chromosome
combinations. Specific inputs for REGENIE are derived from the filtered WES data listed in `project-G6BJF50JJv8p4PjGB9yy7YQ2:file-G86GJ3jJJv8fbXVB9PQ2pjz6`.
Specifically, for burden tests, REGENIE requires:

1. A bgen format file of variants. Unlike for the other tools, we do not create this file during mrcepid-collapsevariants
and instead directly use the filtered bgen files listed in  `project-G6BJF50JJv8p4PjGB9yy7YQ2:file-G86GJ3jJJv8fbXVB9PQ2pjz6`,
filtered to individuals in accordance with inclusion/exclusion lists. This is due to how REGENIE handles variant masks.

2. An annotation file. This file lists the variant ID, identical to that provided to SAIGE-GENE, a gene, and the mask name
provided by the collapsevariants tarball. No header is required:

```text
1:1000000:A:T   ENST00000000001 HC_PTV-MAF_01
1:1000010:T:G   ENST00000000001 HC_PTV-MAF_01
1:1000020:G:A   ENST00000000001 HC_PTV-MAF_01
1:2000000:G:C   ENST00000000002 HC_PTV-MAF_01
1:2000030:A:C   ENST00000000002 HC_PTV-MAF_01
```

3. A gene-to-variant file. This file is almost identical to that provided to SAIGE-GENE, except with additional columns for
chromosome and gene start coordinate, and the list of variant IDs in comma-separated formated:

```text
ENST000000001 1 12345   1:1000000:A:T,1:1000010:T:G,1:1000020:G:A
ENST000000002 2 67890   1:2000000:G:C,1:2000030:A:C,1:2000050:ATC:A,1:2000000:G:GATC 
```

4. A file with one row that provides the name of the mask as listed in file (2) from above. 

```text
HC_PTV-MAF_01   HC_PTV-MAF_01
```

##### Command line example

Step One of REGENIE uses a command like:

```commandline
regenie  \                                                                  
  --step 1  \                                                               # Indicates to REGENIE to run step 1 of 2
  --bed /test/genetics/UKBB_450K_Autosomes_QCd_WBA  \                       # UKBB genetic data filtered to individuals analysed in this specific run
  --covarFile /test/phenotypes_covariates.formatted.txt  \                  # The processed and formated covariate/pheno file from above
  --phenoFile /test/phenotypes_covariates.formatted.txt  \                  # The processed and formated covariate/pheno file from above (required to be listed twice by REGENIE)
  --extract /test/REGENIE_extract.snplist  \                                # A set of variants with MAC > 100 & < ((N_Indv*2) - 100). The latter is because plink does not appear to calculate AC based on the minor allele, but rather on the alternate allele
  --bsize 100  \                                                            # REGENIE computation parameter
  --out /test/fit_out  \                                                    # Name of the file that contains the null model for REGENIE
  --threads 64 \                                                            # Number of threads. Default is to use 64, will be modified by the instance type selected
  --phenoCol <phenoname> \                                                  # Name of the phenotype in covarFile/phenoFile
  --covarColList PC{1:10},age,age_squared,sex \                             # Quantitative covariates to include. Only standard default covariates are listed here.
  --catCovarList wes_batch \                                                # Categorical covariates to include. Only standard default covariates are listed here.
  --bt                                                                      # Indicates running a binary trait. Only included for binary phenotypes
```

Step Two of REGENIE uses a command like:

```commandline
regenie  \                                                                                  
  --step 2  \                                                                               # Indicates to REGENIE to run step 2 of 2    
  --bgen /test/<chromosome>.markers.bgen  \                                                 # QCd bgen file                                    
  --sample /test/<chromosome>.markers.bolt.sample  \                                        # Matching .sample file for file provided by .bgen                                                 
  --covarFile /test/phenotypes_covariates.formatted.txt  \                                  # File identical to that provided to step 1 above                                               
  --phenoFile /test/phenotypes_covariates.formatted.txt  \                                  # File identical to that provided to step 1 above                                                                                                            
  --firth --approx  \                                                                       # Use Firth regression to calculate p. values with the approximation speedup (--approx)            
  --pred /test/fit_out_pred.list  \                                                         # Output file from step 1                        
  --anno-file /test/<tarball_prefix>.<chromosome>.REGENIE.annotationFile.tsv  \             # Annotations file (see inputs section above)                                                                            
  --set-list /test/<tarball_prefix>.<chromosome>.REGENIE.setListFile.tsv  \                 # Set list file (see inputs section above)                                                                       
  --mask-def /test/<tarball_prefix>.<chromosome>.REGENIE.maskfile.tsv  \                    # Mask definition file (see inputs section above)                                                                        
  --aaf-bins 1  \                                                                           # Tells REGENIE to include ALL variants (MAF < 100%) as we define the MAF bin ourselves       
  --vc-tests skato-acat,acato-full  \                                                       # Provide p. values for skat-o and acato-full (see REGENIE full documentation for more information)                       
  --bsize 200  \                                                                            # REGENIE computation parameter
  --threads 1  \                                                                            # Run 1 thread per step 2 job      
  --covarColList PC{1:10},age,age_squared,sex \                                             # Quantitative covariates to include. Only standard default covariates are listed here. We have to include again or REGENIE tries to include the phenotype as a covariate and fails.
  --catCovarList wes_batch \                                                                # Categorical covariates to include. Only standard default covariates are listed here.  
  --out /test/<output_prefix>.<tarball_prefix>.<chromosome>                                 # Name the outfile                                                                
```

REGENIE also (when requested with the `run_marker_tests` flag) performs per-marker tests. These are run identically to step 2
as shown above, without mask and annotation definitions.

#### GLMs

The method to perform Generalised Linear Models (GLMs) as part of this applet is implemented using the [statsmodels](https://www.statsmodels.org/stable/index.html)
python package.

##### Inputs

Variant data is identical to that for [STAAR](#staar). To ingest this data into Python, we convert the R object for 
staar into a sparse matrix using the script `sparseMatrixProcessor.R` in the `resources/usr/bin/` directory. This script 
creates a tab-delimited file of all non-reference genotypes with ENST information like:

```text
FID varID   gt  ENST
1000000 1:12345:A:T 1   ENST000000001
1000001 1:12345:A:T 2   ENST000000001
1000003 1:12367:G:C 1   ENST000000002
```

Where 'gt' is 1 if an individual is heterozygous and 2 if an individual is homozygous.

##### Command line example

This method has no command line as it is run within the applet's source code. For specifics on the method used, please see 
the applet source code and the [statsmodels glm README](https://www.statsmodels.org/stable/glm.html). Briefly, this runs in 
two basic steps:

```python
import pandas as pd
import statsmodels.api as sm

# 1. read file (1) from above as a pandas dataframe:
geno_table = pd.read_csv("lm.tsv",
                         sep = "\t",
                         names = ['FID', 'varID', 'gt', 'ENST'])

genes = ["ENST0000000001", "ENST0000000002"]
is_binary = True
pheno_name = 'T2D'
pheno_covars = pd.DataFrame() # This is actually filled with covariate and phenotype information per-participant

# 2. Loop through all genes from file (lm.tsv):
for gene in genes:
    if is_binary == True:
       family = sm.familes.Binomial()
    else:
       family = sm.families.Gaussian()
    
    sm.GLM.from_formula(pheno_name + ' ~ has_var + sex + age + batch + wes_batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10', 
                        data=pheno_covars,
                        family=family).fit()
```

The 'has_var' variable in the generalised linear model is an additive variable where individuals are coded as 0, 1, 2, 3, ...
depending on the number of variants they have in a given gene/gene set.

### Variant Extraction

This module will extract:

1. The exact phenotype/covariate information in tabular format that was used as input for burden tests
2. Individual-level information describing genotype-sample pairs
3. Variant-level information for all variants included in a test
4. A basic linear model to confirm the correct mask/maf combination was used to generate the results output by this applet

This module uses an IDENTICAL methodology to run GLMs as in the GLM `burden` mode as described above. Please see the
documentation for [that section](#glms) on the basics of running GLM-based burden tests on the RAP.

### PheWAS

This module will run a PheWAS on ALL phenotypes in the provided 'phenotype' file. This module uses an IDENTICAL methodology 
to run GLMs as in the GLM `burden` mode as described above. Please see the documentation for [that section](#glms) 
on the basics of running GLM-based burden tests on the RAP.

## Running on DNANexus

### Inputs

There are a standard set of command-line inputs that may be useful to change for the typical user. Note that the required
inputs change depending on the selected `mode`. Please see the detailed section on [mode](#mode) below for more information.

| input                   | description                                                                                                                                                                                                                                                                                                              |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| association_tarballs    | Hash ID(s) of the output from [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants) that you wish to use for rare variant burden testing. See below for more information.                                                                                                                  |
| tool                    | Tool to use for the burden testing module. **MUST** be one of 'bolt', 'saige', 'staar', 'glm', or 'regenie'. Case must match.                                                                                                                                                                                            |
| mode                    | Mode to run this applet in. **MUST** be one of 'burden', 'extract', or 'phewas'. Case must match.                                                                                                                                                                                                                        |
| gene_ids                | A comma-separated list of genes (either SYMBOL or ENST) to use for phewas or extract mode.                                                                                                                                                                                                                               |
| phenofile               | Phenotype file(s) – see below for more information on the format of these file(s). This can either be a single file or an array of file IDs when running PheWAS mode. See the [DNANexus Documentation](https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications) for what this means. |
| phenoname               | A single phenotype name to run association tests for. This allows for a user to provide a single phenotype file with multiple phenotypes and select one phenotype to run.                                                                                                                                                |
| run_marker_tests        | run SAIGE/BOLT/REGENIE per-marker tests? Note that tests with SAIGE currently take a VERY long time. **[TRUE]**                                                                                                                                                                                                          |
| is_binary               | Is the given trait in the phenofile binary?                                                                                                                                                                                                                                                                              |
| sex                     | Run only one sex or both sexes be run (0 = female, 1 = male, 2 = both) **[2]**?                                                                                                                                                                                                                                          |
| inclusion_list          | List of samples (eids) to include in analysis **[None]**                                                                                                                                                                                                                                                                 |
| exclusion_list          | List of samples (eids) to exclude in analysis **[None]**                                                                                                                                                                                                                                                                 |
| output_prefix           | Prefix to use for naming output tar file of association statistics. Default is to use the file name 'assoc_stats.tar.gz'                                                                                                                                                                                                 |
| covarfile               | File containing additional covariates to correct for when running association tests                                                                                                                                                                                                                                      |
| categorical_covariates  | comma-delimited list of categorical covariates found in covarfile to include in this model                                                                                                                                                                                                                               |
| quantitative_covariates | comma-delimited file of quantitative covariates found in covarfile to include in this model                                                                                                                                                                                                                              |

There are also several command-line inputs that should not need to be changed if running from within application 9905. These
mostly have to do with the underlying inputs to models that are generated by other tools in this pipeline. We have set
sensible defaults for these files and only change them if running from a different set of filtered data.

| input             | description                                                                                                                                | default file (all in `project-G6BJF50JJv8p4PjGB9yy7YQ2`) | 
|-------------------|--------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------|
| bgen_index        | index file with information on filtered and annotated UKBB variants                                                                        | `file-G86GJ3jJJv8fbXVB9PQ2pjz6`                          |
| transcript_index  | Tab-delimited file of information on transcripts expected by runassociationtesting output                                                  | `file-G7xyzF8JJv8kyV7q5z8VV3Vb`                          |
| base_covariates   | base covariates (age, sex, wes_batch, PC1..PC10) file for all WES UKBB participants                                                        | `file-G7PzVbQJJv8kz6QvP41pvKVg`                          |
| bed_file          | plink .bed format file from UKBB genetic data, filtered according to [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms) | `file-GBVYX88J57yK9fjf5qZ0kK61`                          |
| fam_file          | corresponding .fam file for 'bed_file'                                                                                                     | `file-GBVYYK0J57y8f2qZ2yfz09f5`                          |
| bim_file          | corresponding .bim file for 'bed_file'                                                                                                     | `file-GBVYYKQJ57yF1zKk15bZ3Yjb`                          |
| low_MAC_list      | list of low MAC (<100) variants in 'bed_file'                                                                                              | `file-GBVYYVjJ57y17ZyY2y65gFFk`                          |
| sparse_grm        | a sparse GRM for all individuals in 'bed_file' provided by [Bycroft et al.](https://www.nature.com/articles/s41586-018-0579-z)             | `file-GBVYYQjJ57yJ203GKx5KYxJ6`                          |
| sparse_grm_sample | corresponding samples in 'sparse_grm'                                                                                                      | `file-GBVYYV8J57yGZ20kKyFjqQ6P`                          |

#### Mode

As discussed above, this applet can be run in one of three modes with each of these modes requiring slightly different inputs:

1. `burden` - requires the `tool` option. Providing `gene_ids` does nothing.
2. `extract` – requires `gene_ids` UNLESS providing a SNP/GENE tarball from [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants.git)
3. `phewas` – requires `gene_ids` UNLESS providing a SNP/GENE tarball from [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants.git)

**Note:** – tool outputs and functionality is further changed depending on whether the provided tarball was generated using a 
list of variant IDs from mrcepid-collapsevariants.

**Note:** – PheWAS currently **CANNOT** run continuous and categorical phenotypes at the same time. You must provide
separate phenotype files and use the `is_binary` flag accordingly. 

#### Association Tarballs

`association_tarballs` takes either a single DNANexus file-ID that points to a single output from the 'mrcepid-collapsevariants' tool **OR** 
a file with multiple such DNANexus file-IDs. This file is one file-ID per line, like:

```text
file-G7z31B0J6F3ZbpGK6J0y2xxF
file-G7z2zP8JQy2PjKjY97Zvb3z9
file-G7z2pXjJ91Jx1ZJ4335bX52Z
file-G7z2gBjJ9ZpVqP3098X3xK6Y
file-G7z2VZQJ845v91751z2k9v2B
file-G7yv39QJPY6bYJKY9776JJBG
file-G7yqY10JkX6Y4fzvJf80z7P6
file-G7yq4ZjJV8qPjKjY97ZvZJ85
file-G7ypg98Jq133zgzY1yVy8gFz
```

Where each line is a different mask. An example file that includes 16 variant masks (8 variant types at two MAF cutoffs) 
is available at `collapsed_variants_new/variant_mask_list.txt (file-G7zPvZ0JJv8v06j8Gv2ppxpJ)` in project 
`project-G6BJF50JJv8p4PjGB9yy7YQ2`. Individual masks are available in `collapsed_variants_new/`. Please see the 
[high-level documentation](https://github.com/mrcepid-rap#collapsed-variants) for all apps for more information on
pre-collapsed variant files.

##### SNP and GENE Tarballs

Users can also provide a tarball from [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants.git)
that was created using a SNP or GENE list. Only *one* tar can be used at a time when running in SNP or GENE list mode.

#### Phenotypes File

A tab-delimited file with three columns and a header. Column 1 and 2 **MUST** be named FID and IID in that order. Column 3 can
have any non-whitespace name (`pheno_name` below) that represents the phenotype that you want to test. The applet automatically 
uses this name when performing variant testing and creating output. If providing a file with multiple phenotype columns, 
users **MUST** choose that phenotype to test in `burden` mode using the `phenoname` input. `phenoname` must exactly match 
the desired column. For `phewas` mode, `phenoname` is not required and **ALL* phenotypes in the provided file will be tested.
Values for binary traits MUST be either 0/1/NA/NaN while values for continuous traits can be any float/int (e.g. 1.42 / 5) 
or NA/NaN. Individuals with NA/NaN values are automatically excluded during testing.

When running in PheWAS mode, the user can also provide a .JSON compatible file-array of multiple phenotype files on which to
perform PheWAS. Please see the [DNANexus Documentation](https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications) 
for how to do this in the UKBB RAP.

```text
FID IID pheno_name
1000000 1000000 0
1000001 1000001 1
1000002 1000002 1
1000003 1000003 0
1000003 1000003 NA
```

#### Inclusion / Exclusion lists

These files are single-row .txt files with one eid per line. I have created three files already, but any file that 
conforms to the proper format can be used. The following files are already available on the RAP and have already been 
restricted to samples that have WES:

| name                               | file ID                       | description                                                   | list length  |
|------------------------------------|-------------------------------|---------------------------------------------------------------|--------------|
| EXCLUDEFOR_Relateds.txt            | file-GBVYYQ0J57y5gp4j5q7gP2Qg | List of related individuals (all ancestries)                  | 65,575       |
| EXCLUDEFOR_White_Euro_Relateds.txt | file-GBVYYP0J57y3Z5395g7q33z6 | List of related and NON-european ancestry individuals         | 96,271       |
| KEEPFOR_White_Euro.txt             | file-GBVYYPQJ57y5ZV904gzFb873 | List of all European ancestry individuals (including related) | 421,839      |

#### Additional Covariate (Quantitative / Categorical) File

This tool can also take as optional input a file (via `covarfile`) that contains either quantitative or categorical covariates to correct
for when running models. This file takes a similar input structure to the file provided to `phenofile`:

```text
FID IID quantCovar1  catCovar2  quantCovar3
1234567 1234567 0.1 centre1 4
7654321 7654321 1.3 centre2 9
8193498 8193498 0.4 centre1 2
1945892 1945892 NaN NA  NaN
```

Quantitative fields must be either floats or integers. Categorical covariates can be any form so long as they can be coerced
to a string and have relative few levels (I do not know what a reasonable limit on levels is, but assume that runtime will
rapidly increase as a factor of levels in your categorical covariates). Samples with **ANY** NaN/NA/blank fields will be 
skipped by the parser. The user must also provide a comma seperated list of covariates to include. Using the example 
above, all covariates will be included with the following command line:

```
-icovarfile=covars.txt -icategorical_covariates=catCovar2 quantitative_covariates=quantCovar1,quantCovar3
```

#### Gene IDs

Users can provide a comma delimited list of Gene Symbols **OR** Transcript ENST IDs for use during extract/phewas modes.
For example, the following inputs do the exact same thing:

`-igene_ids=BRCA1,BRCA2,PALB` OR `-gene_ids=ENST00000357654,ENST00000380152,ENST00000261584`

Note that if using SYMBOL, the ID **MUST** match exactly and not represent more than one transcript. All valid SYMBOLs
and ENST IDs must also be in the file specified by `transcript_index` (default: `file-G7xyzF8JJv8kyV7q5z8VV3Vb`) 

### Outputs

| output                   | description                            |
|--------------------------|----------------------------------------|
| output_tarball           | Output tarball containing test results |

output_tarball is either named `assoc_results.tar.gz` by default. If the parameter `output_prefix` is provided a file like (set `-ioutput_prefix="PTV"`):

`PTV.assoc_results.tar.gz`

Would be created. This tar.gz file will contain files that are specific to the tool and mode that was requested. More 
information on some of these outputs is given below.

1. burden:
   1. `<file_prefix>.genes.<TOOL>.stats.tsv.gz` (per-gene output)
   2. `<file_prefix>.genes.<TOOL>.stats.tsv.gz.tbi` (per-gene output index)
   3. `<file_prefix>.marker.<TOOL>.stats.tsv.gz` (per-marker output [when requested for BOLT / SAIGE / REGENIE])
   4. `<file_prefix>.marker.<TOOL>.stats.tsv.gz.tbi` (per-marker output index [when requested for BOLT / SAIGE / REGENIE])

    Note that some tools provide additional log/stat files that are not documented here, but are discussed in 
    tool-specific documentation.

2. extract:
   1. `<output_prefix>.<GENE>.variant_table.tsv` (a table of per-variant information similar to that generated by VEP annotation)
   2. `<output_prefix>.<GENE>.carriers_formated.tsv` (a table of individual-variant pairs with genotype)
   3. `<output_prefix>.genes.STAAR_glm.stats.tsv.gz` (a table similar to that for burden testing GLM & STAAR output providing summary results for all genes / phenotypes analysed)
   4. `<output_prefix>.genes.STAAR_glm.stats.tsv.gz.tbi` (the associated index)
   5. `<output_prefix>.phenotypes_covariates.formatted.tsv` (a table of individual/covariate/phenotype information identical to that used by internal methods for burden testing)

3. phewas:
   1. `<output_prefix>.genes.STAAR_glm.stats.tsv.gz` (per-gene/phenotype information identical to that for a burden GLM run)
   2. `<output_prefix>.genes.STAAR_glm.stats.tsv.gz.tbi` (the associated index)

**Note:** extract/phewas outputs will **NOT** be bgzipped or tabix indexed when running in SNP/GENE-list mode!

#### Per-gene output

A tab-delimited, gzipped file named like `<output_prefix>.genes.<tool>.stats.tsv.gz` (where `<output_prefix>` is identical
to that provided to the `output_prefix` input and `<tool>` is the name of the tool requesed by the `tool` input parameter) 
containing per-gene burden tests. An index for easy querying with tabix is also provided (`<output_prefix>.genes.<tool>.stats.tsv.gz.tbi`). 
Columns include those in the standard tool output. Additional columns contain per-gene information derived from in the
file `transcripts.tsv.gz (file-G7xyzF8JJv8kyV7q5z8VV3Vb)` in project `project-G6BJF50JJv8p4PjGB9yy7YQ2`. These columns include:
   
| column name       | description                                                                                                                                                                                                                                                                   |
|-------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ENST              | ENSMBL ENST ID. Will normally be the ENST ID that corresponds to MANE transcript or the ENSEMBL canonical transcript except in rare cases. When running in SNP or GENE list mode, this column will contain the dummy value of ENST0000000000 or ENST99999999999, respectively |
| chrom             | chromosome of this gene *without* the 'chr' prefix                                                                                                                                                                                                                            |
| start             | transcription start coordinate in hg38                                                                                                                                                                                                                                        |
| end               | transcription stop coordinate in hg38                                                                                                                                                                                                                                         |
| ENSG              | ENSEMBL ENSG corresponding to ENST                                                                                                                                                                                                                                            |
| MANE              | MANE v0.93 transcript                                                                                                                                                                                                                                                         |
| transcript length | end - start                                                                                                                                                                                                                                                                   |
| SYMBOL            | HGNC gene name                                                                                                                                                                                                                                                                |
| CANONICAL         | Is ENST the ENSEMBL canonical transcript?                                                                                                                                                                                                                                     |
| BIOTYPE           | Should *always* be protein_coding                                                                                                                                                                                                                                             |
| cds_length        | translation stop - translation start accounting for intron length                                                                                                                                                                                                             |
| coord             | formatted 'chrom:start-end' with chr prefix for aid in lookup on databases                                                                                                                                                                                                    |
| manh.pos          | relative position in the genome on a scale of 0-1 (chr1:1 = 0, chrY:14522573 = 1) for easy manhattan plot creation                                                                                                                                                            |

   **Note:** If a transcript ID is NOT found in `transcripts.tsv.gz`, this information **WILL NOT** be included in final
   output.

   Additional columns derived from the prefix of input tarfiles from `collapsevariants` will also be included. If the 
   tar is named like "HC_PTV-MAF_01", these additional columns will be headed as 'MASK' and 'MAF' respectively. This 
   functionality is trigged by the second value delimited by '-' having the prefix of 'MAF' or 'AC' (e.g. MAF_01, MAF_1,
   MAF_005, etc). Otherwise, will be labelled as `var1`, `var2`, etc. as DELIMITED by '-'. For example, if the tarfile 
   name is "Foo-Bar.tar.gz", column `var1` will include "Foo" and column `var2` will include "Bar". The software 
   currently does not have a method for naming these columns otherwise except for the special case as mentioned above.
  
#### Extract and PheWAS mode outputs

In addition to the fields listed above in [per-gene output](#per-gene-output), Extract and PheWAS modes have the following additional columns:

| column name            | binary only? | included in SNP/GENE masks | description                                                                                                              |
|------------------------|--------------|----------------------------|--------------------------------------------------------------------------------------------------------------------------|
| pheno                  | FALSE        | TRUE                       | Name of the phenotype tested in this column                                                                              |
| p_val_init             | FALSE        | TRUE                       | Initial p_value from GLM (**not used in extract/phewas modes**)                                                          |
| n_car                  | FALSE        | TRUE                       | Number of INDIVIDUALS with at least one variant in ENST                                                                  |
| cMAC                   | FALSE        | TRUE                       | Cumulative Minor Allele Count (i.e. taking into account individuals with multiple variants) for transcript given by ENST |
| n_model                | FALSE        | TRUE                       | Total number of individuals included in this model                                                                       |
| p_val_full             | FALSE        | TRUE                       | Full p_value from GLM                                                                                                    |
| effect                 | FALSE        | TRUE                       | Beta / log Odds Ratio from GLM                                                                                           |
| std_err                | FALSE        | TRUE                       | Standard error from GLM                                                                                                  |
| model_run              | FALSE        | TRUE                       | Did this column actually run through GLM (TRUE/FALSE)? Values w/FALSE should have NA for p_value/effect/std_err columns  |
| n_noncar_affected      | **TRUE**     | TRUE                       | Number of non-carriers of variants in ENST are affected by pheno                                                         |
| n_noncar_unaffected    | **TRUE**     | TRUE                       | Number of non-carriers of variants in ENST are NOT affected by pheno                                                     |
| n_car_affected         | **TRUE**     | TRUE                       | Number of carriers of variants in ENST are affected by pheno                                                             |
| n_car_unaffected       | **TRUE**     | TRUE                       | Number of carriers of variants in ENST are not affected by pheno                                                         |
| MASK                   | FALSE        | **FALSE**                  | Mask name from the original .tar.gz file. Not output when running SNP/GENE-list 'masks'                                  |
| MAF                    | FALSE        | **FALSE**                  | MAF cutoff for variants from the original .tar.gz file. Not output when running SNP/GENE-list 'masks'                    | 
| n_var                  | FALSE        | TRUE                       | Overall number of variants in ENST                                                                                       |
| relatedness.correction | FALSE        | TRUE                       | Was STAAR run with a GRM (TRUE/FALSE)? If FALSE number of individuals with pheno was likely too small to converge        |
| staar.O.p              | FALSE        | TRUE                       | STAAR omnbibus p. value                                                                                                  |
| staar.SKAT.p           | FALSE        | TRUE                       | STAAR SKAT p. value                                                                                                      |
| staar.burden.p         | FALSE        | TRUE                       | STAAR burden test p. value                                                                                               |
| staar.ACAT.p           | FALSE        | TRUE                       | STAAR ACAT p. value                                                                                                      |

'Binary only?' and 'included in SNP/GENE masks' indicate if the given column is in the provided table when running 
binary traits and SNP/GENE-lists, respectively. 

#### Per-marker output

A tab-delimited, gzipped file named like `<output_prefix>.markers.<tool>.stats.tsv.gz` (where `<output_prefix>` is identical
to that provided to the `output_prefix` input and `<tool>` is the name of the tool requesed by the `tool` input parameter) 
containing per-marker burden tests. An index for easy querying with tabix is also provided 
(`<output_prefix>.genes.BOLT.stats.tsv.gz.tbi`). Additional columns contain per-marker information contained in the
file `450k_vep.sorted.tsv.gz (file-G857Z4QJJv8x7GXfJ3y5v1qV)` in project `project-G6BJF50JJv8p4PjGB9yy7YQ2`. These columns 
are identical to those provided by [mrcepid-annotatecadd](https://github.com/mrcepid-rap/mrcepid-annotatecadd#outputs).
Note that this output is only produced for BOLT, SAIGE, or REGENIE, where requested.

### Command line example

There are two ways to acquire this applet:

1. As an **applet** – clone the repository from github and `dx build` an APPLET into your own workspace. If this is your first time doing 
this within a project other than "MRC - Variant Filtering", please see our organisational documentation on how to download
and build this app on the DNANexus Research Access Platform:

https://github.com/mrcepid-rap

2. As an **app** – use the app that has been provided in the DNANexus global namespace. This will ensure you are always using
the latest version and keeps you from having to manually update your local version. To be able to access this app, you
will need to be an authorised member of `org-mrc_epid_group_1_2`. Please contact Eugene Gardner if you would like to
be added!

**Note:** All commands below have been provided as if using option (2) above!

This app/applet can then be run in a few different ways:

1. Run a BOLT burden test:

```commandline
dx run app-mrcepid-runassociationtesting --priority low --destination results/ \ 
        -imode=burden
        -itool=bolt
        -iassociation_tarballs=file-G7zPvZ0JJv8v06j8Gv2ppxpJ \
        -iphenofile=file-G5PFv00JXk80Qfx04p9X56X5 \
        -iis_binary=true \
        -iinclusion_list=file-G6qXvjjJ2vfQGPp04ZGf6ygj
        -ioutput_prefix="T2D.bolt" \
        -isex=2
```

**Note:** The above command line DOES NOT have an inclusion/exclusion list. Be sure to consider the individuals that you
keep for your model!

**Note:** boolean parameters must be provided in json-compatible case ('true' or 'false').

2. Extract HC PTVs and samples for the gene BRCA2 from a single tarball:

```commandline
dx run app-mrcepid-runassociationtesting --priority low --destination results/ \ 
        -imode=extract
        -igene_ids=BRCA2
        -iassociation_tarballs=file-G7yYqG0JKx5xg1qP8B6zq0z4 \
        -iphenofile=file-G5PFv00JXk80Qfx04p9X56X5 \
        -iis_binary=true \
        -iinclusion_list=file-G6qXvjjJ2vfQGPp04ZGf6ygj
        -ioutput_prefix="HC_PTV.T2D" \
        -isex=2
```

3. Run a PheWAS on multiple cancer types for multiple masks for BRCA2:

```commandline
dx run app-mrcepid-runassociationtesting --priority low --destination results/ \ 
        -imode=phewas
        -igene_ids=BRCA2
        -iassociation_tarballs=file-G7zPvZ0JJv8v06j8Gv2ppxpJ \
        -iphenofile=file-G8v6YQ8JYVk8XjV31f6zB5z4 \
        -iis_binary=true \
        -iinclusion_list=file-G6qXvjjJ2vfQGPp04ZGf6ygj
        -ioutput_prefix="T2D.phewas" \
        -isex=2
```

Brief I/O information can also be retrieved on the command line:

```commandline
dx run app-mrcepid-runassociationtesting --help
```

#### Selecting an Instance Type

I have set a sensible (and tested) default for compute resources on DNANexus for running burden tests. This is baked into the json used for building
the app (at `dxapp.json`) so setting an instance type when running any of these tools is unnecessary. This current 
default is for a mem3_ssd1_v2_x64 instance (64 CPUs, 512 Gb RAM, 2400Gb storage).

When using the other modules in this applet (`extract` or `phewas`) it may be pertinent to ease instance requirements. Any 
of the mem3_ssd1_v2_xYY instance types ae suitable for this applet, where YY represents a modified CPU requirement. As 
an example, when running `extract` for a single gene, it is only necessary to have 8 cpus (e.g. instance mem3_ssd1_v2_x8).

#### Batch Running

This applet is not compatible with batch running.

