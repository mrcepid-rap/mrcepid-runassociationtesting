# RunAssociationTesting (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
  * [Changelog](#changelog)
  * [Background](#background)
  * [Dependencies](#dependencies)
    + [Docker](#docker)
    + [Resource Files](#resource-files)
- [Methodology](#methodology)
  * [Covariate processing](#covariate-processing)
- [Running on DNANexus](#running-on-dnanexus)
  * [Inputs](#inputs)
    + [DNANexus Inputs](#dnanexus-inputs)
    + [Command-Line Inputs](#command-line-inputs)
    + [Default Modules](#default-modules)
      + [Loading New Modules](#loading-new-modules)
    + [Phenotypes File](#phenotypes-file)
    + [Inclusion / Exclusion lists](#inclusion--exclusion-lists)
    + [Additional Covariate (Quantitative / Categorical) File](#additional-covariate-quantitative--categorical-file)
    + [Association Tarballs](#association-tarballs)
  * [Outputs](#outputs)
  * [Command line example](#command-line-example)
    + [Selecting an Instance Type](#selecting-an-instance-type)

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

* v1.3.0
  * Update makes major changes for how runassociationtesting is run:
    * These changes are largely invisible to the user and mostly affect the underlying code, _except_ for how inputs are provided to the set of tools. See [inputs](#inputs) for more information 
    * The 'runassociationtesting' app/applet is now an entrypoint for individual modules and provides a set of interfaces to implement different functionalities (e.g. burden/extract/etc.)
    * This means individual modules are now kept and maintained as separate repositories. A list of default modules is included as part of this repository in the [inputs section](#default-modules)
    * Users can now develop modules that use the interfaces provided as part of this repository to implement different association tests. Please see the [developer README](https://github.com/mrcepid-rap/mrcepid-runassociationtesting/blob/main/Readme.developer.md) included in this repo for more information

* v1.2.4
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
perform testing. This module enables several functions via a `-imode` flag. Please see the [command-line examples](#command-line-example)
below for more information on how this works. 

### Dependencies

#### Docker

This applet uses [Docker](https://www.docker.com/) to supply dependencies to the underlying AWS instance
launched by DNANexus. The Dockerfile used to build dependencies is available as part of the MRCEpid organisation at:

https://github.com/mrcepid-rap/dockerimages/blob/main/burdentesting.Dockerfile

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

This applet provides an interface to implement various genetic association tests and burden analyses. As part of this 
applet, it _always_ provides and processes a generic set of necessary data.

This applet proceeds in two basic steps:

1. Covariate/phenotype processing and individual-level filtering
2. Starting the requested module provided by `-imode=`

The user selects the variant test/mode that they want to run at runtime. Please see the [inputs](#inputs)
section for more information.

### Covariate processing

This applet uses at least two **TAB-DELIMITED** files for providing phenotypes/covariates:

1. A phenotype file to test for associations. This can be a single file with multiple phenotypes or multiple files with 
   a single phenotype. A specific phenotype to test from a file with multiple phenotypes can be specified via the command-line. 

Please see [inputs](#inputs) for how this file is structured.

2. Standard covariates to control for during association testing:

A default set of standard covariates is used by this applet.

In brief, this file contain a raw version of the database output from the UK Biobank data 
[implemented on the UKBiobank RAP](https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-phenotypic-data-as-a-file).
Currently, these are (with [showcase](https://biobank.ndph.ox.ac.uk/showcase/) ID if applicable):

* Age at assessment - [21003](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21003)
* Genetic sex - [22001](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22001)
* The First 10 genetic principal components - [22009.1-10](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22009)
* WES Batch - Generated by pulling .fam files for the respective UKBB WES batches and extracting individuals unique to each fam file
* Array QC Passed – A binary value (0/1) indicating that the given participant passed array and imputation QC as defined in the [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms) applet

These covariates are automatically included when performing most burden testing. Sex is only included as a covariate
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
   
Individual tools can then be selected with their own specific inputs and data processing according to [user input](#inputs).
Please see those individual module's README for more information.

## Running on DNANexus

### Inputs

This tool uses two 'styles' of inputs documented separately in each section below:

1. [DNANexus Inputs](#dnanexus-inputs): inputs using standard DNANexus conventions – e.g. `-iinput_name=input`
2. [Default Inputs](#command-line-inputs): inputs using command-line argparse conventions provided to the 'input_args' DNANexus input – e.g. `--input_name input`

#### DNANexus Inputs

This applet provides instructions for running association tests in two ways. It provides exactly two arguments that i) control
the type of analysis that is done via the requested module (`mode`) and ii) the name of the output (`output_prefix`). It
then uses a third input to provide information to the default set of commands that all modules provide and the requested 
module itself (`input_args`).

| input         | description                                                                                                              |
|---------------|--------------------------------------------------------------------------------------------------------------------------|
| mode          | Mode to run this applet in. **MUST** be one of 'burden', 'extract', or 'phewas'. Case must match.                        |
| output_prefix | Prefix to use for naming output tar file of association statistics. Default is to use the file name 'assoc_stats.tar.gz' |
| input_args    | Additional inputs to control association tests. See below for what this means.                                           |

#### Command-Line Inputs

These are a standard set of command-line inputs required for all modules provided using the `input_args` input described 
above. For how these inputs are provided, see the [examples below](#running-on-dnanexus). If the option is not required,
defaults for each option are provided in **[bold brackets]**. Boolean options are flags, and change to 'true' when
provided. 

| input                   | Boolean? | Required? | description                                                                                                                                                                          |
|-------------------------|----------|-----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| phenofile               | False    | **True**  | Phenotype file(s) – see below for more information on the format of these file(s). This can either be a single file or a space-delimited list of file IDs.                           |
| phenoname               | False    | False     | A single phenotype name to run association tests for. This allows for a user to provide a single phenotype file with multiple phenotypes and select one phenotype to run **[None]**. |
| is_binary               | **True** | False     | Is the given trait in the phenofile binary **[false]**?                                                                                                                              |
| sex                     | False    | False     | Run only one sex or both sexes be run (0 = female, 1 = male, 2 = both) **[2]**?                                                                                                      |
| inclusion_list          | False    | False     | File containing list of samples (eids) to include in analysis **[None]**                                                                                                             |
| exclusion_list          | False    | False     | File containing list of samples (eids) to exclude in analysis **[None]**                                                                                                             |
| covarfile               | False    | False     | File containing additional covariates to correct for when running association tests                                                                                                  |
| categorical_covariates  | False    | False     | space-delimited list of categorical covariates found in covarfile to include in this model                                                                                           |
| quantitative_covariates | False    | False     | space-delimited file of quantitative covariates found in covarfile to include in this model                                                                                          |
| transcript_index        | False    | **True**  | Tab-delimited file of information on transcripts expected by runassociationtesting output                                                                                            |
| base_covariates         | False    | **True**  | base covariates (age, sex, wes_batch, PC1..PC10) file for all WES UKBB participants                                                                                                  |

#### Default Modules

As discussed above, this applet can be run in one of three modes with each of these modes requiring slightly different 
inputs. A default set of modules is listed below. For specific inputs for these modules, please see their individual
README documents located within their respective GitHub repositories.

1. `burden` - https://github.com/mrcepid-rap/mrcepid-runassociationtesting-burden
2. `extract` – https://github.com/mrcepid-rap/mrcepid-runassociationtesting-extract
3. `phewas` – https://github.com/mrcepid-rap/mrcepid-runassociationtesting-phewas

##### Loading New Modules

This applet will only take as input modules that have been provided when this applet was build on the DNANexus platform.
If you would like to add additional modules, you will need to tell the applet to include it via the dxapp.json. Please 
see the [developer README](https://github.com/mrcepid-rap/mrcepid-runassociationtesting/blob/main/Readme.developer.md)
for how this works.

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
above, all covariates will be included with the following command line provided to `-iinput_args=`:

```
--covarfile covars.txt --categorical_covariates catCovar2 quantitative_covariates quantCovar1 quantCovar3
```

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

### Outputs

| output                   | description                            |
|--------------------------|----------------------------------------|
| output_tarball           | Output tarball containing test results |

output_tarball is either named `assoc_results.tar.gz` by default. If the parameter `output_prefix` is provided a file like (set `-ioutput_prefix="PTV"`):

`PTV.assoc_results.tar.gz`

Would be created. This tar.gz file will contain files that are specific to the tool and mode that was requested. More 
information on some of these outputs is given in the respective modules' README.

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

This app/applet can then be run in a few different ways (`file-1234567890ABCDEFGHIJKLMN` is a placeholder, make sure to 
provide the actual required input!):

1. Run a BOLT burden test:

```commandline
dx run app-mrcepid-runassociationtesting --priority low --destination results/ \ 
        -imode=burden
        -ioutput_prefix="T2D.bolt" \
        -iinput_args=' \
                --tool bolt \
                --run_marker_tests \
                --association_tarballs file-1234567890ABCDEFGHIJKLMN \
                --phenofile file-1234567890ABCDEFGHIJKLMN \
                --is_binary \
                --inclusion_list file-1234567890ABCDEFGHIJKLMN \
                --sex 2 \
                --transcript_index file-1234567890ABCDEFGHIJKLMN \
                --base_covariates file-1234567890ABCDEFGHIJKLMN \
                --bgen_index file-1234567890ABCDEFGHIJKLMN \
                --bed_file file-1234567890ABCDEFGHIJKLMN \
                --fam_file file-1234567890ABCDEFGHIJKLMN \
                --bim_file file-1234567890ABCDEFGHIJKLMN \
                --low_MAC_list file-1234567890ABCDEFGHIJKLMN \
                --sparse_grm file-1234567890ABCDEFGHIJKLMN \
                --sparse_grm_sample file-1234567890ABCDEFGHIJKLMN`
```

**Note:** The above command line DOES NOT have an inclusion/exclusion list. Be sure to consider the individuals that you
keep for your model!

**Note:** boolean parameters are simple flags and do not require actual inputs (e.g. 'true' or 'false').

2. Extract HC PTVs and samples for the gene BRCA2 from a single tarball:

```commandline
dx run app-mrcepid-runassociationtesting --priority low --destination results/ \ 
        -imode=extract
        -ioutput_prefix="HC_PTV.T2D" \
        -iinput_args='
                --gene_ids BRCA2 CHEK2 \
                --association_tarballs file-1234567890ABCDEFGHIJKLMN \
                --phenofile file-1234567890ABCDEFGHIJKLMN \
                --is_binary \
                --inclusion_list file-1234567890ABCDEFGHIJKLMN \
                --sex 2 \
                --transcript_index file-1234567890ABCDEFGHIJKLMN \
                --base_covariates file-1234567890ABCDEFGHIJKLMN \
                --bgen_index file-1234567890ABCDEFGHIJKLMN \'
```

3. Run a PheWAS on multiple cancer types for multiple masks for BRCA2 (note that `--phenofile` includes MULTIPLE phenotypes in this example!):

```commandline
dx run app-mrcepid-runassociationtesting --priority low --destination results/ \ 
        -imode=phewas
        -ioutput_prefix="cancer.PHEWAS" \
        -iinput_args='
                --gene_ids BRCA2 \
                --association_tarballs file-1234567890ABCDEFGHIJKLMN \
                --phenofile file-1234567890ABCDEFGHIJKLMN \
                --is_binary \
                --inclusion_list file-1234567890ABCDEFGHIJKLMN \
                --sex 2 \
                --transcript_index file-1234567890ABCDEFGHIJKLMN \
                --base_covariates file-1234567890ABCDEFGHIJKLMN'
```

Brief I/O information can also be retrieved on the command line:

1. For DNANexus-style options:

```commandline
dx run app-mrcepid-runassociationtesting --help
```

2. For commandline-style options:

```commandline
dx run app-mrcepid-runassociationtesting -imode='burden' -ioutput_prefix='test' -iinput_args='--help'
```


#### Selecting an Instance Type

I have set a sensible (and tested) default for compute resources on DNANexus for running burden tests. This is baked into the json used for building
the app (at `dxapp.json`) so setting an instance type when running any of these tools is unnecessary. This current 
default is for a mem3_ssd1_v2_x64 instance (64 CPUs, 512 Gb RAM, 2400Gb storage).

When using the other modules in this applet it may be pertinent to ease instance requirements. Any 
of the mem3_ssd1_v2_xYY instance types ae suitable for this applet, where YY represents a modified CPU requirement. As 
an example, when running `extract` for a single gene, it is only necessary to have 8 cpus (e.g. instance mem3_ssd1_v2_x8).