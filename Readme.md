# RunAssociationTesting (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
    * [Background](#background)
        * [BOLT-LMM](#bolt-lmm)
        * [SAIGE-GENE](#saige-gene)
        * [STAAR](#staar)
        * [Generalised Linear Models (GLMs)](#generalised-linear-models--glms-)
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
        + [Association Tarballs](#association-tarballs)
        + [Phenotypes File](#phenotypes-file)
        + [Inclusion / Exclusion lists](#inclusion---exclusion-lists)
        + [Additional Covariate (Quantitative / Categorical) File](#additional-covariate--quantitative---categorical--file)
    * [Outputs](#outputs-4)
    * [Command line example](#command-line-example-2)
        + [Runtime Examples, System Requirements, and Output Expectations](#runtime-examples--system-requirements--and-output-expectations)
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
* [qctool](https://www.well.ox.ac.uk/~gav/qctool_v2/index.html)
* [bgenix](https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md)
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

This applet is step 6 (mrc-runassociationtesting) of the rare variant testing pipeline developed by Eugene Gardner for the UKBiobank
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

This applet uses at least two **TAB-DELIMITED** files for providing phenotypes/covariates:

1. A phenotype to test for an association to rare variants

Please see [inputs](#inputs) for how this file is structured.   

2. Standard covariates to control for during association testing:

The standard covariates file is currently hard-coded into the applet. This file has the DNANexus hash-id of: `file-G7PzVbQJJv8kz6QvP41pvKVg`

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
   
Individual rare variant burden tests are then run according to [user input](#inputs). 

### BOLT

BOLT roughly proceeds in two steps. First, BOLT computes a genetic relatedness matrix (GRM) from genetic data curated during
the previous step of this applet. Next, it tests for the association a variant provided in a bgen file for an association
with the phenotype of interest. Normally, this bgen file is a collection of phased SNP data from genotying chips. Here,
we have created a "dummy" bgen file that instead encodes presence/absence of a qualifying variant per-individual. This dummy
variable is then used by BOLT to test for an association with a given phenotype. We test for both per-gene and per-marker 
(i.e. SNV/InDels) association with our phenotype of interest.

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
chr1 ENST000000001 0 1000000
chr1 ENST000000002 0 2000000
```

This file is then converted to .bed using plink2 and then .bgen 1.2 8-bit, ref-last format using plink2:

```commandline
plink --make-bed --file <file_prefix>.BOLT --out <file_prefix>.BOLT
plink2 --export bgen-1.2 'bits='8 --bfile <file_prefix>.BOLT --out <file_prefix>.BOLT"
```

The above steps are done per-chromosome and provided via the `--bgenSamplesFileList` argument as described below.

To perform per-marker tests, we simply convert the VCF file used for [SAIGE](#saige-gene) into bgen format using plink2, 
and run exactly the same command as described below. Inputs and outputs are essentially identical. 

#### Command Line Example

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

#### Outputs

Four files:

1. A tab-delimited, gzipped file named like `<output_prefix>.genes.BOLT.stats.tsv.gz` (where `<output_prefix>` is identical
   to that provided to the `output_prefix` input) containing per-gene burden tests. An index for easy querying with tabix
   is also provided (`<output_prefix>.genes.BOLT.stats.tsv.gz.tbi`). Columns include those in the standard 
   BOLT [output](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html#x1-470008) run with `--lmmInfOnly` 
   (i.e. output *excludes* the `P_BOLT_LMM` column). Additional columns contain per-gene information derived from in the
   file `transcripts.tsv.gz (file-G7xyzF8JJv8kyV7q5z8VV3Vb)` in project `project-G6BJF50JJv8p4PjGB9yy7YQ2`. These columns include:
   
| column name       | description |
| ----------------- | ----------- |
| ENST              | ENSMBL ENST ID. Will normally be the ENST ID that corresponds to MANE transcript or the ENSEMBL canonical transcript except in rare cases |
| chrom             | chromosome of this gene *without* the 'chr' prefix |
| start             | transcription start coordinate in hg38 |
| end               | transcription stop coordinate in hg38 |
| ENSG              | ENSEMBL ENSG corresponding to ENST |
| MANE              | MANE v0.93 transcript |
| transcript length | end - start |
| SYMBOL            | HGNC gene name |
| CANONICAL         | Is ENST the ENSEMBL canonical transcript? |
| BIOTYPE           | Should *always* be protein_coding |
| cds_length        | translation stop - translation start accounting for intron length |
| coord             | formatted 'chrom:start-end' with chr prefix for aid in lookup on databases |
| manh.pos          | relative position in the genome on a scale of 0-1 (chr1:1 = 0, chrY:14522573 = 1) for easy manhattan plot plotting |

   **BIG NOTE:** If a transcript ID is NOT found in `transcripts.tsv.gz`, this information WILL NOT be included in final
   output.

   Additional columns derived from the prefix of input tarfiles from `collapsevariants` will also be included. These columns 
   will be labelled as `var1`, `var2`, etc. DELIMITED by '-'. For example, if the tarfile name is "HC_PTV-MAF_01.tar.gz", 
   column `var1` will include "HC_PTV" and column `var2` will include "MAF_01". The software currently does not have a method
   for naming these columns otherwise.
  
2. A tab-delimited, gzipped file named like `<output_prefix>.markers.BOLT.stats.tsv.gz` (where `<output_prefix>` is identical
   to that provided to the `output_prefix` input) containing per-marker burden tests. An index for easy querying with tabix
   is also provided (`<output_prefix>.genes.BOLT.stats.tsv.gz.tbi`). Columns include those in the standard
   BOLT [output](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html#x1-470008) run with `--lmmInfOnly`
   (i.e. output *excludes* the `P_BOLT_LMM` column). Additional columns contain per-marker information contained in the
   file `450k_vep.sorted.tsv.gz (file-G857Z4QJJv8x7GXfJ3y5v1qV)` in project `project-G6BJF50JJv8p4PjGB9yy7YQ2`. These columns 
   are identical to those provided by [mrcepid-annotatecadd](https://github.com/mrcepid-rap/mrcepid-annotatecadd#outputs).
   
3. The standard BOLT log file named like `<output_prefix>.BOLT.log` (where `<output_prefix>` is identical
   to that provided to the `output_prefix` input).
   
4. The genotype-based results from the first step of BOLT named like `<output_prefix>.bolt.stats.gz` (where `<output_prefix>` 
   is identical to that provided to the `output_prefix` input).

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
ENST000000001 1:1000000_A/T 1:1000010_T/G   1:1000020_G/A
ENST000000002 1:2000000_G/C 1:2000030_A/C   1:2000050_ATC/A 1:2000000_G/GATC 
```

Files are created per-chromosome to enable faster parallelization on DNA Nexus.

#### Command Line Example

SAIGE-GENE proceedes in two steps:

1. Fitting the null GLMM (done for the whole genome):

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
          --traitType=binary/quantitative \                                  # set according to type of trait as provided by user
          --useSparseGRMtoFitNULL=TRUE                                       # Use the sparse GRM to fit the null model rather than estimating the variance ratio from genetic data (massively improves runtime)
```

**Note:** I have shortened the name of the sparseGRMFile for readability (see source code).

2. Performing rare variant burden tests (done per chromosome):

```commandline
step2_SPAtests.R
          --vcfFile=saige_input.vcf.gz \                                   # Input vcf file I document above
          --vcfField=GT \                                                  # Hardcoded INFO field to check for presence absence in the vcf file
          --GMMATmodelFile=SAIGE_OUT.rda \                                 # File generated by step1 above
          --varianceRatioFile=SAIGE_OUT_cate.varianceRatio.txt \           # File generated by step1 above
          --LOCO=FALSE \                                                   # Should we do leave-one-chrom-out (NO!)
          --SAIGEOutputFile=<file_prefix>.<chr>.SAIGE_OUT.SAIGE.gene.txt \ # Output file from this step
          --groupFile=<file_prefix>.<chr>.SAIGE.groupFile.txt \            # Input groupFile I document above
          --sparseSigmaFile=SAIGE_OUT_cate.sparseSigma.mtx \               # File generated by step1 above
          --IsSingleVarinGroupTest=TRUE \                                  # Should we also test individual variants (YES!)
          --MACCutoff_to_CollapseUltraRare=0.5 \                           # Minimum allele count to include a variant. 0.5 means we include all variants (even singletons)
          --IsOutputHetHomCountsinCaseCtrl " \                             # Output Het/Hom counts to help calculate MAC
          --IsOutputNinCaseCtrl " \                                        # Output number of case/controls
          --maxMAFforGroupTest=1"                                          # Include all variants (we define our own cutoffs) 
```

**Note:** I have shortened the name of the sparseSigmaFile for readability (see source code).

#### Outputs

!!!TO DO!!!

### STAAR

STAAR is unique here in that it does not have a command-line tool associated with it and is instead provided as an R package.
To run STAAR, we created two R wrapper scripts that 

1. Generates a null model (`runSTAAR_Null.R`).
2. Ingests the data required to run burden tests, and formats the phenotypic and genetic data into tables / matrices that are compatible with
the STAAR package (`runSTAAR_Genes.R`).
   
These commented scripts are provided in this repository at:

`resources/usr/bin/`

#### Inputs

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
[provided to SAIGE-GENE](#saige-gene).

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

#### Command line example

This example command-line is to run both scripts that we have created for this applet for one chromosome. All commands are
space delimited:

```commandline
Rscript runSTAAR_Null.R /
          phenotypes_covariates.formatted.txt \               # Formated phenotypes/covariates to derive the NULL STAAR model
          <pheno_name> \                                      # phenotype name extracted from the provided phenotype file
          <is_binary> \                                       # Is the trait binary (false for quantitative)?
          <quant_covars> \                                    # Names of additional quantitative covariates to include in model (NULL for none)
          <catagorical_covars>                               # Names of additional catagorical covariates to include in model (NULL for none)

Rscript runSTAAR_Genes.R /                                     
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

#### Outputs

We output a single tab-delimited file for all gene/mask/MAF combinations (`<output_prefix>.STAAR_results.tsv.gz`). A tabix
index for easy querying is also provided (`*.tbi`). Per-transcript information is identical to that described [for BOLT](#outputs), 
above. STAAR-specific columns include:

1. n.samps – number of samples run through STAAR
2. pheno – name of the phenotype from the file provided to `pheno_file` 
3. staar.O.p – Overall p. value
4. staar.SKAT.p – SKAT p. value
5. staar.burden.p – Burden p. value
6. staar.ACAT.p ACAT-V p. value
7. n.var – Number of variants considered for this gene
8. cMAC – Cumulative minor allele count of all variants considered for this gene

We also output the null model in R RDS format (`<output_prefix>.STAAR_null.rds`).

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
1000000 ENST000000001 0/1
1000001 ENST000000001 0/1
1000003 ENST000000002 0/1
```

2. A list of genes to process:

```text
ENST000000001
ENST000000002
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

genes = ["ENST0000000001", "ENST0000000002"]
is_binary = True
pheno_name = 'T2D'
pheno_covars = pd.DataFrame() # This is actually filled with covariate and phenotype information per-participant

# 2. Loop through all genes from file (2):
for gene in genes:
    if is_binary == True:
       family = sm.familes.Binomial()
    else:
       family = sm.families.Gaussian()
    
    sm.GLM.from_formula(pheno_name + ' ~ has_var + sex + age + batch + wes_batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10', 
                        data=pheno_covars,
                        family=family).fit()
```

#### Outputs

!!!TO DO!!!

## Running on DNANexus

### Inputs

There are a standard set of command-line inputs that may be useful to change for the typical user. 

|input|description             |
|---- |------------------------|
|association_tarballs  | Hash ID(s) of the output from [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants) that you wish to use for rare variant burden testing. See below for more information. |
|phenofile | Phenotypes file – see below for more information on the format of this file |
|run_bolt    | run BOLT? **[FALSE]** |
|run_saige   | run SAIGE? **[FALSE]** |
|run_staar   | run STAAR? **[FALSE]** |
|run_linear_model   | run GLMs? **[FALSE]** |
|run_all | run all of the above **[FALSE]** |
|run_marker_tests| run SAIGE per-marker tests? **[FALSE]**|
|is_binary   | Is the given trait in the phenofile binary? |
|sex   | Run only one sex or both sexes be run (0 = female, 1 = male, 2 = both) **[2]**? |
|inclusion_list | List of samples (eids) to include in analysis **[None]** |
|exclusion_list | List of samples (eids) to exclude in analysis **[None]** |
|output_prefix   | Prefix to use for naming output tar file of association statistics. Default is to use the file name 'assoc_stats.tar.gz' |
|covarfile | File containing additional covariates to correct for when running association tests |
|categorical_covariates | comma-delimited list of categorical covariates found in covarfile to include in this model |
|quantitative_covariates | comma-delimited file of quantitative covariates found in covarfile to include in this model |

There are also several command-line inputs that should not need to be changed if running from within application 9905. These
mostly have to do with the underlying inputs to models that are generated by other tools in this pipeline. We have set
sensible defaults for these files and only change them if running from a different set of filtered data.

|input              |description             | default file (all in `project-G6BJF50JJv8p4PjGB9yy7YQ2`) | 
|-------------------|------------------------| ------- |
| bgen_index        | index file with information on filtered and annotated UKBB variants      | `file-G86GJ3jJJv8fbXVB9PQ2pjz6` |
| transcript_index  | Tab-delimited file of information on transcripts expected by runassociationtesting output | `file-G7xyzF8JJv8kyV7q5z8VV3Vb` |
| fam_file          | base covariates (age, sex, wes_batch, PC1..PC10) file for all WES UKBB participants       | `file-G7PzVbQJJv8kz6QvP41pvKVg` |
| bed_file          | plink .bed format file from UKBB genetic data, filtered according to [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms)  | `file-G6qXq38J2vfXKBz44Z8bxf5V` |
| fam_file          | corresponding .fam file for 'bed_file'                                   | `file-G6qXvg0J2vffF7Y44VFJb7jB` |
| bim_file          | corresponding .bim file for 'bed_file'                                   | `file-G6qXvgQJ2vfqY3z64x7x8jPq` |
| low_MAC_list      | list of low MAC (<100) variants in 'bed_file'                            | `file-G6qXvq8J2vfYY6y64Qy23v5b` |
| sparse_grm        | a sparse GRM for all individuals in 'bed_file' created by SAIGE 'Step0'  | `file-G6xK1xjJb3g9vz0p236vZ9F5` |
| sparse_grm_sample | corresponding samples in 'sparse_grm'                                    | `file-G6xK1z8Jb3g1Fp4j2311JKYF` |

#### Association Tarballs

`association_tarballs` takes a list file of DNANexus file-IDs that point to multiple outputs from the 'mergecollapsevariants'. 
This file is one file-ID per line, like:

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

#### Phenotypes File

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

#### Inclusion / Exclusion lists

These files are single-row .txt files with one eid per line. I have created three files already, but any file that 
conforms to the proper format can be used. The following files are already available on the RAP and have already been 
restricted to samples that have WES:

| name | file ID | description | list length |
| ---- | ------- | ----------- | ----------- |
| EXCLUDEFOR_Relateds.txt | file-G6qXvkjJ2vfY7yF74VVPG7xg | List of related individuals (all ancestries) | 65,575 |
| EXCLUDEFOR_White_Euro_Relateds.txt | file-G6qXvj8J2vfz37qZ4ZPf7p2Z | List of related and NON-european ancestry individuals | 96,271 |
| KEEPFOR_White_Euro.txt | file-G6qXvjjJ2vfQGPp04ZGf6ygj | List of all European ancestry individuals (including related) | 421,839 |

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

**NOTE:** Categorical covariates are never corrected for when running SAIGE!

### Outputs

|output                 | description       |
|-----------------------|-------------------|
|output_tarball         |  Output tarball containing test results  |

output_tarball is either named `assoc_stats.tar.gz` by default. If the parameter `output_prefix` is provided a file like (set `-ioutput_prefix="PTV"`):

`PTV.assoc_stats.tar.gz`

Would be created. This tar.gz file will contain files that are specific to the tool that was requested. 

1. `<file_prefix>.SAIGE_OUT.SAIGE.gene.txt` (SAIGE-GENE output)
2. `<file_prefix>.SAIGE_OUT.SAIGE.gene.txt_single` (SAIGE-GENE per-marker output)
3. `<file_prefix>.SAIGE_OUT.rda` (SAIGE-GENE null model file)
4. `<file_prefix>.SAIGE_OUT_cate.varianceRatio.txt` (SAIGE-GENE null model file)
5. `<file_prefix>.STAAR_results.tsv` (STAAR output)
6. `<file_prefix>.STAAR_null.rds` (STAAR null model in .rds format)
7. `<file_prefix>.genes.BOLT.stats.tsv.gz` (BOLT per-gene output)
8. `<file_prefix>.genes.BOLT.stats.tsv.gz.tbi` (BOLT per-gene output index)
9. `<file_prefix>.marker.BOLT.stats.tsv.gz` (BOLT per-marker output)
10. `<file_prefix>.marker.BOLT.stats.tsv.gz.tbi` (BOLT per-marker output index)
11. `<file_prefix>.stats.gz` (BOLT raw genotype output)
12. `<file_prefix>.lm_stats.tsv` (GLM output)

Please see each tool's respective output section for more information on these outputs.

### Command line example

If this is your first time running this applet within a project other than "MRC - Variant Filtering", please see our
organisational documentation on how to download and build this app on the DNANexus Research Access Platform:

https://github.com/mrcepid-rap

This applet can be run in a few different ways:

1. Run one tool/method at a time:

```commandline
dx run mrcepid-runassociationtesting --priority low --destination results/ \ 
        -iassociation_tarballs=file-G7zPvZ0JJv8v06j8Gv2ppxpJ \
        -iphenofile=file-G5PFv00JXk80Qfx04p9X56X5 \
        -iis_binary=true \
        -ioutput_prefix="T2D.bolt" \
        -isex=2 \
        -irun_bolt=true
```

**Note**: The above command line DOES NOT have an inclusion/exclusion list. Be sure to consider the individuals that you
keep for your model!

2. Run all tools/methods:

```commandline
dx run mrcepid-runassociationtesting --priority low --destination results/ \ 
        -iassociation_tarballs=file-G7zPvZ0JJv8v06j8Gv2ppxpJ \
        -iphenofile=file-G5PFv00JXk80Qfx04p9X56X5 \
        -iis_binary=true \
        -ioutput_prefix="T2D.all" \
        -isex=2 \
        -irun_all=true
```

**Note:** boolean parameters must be provided in json-compatible case ('true' or 'false').

Brief I/O information can also be retrieved on the command line:

```commandline
dx run mrcepid-runassociationtesting --help
```

I have set a sensible (and tested) default for compute resources on DNANexus for running BOLT and SAIGE. This is baked into the json used for building
the app (at `dxapp.json`) so setting an instance type when running either of those tools or all tools at once is unnecessary.
This current default is for a mem3_ssd1_v2_x64 instance (64 CPUs, 512 Gb RAM, 2400Gb storage).

#### Runtime Examples, System Requirements, and Output Expectations

| Tool | Per-gene p. value | Per-variant p. value | Accurate β / OR | Runtime§ | Cost§ |
| ---- | ----------------- | -------------------- | --------------- | ------- | ---- |
| BOLT | <span style="color:green">**TRUE**</span> | <span style="color:green">**TRUE**</span> | <span style="color:green">**TRUE**</span> / <span style="color:red">**FALSE**</span> | 17hr21m | £7.33 |
| SAIGE | <span style="color:green">**TRUE**</span> | <span style="color:orange">**POSSIBLE**</span> | <span style="color:green">**TRUE**</span> / <span style="color:red">**FALSE**</span> | unk | unk |
| STAAR | <span style="color:green">**TRUE**</span> | <span style="color:red">**FALSE**</span> | <span style="color:red">**FALSE**</span> / <span style="color:red">**FALSE**</span> | 10hr53m | £4.60 |
| GLM | <span style="color:green">**TRUE**</span> | <span style="color:red">**FALSE**</span> | <span style="color:green">**TRUE**</span> / <span style="color:green">**TRUE**</span> | unk | unk |

§Numbers are representative runtimes and cost for a real analysis testing the association of rare variant burden on having type II diabetes.

#### Batch Running

This applet is not compatible with batch running.

