import os
import csv
import dxpy
import gzip
import subprocess
import pandas as pd
import pandas.core.series

from association_pack import AssociationPack


# This function runs a command on an instance, either with or without calling the docker instance we downloaded
# By default, commands are not run via Docker, but can be changed by setting is_docker = True
# Also, by default, standard out is not saved, but can be modified with the 'stdout_file' parameter.
# print_cmd is for internal debugging purposes when testing new code
def run_cmd(cmd: str, is_docker: bool = False, stdout_file: str = None, print_cmd = False) -> None:

    # -v here mounts a local directory on an instance (in this case the home dir) to a directory internal to the
    # Docker instance named /test/. This allows us to run commands on files stored on the AWS instance within Docker.
    # This looks slightly different from other versions of this command I have written as I needed to write a custom
    # R script to run STAAR. That means we have multiple mounts here to enable this code to find the script.
    if is_docker:
        cmd = "docker run " \
              "-v /home/dnanexus:/test " \
              "-v /usr/bin/:/prog " \
              "egardner413/mrcepid-associationtesting " + cmd

    if print_cmd:
        print(cmd)

    # Standard python calling external commands protocol
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stdout_file is not None:
        with open(stdout_file, 'w') as stdout_writer:
            stdout_writer.write(stdout.decode('utf-8'))
        stdout_writer.close()

    # If the command doesn't work, print the error stream and close the AWS instance out with 'dxpy.AppError'
    if proc.returncode != 0:
        print("The following cmd failed:")
        print(cmd)
        print("STDOUT follows\n")
        print(stdout.decode('utf-8'))
        print("STDERR follows\n")
        print(stderr.decode('utf-8'))
        raise dxpy.AppError("Failed to run properly...")


# This is to generate a global CHROMOSOMES variable for parallelisation
def get_chromosomes(is_snp_tar: bool = False, is_gene_tar: bool = False):

    if is_snp_tar:
        chromosomes = list(['SNP'])
    elif is_gene_tar:
        chromosomes = list(['GENE'])
    else:
        chromosomes = list(range(1,23)) # Is 0 based on the right coordinate...? (So does 1..22)
        chromosomes.extend(['X'])
        chromosomes = list(map(str, chromosomes))

    return chromosomes


# This downloads and process a bgen file when requested
def process_bgen_file(chrom_bgen_index: dict, chromosome: str) -> None:

    # First we have to download the actual data
    bgen_index = dxpy.DXFile(chrom_bgen_index['index'])
    bgen_sample = dxpy.DXFile(chrom_bgen_index['sample'])
    bgen = dxpy.DXFile(chrom_bgen_index['bgen'])
    vep = dxpy.DXFile(chrom_bgen_index['vep'])
    dxpy.download_dxfile(bgen_index.get_id(), "filtered_bgen/" + chromosome + ".filtered.bgen.bgi")
    dxpy.download_dxfile(bgen_sample.get_id(), "filtered_bgen/" + chromosome + ".filtered.sample")
    dxpy.download_dxfile(bgen.get_id(), "filtered_bgen/" + chromosome + ".filtered.bgen")
    dxpy.download_dxfile(vep.get_id(), "filtered_bgen/" + chromosome + ".filtered.vep.tsv.gz")

    # And then do the filtering...
    cmd = "plink2 --threads 4 --bgen /test/filtered_bgen/" + chromosome + ".filtered.bgen 'ref-last' " \
            "--sample /test/filtered_bgen/" + chromosome + ".filtered.sample " \
            "--export bgen-1.2 'bits='8 " \
            "--out /test/" + chromosome + ".markers " \
            "--keep /test/SAMPLES_Include.txt"
    run_cmd(cmd, True)

    # The sample file output by plink2 is a disaster, so fix it here:
    os.rename(chromosome + '.markers.sample', chromosome + '.old')
    with open(chromosome + '.old', 'r') as samp_file:
        fixed_samp_bolt = open(chromosome + '.markers.bolt.sample', 'w')
        for line in samp_file:
            line = line.rstrip().split(" ")
            if line[0] == 'ID_1':
                fixed_samp_bolt.write('ID_1 ID_2 missing sex\n')
            elif line[3] == 'D':
                fixed_samp_bolt.write('0 0 0 D\n')
            else:
                fixed_samp_bolt.write("%s %s 0 NA\n" % (line[1], line[1]))
        samp_file.close()
        fixed_samp_bolt.close()

    # And finally index the file
    cmd = "bgenix -index -g /test/" + chromosome + '.markers.bgen'
    run_cmd(cmd, True)


def build_transcript_table() -> pandas.DataFrame:

    transcripts_table = pd.read_csv(gzip.open('transcripts.tsv.gz', 'rt'), sep = "\t")
    transcripts_table = transcripts_table.rename(columns={'#chrom':'chrom'})
    transcripts_table = transcripts_table.set_index('ENST')
    return(transcripts_table)


def get_gene_id(gene_id: str, transcripts_table: pandas.DataFrame) -> pandas.core.series.Series:

    if 'ENST' in gene_id:
        print("gene_id ??? " + gene_id + " ??? looks like an ENST value... validating...")
        try:
            gene_info = transcripts_table.loc[gene_id]
            print("Found one matching ENST (%s - %s)... proceeding..." % (gene_id, gene_info['coord']))
        except KeyError:
            print("Did not find a transcript with ENST value %s... terminating..." % gene_id)
    else:
        print("gene_id ??? " + gene_id + " ??? does not look like an ENST value, searching for symbol instead...")
        found_rows = transcripts_table[transcripts_table['SYMBOL'] == gene_id]
        if len(found_rows) == 1:
            found_enst = found_rows.index[0]
            gene_info = transcripts_table.loc[found_enst]
            print("Found one matching ENST (%s - %s) for SYMBOL %s... proceeding..." % (found_enst, gene_info['coord'], gene_id))
        elif len(found_rows) > 1:
            gene_info = None
            print("Found %i ENST IDs (%s) for SYMBOL %s... Please re-run using exact ENST to ensure consistent results..." % (len(found_rows), ','.join(found_rows.index.to_list()), gene_id))
            raise Exception("Multiple ENST IDs")
        else:
            gene_info = None
            print("Did not find an associated ENST ID for SYMBOL %s... Please re-run after checking SYMBOL/ENST used..." % (gene_id))
            raise Exception("ENST ID not found")

    return gene_info


def process_snp_or_gene_tar(is_snp_tar, is_gene_tar, tarball_prefix) -> tuple:

    file_prefix = None
    if is_snp_tar:
        print("Running in SNP mode...")
        file_prefix = 'SNP'
        gene_id = 'ENST00000000000'
    elif is_gene_tar:
        print("Running in GENE mode...")
        file_prefix = 'GENE'
        gene_id = 'ENST99999999999'

    # Get the chromosomes that are represented in the SNP/GENE tarball
    # This variants_table file will ONLY have chromosomes in it represented by a provided SNP/GENE list
    chromosomes = set()
    sparse_matrix = csv.DictReader(
        open(tarball_prefix + '.' + file_prefix + '.variants_table.STAAR.tsv', 'r'),
        delimiter="\t",
        quoting=csv.QUOTE_NONE)

    for row in sparse_matrix:
        chromosomes.add(str(row['chrom']))

    # And filter the relevant SAIGE file to just the individuals we want so we can get actual MAC
    cmd = "bcftools view --threads 4 -S /test/SAMPLES_Include.txt -Ob -o /test/" + \
          tarball_prefix + "." + file_prefix + ".saige_input.bcf /test/" + \
          tarball_prefix + "." + file_prefix + ".SAIGE.bcf"
    run_cmd(cmd, True)

    # Build a fake gene_info that can feed into the other functions in this class
    gene_info = pd.Series({'chrom': file_prefix, 'SYMBOL': file_prefix})
    gene_info.name = gene_id

    return gene_info, chromosomes


def process_linear_model_outputs(association_pack: AssociationPack, gene_infos: list) -> list:

    if gene_infos is not None:
        valid_genes = []
        for gene_info in gene_infos:
            valid_genes.append(gene_info.name)
    else:
        valid_genes = None

    if association_pack.is_snp_tar:
        os.rename(association_pack.output_prefix + '.lm_stats.tmp',
                  association_pack.output_prefix + '.SNP.glm.stats.tsv')
        outputs = [association_pack.output_prefix + '.SNP.glm.stats.tsv']
    elif association_pack.is_gene_tar:
        os.rename(association_pack.output_prefix + '.lm_stats.tmp',
                  association_pack.output_prefix + '.GENE.glm.stats.tsv')
        outputs = [association_pack.output_prefix + '.GENE.glm.stats.tsv']
    else:
        # read in the GLM stats file:
        glm_table = pd.read_csv(open(association_pack.output_prefix + ".lm_stats.tmp", 'r'), sep = "\t")

        # Now process the gene table into a useable format:
        # First read in the transcripts file
        transcripts_table = pd.read_csv(gzip.open('transcripts.tsv.gz', 'rt'), sep = "\t")
        transcripts_table = transcripts_table.rename(columns={'#chrom':'chrom'})
        transcripts_table = transcripts_table.set_index('ENST')
        transcripts_table = transcripts_table[transcripts_table['fail'] == False]
        transcripts_table = transcripts_table.drop(columns=['syn.count','fail.cat','fail'])

        # Limit to genes we care about if running only a subset:
        if valid_genes is not None:
            transcripts_table = transcripts_table.loc[valid_genes]

        # Test what columns we have in the 'SNP' field so we can name them...
        field_one = glm_table.iloc[0]
        field_one = field_one['maskname'].split("-")
        field_names = []
        if len(field_one) == 2: # This could be the standard naming format... check that column [1] is MAF/AC
            if 'MAF' in field_one[1] or 'AC' in field_one[1]:
                field_names.extend(['MASK','MAF'])
            else: # This means we didn't hit on MAF/AC in column [2] and a different naming convention is used...
                field_names.extend(['var1','var2'])
        else:
            for i in range(2,len(field_one) + 1):
                field_names.append('var%i' % i)

        # Process the 'SNP' column into separate fields and remove
        glm_table[field_names] = glm_table['maskname'].str.split("-",expand = True)
        glm_table = glm_table.drop(columns=['maskname'])

        # Now merge the transcripts table into the gene table to add annotation and the write
        glm_table = pd.merge(transcripts_table, glm_table, on='ENST', how="left")
        with open(association_pack.output_prefix + '.genes.glm.stats.tsv', 'w') as gene_out:

            # Sort just in case
            glm_table = glm_table.sort_values(by=['chrom','start','end'])

            glm_table.to_csv(path_or_buf=gene_out, index = False, sep="\t", na_rep='NA')
            gene_out.close()

            # And bgzip and tabix...
            cmd = "bgzip /test/" + association_pack.output_prefix + '.genes.glm.stats.tsv'
            run_cmd(cmd, True)
            cmd = "tabix -S 1 -s 2 -b 3 -e 4 /test/" + association_pack.output_prefix + '.genes.glm.stats.tsv.gz'
            run_cmd(cmd, True)

        outputs = [association_pack.output_prefix + '.genes.glm.stats.tsv.gz',
                   association_pack.output_prefix + '.genes.glm.stats.tsv.gz.tbi']

    return outputs


def merge_glm_staar_runs(association_pack: AssociationPack) -> list:

    if association_pack.is_non_standard_tar:

        if association_pack.is_gene_tar:
            prefix = 'GENE'
        elif association_pack.is_snp_tar:
            prefix = 'SNP'

        glm_table = pd.read_csv(association_pack.output_prefix + '.' + prefix + '.glm.stats.tsv', sep='\t')
        staar_table = pd.read_csv(association_pack.output_prefix + '.' + prefix + '.STAAR.stats.tsv', sep='\t')

        # We need to pull the ENST value out of the STAAR table 'SNP' variable while trying to avoid any bugs due to file/pheno names
        field_one = staar_table.iloc[0]
        field_one = field_one['SNP'].split("-")
        field_names = []
        if 'ENST' in field_one[0]:
            field_names.append('ENST')
            for i in range(2, len(field_one) + 1):
                field_names.append('var%i' % i)
        else:  # This means we didn't hit on MAF/AC in column [2] and a different naming convention is used...
            raise dxpy.AppError('ENST value is not in the first column of the STAAR table... Error!')

        staar_table[field_names] = staar_table['SNP'].str.split("-", expand=True)

        # Select STAAR columns we need to merge in/match on
        staar_table = staar_table[['ENST', 'pheno', 'n_var', 'relatedness.correction', 'staar.O.p', 'staar.SKAT.p', 'staar.burden.p', 'staar.ACAT.p']]

        final_table = glm_table.merge(right=staar_table, on=['ENST', 'pheno'])

        with open(association_pack.output_prefix + '.' + prefix + '.STAAR_glm.stats.tsv', 'w') as gene_out:

            final_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')
            gene_out.close()

        outputs = [association_pack.output_prefix + '.' + prefix + '.STAAR_glm.stats.tsv']

    else:
        glm_table = pd.read_csv(gzip.open(association_pack.output_prefix + '.genes.glm.stats.tsv.gz', 'rb'), sep='\t')
        staar_table = pd.read_csv(gzip.open(association_pack.output_prefix + '.genes.STAAR.stats.tsv.gz', 'rb'), sep='\t')

        # Select STAAR columns we need to merge in/match on
        staar_table = staar_table[['ENST', 'MASK', 'MAF', 'pheno', 'n_var', 'relatedness.correction', 'staar.O.p', 'staar.SKAT.p', 'staar.burden.p', 'staar.ACAT.p']]

        final_table = glm_table.merge(right = staar_table, on = ['ENST','MASK','MAF','pheno'])

        with open(association_pack.output_prefix + '.genes.STAAR_glm.stats.tsv', 'w') as gene_out:

            # Sort just in case
            final_table = final_table.sort_values(by=['chrom', 'start', 'end'])

            final_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')
            gene_out.close()

            # And bgzip and tabix...
            cmd = "bgzip /test/" + association_pack.output_prefix + '.genes.STAAR_glm.stats.tsv'
            run_cmd(cmd, True)
            cmd = "tabix -S 1 -s 2 -b 3 -e 4 /test/" + association_pack.output_prefix + '.genes.STAAR_glm.stats.tsv.gz'
            run_cmd(cmd, True)

        outputs = [association_pack.output_prefix + '.genes.STAAR_glm.stats.tsv.gz',
                   association_pack.output_prefix + '.genes.STAAR_glm.stats.tsv.gz.tbi']

    return outputs