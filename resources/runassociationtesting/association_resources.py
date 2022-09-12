import os
import csv
import dxpy
import gzip
import subprocess
import pandas as pd
import pandas.core.series


# This function runs a command on an instance, either with or without calling the docker instance we downloaded
# By default, commands are not run via Docker, but can be changed by setting is_docker = True
# Also, by default, standard out is not saved, but can be modified with the 'stdout_file' parameter.
# print_cmd is for internal debugging purposes when testing new code
def run_cmd(cmd: str, is_docker: bool = False, stdout_file: str = None, print_cmd=False) -> None:

    # -v here mounts a local directory on an instance (in this case the home dir) to a directory internal to the
    # Docker instance named /test/. This allows us to run commands on files stored on the AWS instance within Docker.
    # This looks slightly different from other versions of this command I have written as I needed to write a custom
    # R script to run STAAR. That means we have multiple mounts here to enable this code to find the script.
    if is_docker:
        cmd = "docker run " \
              "-v /home/dnanexus:/test " \
              "-v /usr/bin/:/prog " \
              "egardner413/mrcepid-burdentesting " + cmd

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
        chromosomes = list(range(1, 23))  # Is 0 based on the right coordinate...? (So does 1..22)
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
    dxpy.download_dxfile(bgen_index.get_id(), f'filtered_bgen/{chromosome}.filtered.bgen.bgi')
    dxpy.download_dxfile(bgen_sample.get_id(), f'filtered_bgen/{chromosome}.filtered.sample')
    dxpy.download_dxfile(bgen.get_id(), f'filtered_bgen/{chromosome}.filtered.bgen')
    dxpy.download_dxfile(vep.get_id(), f'filtered_bgen/{chromosome}.filtered.vep.tsv.gz')

    # And then do the filtering...
    cmd = f'plink2 --threads 4 --bgen /test/filtered_bgen/{chromosome}.filtered.bgen "ref-last" ' \
          f'--sample /test/filtered_bgen/{chromosome}.filtered.sample ' \
          f'--export bgen-1.2 "bits="8 ' \
          f'--out /test/{chromosome}.markers ' \
          f'--keep /test/SAMPLES_Include.txt'
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
                fixed_samp_bolt.write(f'{line[1]} {line[1]} 0 NA\n')
        samp_file.close()
        fixed_samp_bolt.close()

    # And finally index the file
    cmd = f'bgenix -index -g /test/{chromosome}.markers.bgen'
    run_cmd(cmd, True)


# Build the pandas DataFrame of transcripts
def build_transcript_table() -> pandas.DataFrame:

    transcripts_table = pd.read_csv(gzip.open('transcripts.tsv.gz', 'rt'), sep="\t")
    transcripts_table = transcripts_table.set_index('ENST')
    transcripts_table = transcripts_table[transcripts_table['fail'] == False]
    transcripts_table = transcripts_table.drop(columns=['syn.count', 'fail.cat', 'fail'])
    return transcripts_table


def get_gene_id(gene_id: str, transcripts_table: pandas.DataFrame) -> pandas.core.series.Series:

    if 'ENST' in gene_id:
        print("gene_id – " + gene_id + " – looks like an ENST value... validating...")
        try:
            gene_info = transcripts_table.loc[gene_id]
            print(f'Found one matching ENST ({gene_id} - {gene_info["coord"]})... proceeding...')
        except KeyError:
            raise dxpy.AppError(f'Did not find a transcript with ENST value {gene_id}... terminating...')
    else:
        print("gene_id – " + gene_id + " – does not look like an ENST value, searching for symbol instead...")
        found_rows = transcripts_table[transcripts_table['SYMBOL'] == gene_id]
        if len(found_rows) == 1:
            found_enst = found_rows.index[0]
            gene_info = transcripts_table.loc[found_enst]
            print(f'Found one matching ENST ({found_enst} - {gene_info["coord"]}) for SYMBOL {gene_id}... '
                  f'proceeding...')
        elif len(found_rows) > 1:
            raise dxpy.AppError(f'Found {len(found_rows)} ENST IDs ({",".join(found_rows.index.to_list())} for SYMBOL '
                                f'{gene_id}... Please re-run using exact ENST to ensure consistent results...')
        else:
            raise dxpy.AppError(f'Did not find an associated ENST ID for SYMBOL {gene_id}... '
                                f'Please re-run after checking SYMBOL/ENST used...')

    return gene_info


def process_snp_or_gene_tar(is_snp_tar, is_gene_tar, tarball_prefix) -> tuple:

    if is_snp_tar:
        print("Running in SNP mode...")
        file_prefix = 'SNP'
        gene_id = 'ENST00000000000'
    elif is_gene_tar:
        print("Running in GENE mode...")
        file_prefix = 'GENE'
        gene_id = 'ENST99999999999'
    else:
        raise dxpy.AppError('There is no way you should see this error (process_snp_or_gene_tar)')

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
    cmd = f'bcftools view --threads 4 -S /test/SAMPLES_Include.txt -Ob -o /test/' \
          f'{tarball_prefix}.{file_prefix}.saige_input.bcf /test/' \
          f'{tarball_prefix}.{file_prefix}.SAIGE.bcf'
    run_cmd(cmd, True)

    # Build a fake gene_info that can feed into the other functions in this class
    gene_info = pd.Series({'chrom': file_prefix, 'SYMBOL': file_prefix})
    gene_info.name = gene_id

    return gene_info, chromosomes
