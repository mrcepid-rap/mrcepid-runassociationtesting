import os
import dxpy
import gzip
import subprocess
import pandas as pd
import pandas.core.series


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
        print("STDERROR follows\n")
        print(stderr.decode('utf-8'))
        raise dxpy.AppError("Failed to run properly...")


# This is to generate a global CHROMOSOMES variable for parallelisation
def get_chromosomes(is_snp_tar: bool = False):

    if is_snp_tar:
        chromosomes = list(['SNP'])
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
        fixed_samp_saige = open(chromosome + '.markers.saige.sample', 'w')
        for line in samp_file:
            line = line.rstrip().split(" ")
            if line[0] == 'ID_1':
                fixed_samp_bolt.write('ID_1 ID_2 missing sex\n')
            elif line[3] == 'D':
                fixed_samp_bolt.write('0 0 0 D\n')
            else:
                fixed_samp_bolt.write("%s %s 0 NA\n" % (line[1], line[1]))
                fixed_samp_saige.write("%s 0 NA\n" % (line[1]))
        samp_file.close()
        fixed_samp_bolt.close()
        fixed_samp_saige.close()

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
        print("gene_id – " + gene_id + " – looks like an ENST value... validating...")
        try:
            gene_info = transcripts_table.loc[gene_id]
            print("Found one matching ENST (%s - %s)... proceeding..." % (gene_id, gene_info['coord']))
        except KeyError:
            print("Did not find a transcript with ENST value %s... terminating..." % gene_id)
    else:
        print("gene_id – " + gene_id + " – does not look like an ENST value, searching for symbol instead...")
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

    return gene_info