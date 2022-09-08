from runassociationtesting.association_resources import *


def process_linear_model_outputs(output_prefix: str, is_snp_tar: bool = False, is_gene_tar: bool = False,
                                 gene_infos: list = None) -> list:

    if gene_infos is not None:
        valid_genes = []
        for gene_info in gene_infos:
            valid_genes.append(gene_info.name)
    else:
        valid_genes = None

    if is_snp_tar:
        os.rename(output_prefix + '.lm_stats.tmp',
                  output_prefix + '.SNP.glm.stats.tsv')
        outputs = [output_prefix + '.SNP.glm.stats.tsv']
    elif is_gene_tar:
        os.rename(output_prefix + '.lm_stats.tmp',
                  output_prefix + '.GENE.glm.stats.tsv')
        outputs = [output_prefix + '.GENE.glm.stats.tsv']
    else:
        # read in the GLM stats file:
        glm_table = pd.read_csv(open(output_prefix + ".lm_stats.tmp", 'r'), sep="\t")

        # Now process the gene table into a useable format:
        # First read in the transcripts file
        transcripts_table = build_transcript_table()

        # Limit to genes we care about if running only a subset:
        if valid_genes is not None:
            transcripts_table = transcripts_table.loc[valid_genes]

        # Test what columns we have in the 'SNP' field so we can name them...
        field_one = glm_table.iloc[0]
        field_one = field_one['maskname'].split("-")
        field_names = []
        if len(field_one) == 2:  # This could be the standard naming format... check that column [1] is MAF/AC
            if 'MAF' in field_one[1] or 'AC' in field_one[1]:
                field_names.extend(['MASK', 'MAF'])
            else:  # This means we didn't hit on MAF/AC in column [2] and a different naming convention is used...
                field_names.extend(['var1', 'var2'])
        else:
            for i in range(2, len(field_one) + 1):
                field_names.append('var%i' % i)

        # Process the 'SNP' column into separate fields and remove
        glm_table[field_names] = glm_table['maskname'].str.split("-", expand=True)
        glm_table = glm_table.drop(columns=['maskname'])

        # Now merge the transcripts table into the gene table to add annotation and the write
        glm_table = pd.merge(transcripts_table, glm_table, on='ENST', how="left")
        with open(output_prefix + '.genes.glm.stats.tsv', 'w') as gene_out:

            # Sort just in case
            glm_table = glm_table.sort_values(by=['chrom', 'start', 'end'])

            glm_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')
            gene_out.close()

            # And bgzip and tabix...
            cmd = "bgzip /test/" + output_prefix + '.genes.glm.stats.tsv'
            run_cmd(cmd, True)
            cmd = "tabix -S 1 -s 2 -b 3 -e 4 /test/" + output_prefix + '.genes.glm.stats.tsv.gz'
            run_cmd(cmd, True)

        outputs = [output_prefix + '.genes.glm.stats.tsv.gz',
                   output_prefix + '.genes.glm.stats.tsv.gz.tbi']

    return outputs


def merge_glm_staar_runs(output_prefix: str, is_snp_tar: bool = False, is_gene_tar: bool = False) -> list:

    if is_gene_tar or is_snp_tar:
        if is_gene_tar:
            prefix = 'GENE'
        elif is_snp_tar:
            prefix = 'SNP'
        else:
            raise dxpy.AppError('There should be no way to see this message (error in "merge_glm_staar_runs()")...')

        glm_table = pd.read_csv(output_prefix + '.' + prefix + '.glm.stats.tsv', sep='\t')
        staar_table = pd.read_csv(output_prefix + '.' + prefix + '.STAAR.stats.tsv', sep='\t')

        # We need to pull the ENST value out of the STAAR table 'SNP' variable while
        # trying to avoid any bugs due to file/pheno names
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
        staar_table = staar_table[['ENST', 'pheno', 'n_var', 'relatedness.correction', 'staar.O.p', 'staar.SKAT.p',
                                   'staar.burden.p', 'staar.ACAT.p']]

        final_table = glm_table.merge(right=staar_table, on=['ENST', 'pheno'])

        with open(output_prefix + '.' + prefix + '.STAAR_glm.stats.tsv', 'w') as gene_out:

            final_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')
            gene_out.close()

        outputs = [output_prefix + '.' + prefix + '.STAAR_glm.stats.tsv']

    else:
        glm_table = pd.read_csv(gzip.open(output_prefix + '.genes.glm.stats.tsv.gz', 'rb'), sep='\t')
        staar_table = pd.read_csv(gzip.open(output_prefix + '.genes.STAAR.stats.tsv.gz', 'rb'), sep='\t')

        # Select STAAR columns we need to merge in/match on
        staar_table = staar_table[['ENST', 'MASK', 'MAF', 'pheno', 'n_var', 'relatedness.correction', 'staar.O.p',
                                   'staar.SKAT.p', 'staar.burden.p', 'staar.ACAT.p']]

        final_table = glm_table.merge(right=staar_table, on=['ENST', 'MASK', 'MAF', 'pheno'])

        with open(output_prefix + '.genes.STAAR_glm.stats.tsv', 'w') as gene_out:

            # Sort just in case
            final_table = final_table.sort_values(by=['chrom', 'start', 'end'])

            final_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')
            gene_out.close()

            # And bgzip and tabix...
            cmd = "bgzip /test/" + output_prefix + '.genes.STAAR_glm.stats.tsv'
            run_cmd(cmd, True)
            cmd = "tabix -S 1 -s 2 -b 3 -e 4 /test/" + output_prefix + '.genes.STAAR_glm.stats.tsv.gz'
            run_cmd(cmd, True)

        outputs = [output_prefix + '.genes.STAAR_glm.stats.tsv.gz',
                   output_prefix + '.genes.STAAR_glm.stats.tsv.gz.tbi']

    return outputs
