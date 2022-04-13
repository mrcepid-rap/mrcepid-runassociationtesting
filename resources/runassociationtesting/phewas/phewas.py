import dxpy

from ..association_pack import AssociationPack
from ..association_resources import *
from ..tool_runners.glm_runner import GLMRunner


class PheWAS:

    def __init__(self, association_pack: AssociationPack):

        self._association_pack = association_pack

        # If running a SNP-list
        gene_ENST_to_run = []
        if self._association_pack.is_snp_tar:
            gene_ENST_to_run.append('ENST00000000000')
        else:
            genes_to_run = self._association_pack.gene_ids
            transcripts_table = build_transcript_table()
            for gene in genes_to_run:
                gene_info = get_gene_id(gene, transcripts_table)
                gene_ENST_to_run.append(gene_info.name)

        glm_run = GLMRunner(association_pack, genes_to_run=gene_ENST_to_run)
        self.outputs = glm_run.outputs