import utils
from config.context import Context
from models.fasta import Fasta


def make_annotated_unique_prots_fasta(ctx: Context) -> None:
    """
    Create the annotated unique proteins FASTA file.

    Selects protein sequences from the active genomes that are not identified
    as paralogs and do not belong to any orthologous group. The selected
    sequences are written to the annotated unique proteins FASTA file.

    Args:
        ctx: Pipeline context containing input/output paths, active genome IDs,
            paralog IDs, and other runtime information.
    """
    ogs = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)
    seq_in_group_ids = set()
    for og in ogs:
        seq_in_group_ids.update(og.ids)

    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)
    annotated_unique_prots_fasta = Fasta(ctx.paths.annotated_unique_prots_fasta)

    seqs = []
    for seq in proteomes_fasta.seqs:
        if seq.genome_id not in ctx.runtime.active_genome_ids:
            continue
        if seq.id not in ctx.runtime.paralog_ids and seq.id not in seq_in_group_ids:
            seqs.append(seq)

    annotated_unique_prots_fasta.add_seqs(*seqs)
