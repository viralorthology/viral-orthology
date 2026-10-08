import logging

import engines
import utils
from cli.ui import UI
from config.context import Context
from engines.hmmer import HMMHit
from models.fasta import Fasta
from models.ortholog_group import OrthologGroup
from models.seq import Seq

logger = logging.getLogger(__name__)


def hmm_search(ctx: Context) -> None:
    """
    Search for new orthologs using HMM profiles.

    Searches for orthologous proteins among annotated and predicted unique
    proteins using HMM profiles associated with the existing ortholog groups.
    Annotated proteins are searched iteratively until no additional proteins
    can be assigned, while predicted proteins are processed in a single round.
    """
    ctx.ui.show("Searching for new orthologs using HMMs...")

    annotated_unique_prots_fasta = Fasta(ctx.paths.annotated_unique_prots_fasta)
    predicted_unique_prots_fasta = Fasta(ctx.paths.predicted_unique_prots_fasta)

    ogs = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)
    active_genomes_count = len(ctx.runtime.active_genome_ids)

    while True:
        while True:
            annotated_prots_added = search_and_add_prots_with_hmm(
                ctx.ui,
                [og for og in ogs if og.n_seqs < active_genomes_count],
                annotated_unique_prots_fasta,
                ctx.args.tool_args["hmmsearch"],
            )

            if not annotated_prots_added:
                break

        predicted_prots_added = search_and_add_prots_with_hmm(
            ctx.ui,
            [og for og in ogs if og.n_seqs < active_genomes_count],
            predicted_unique_prots_fasta,
            ctx.args.tool_args["hmmsearch"],
        )

        if not predicted_prots_added:
            break


def search_and_add_prots_with_hmm(
    ui: UI,
    ogs: list[OrthologGroup],
    unique_prots_fasta: Fasta,
    params: str,
) -> bool:
    """
    Assign unique proteins to ortholog groups based on HMM search hits.

    Searches each incomplete ortholog group against the provided protein
    sequences and adds each protein to the ortholog group corresponding
    to its best HMM hit. Added proteins are removed from the unique proteins FASTA.

    Args:
        ui: User interface used to display progress.
        ogs: Ortholog groups to search against.
        unique_prots_fasta: FASTA file containing candidate protein sequences.
        params: Parameters passed to the HMM search tool.

    Returns:
        True if at least one protein was added to an ortholog group,
        otherwise False.
    """
    all_hits: dict[HMMHit, OrthologGroup] = {}

    for og in ui.progress_bar(ogs):
        hits = engines.hmm_search(og, unique_prots_fasta, params)
        all_hits.update({hit: og for hit in hits})

    sorted_hits = dict(sorted(all_hits.items(), key=lambda item: item[0].evalue))

    seqs_by_seq_id = {seq.id: seq for seq in unique_prots_fasta.seqs}
    added_gene_ids = _add_proteins_to_ogs(sorted_hits, seqs_by_seq_id)

    if not added_gene_ids:
        return False

    unique_prots_fasta.remove_seqs(*added_gene_ids)

    return True


def _add_proteins_to_ogs(
    sorted_hits: dict[HMMHit, OrthologGroup],
    seqs_by_seq_id: dict[str, Seq],
) -> set[str]:
    """
    Add proteins to ortholog groups based on their best HMM hits.

    Processes HMM hits in ascending E-value order so that each protein is
    assigned to the ortholog group associated with its best hit. A protein
    is not added if its genome is already represented in the target
    ortholog group.

    Returns:
        Set of protein IDs that were successfully added to ortholog groups.
    """
    added_gene_ids = set()
    for hit, og in sorted_hits.items():
        if hit.seq_id in added_gene_ids:
            continue
        if hit.genome_id in og.genome_ids:
            continue

        seq = seqs_by_seq_id[hit.seq_id]
        og.add_seqs(seq)
        added_gene_ids.add(seq.id)

        logger.info(
            "Protein %s added to %s based on HMM profile search, with an evalue of %s",
            seq.id,
            og.path.name,
            hit.evalue,
        )

    return added_gene_ids
