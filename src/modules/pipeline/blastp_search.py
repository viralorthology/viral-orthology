import logging
from pathlib import Path

import engines
import utils
from config.context import Context
from engines.blast import BlastHit
from models.fasta import Fasta
from models.ortholog_group import OrthologGroup
from models.seq import Seq
from models.tmp_fasta import TmpFasta

logger = logging.getLogger(__name__)


def blastp_search(ctx: Context, unique_prots_fasta_path: Path) -> None:
    """
    Search for new orthologs using BlastP.

    Searches unique proteins against the existing ortholog groups and
    assigns proteins to groups based on their best BlastP hit. Proteins
    are added only when their genome is not already represented in the
    corresponding ortholog group.
    """
    ctx.ui.show("Searching for new orthologs using BlastP...")

    unique_prots_fasta = Fasta(unique_prots_fasta_path)
    ogs = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)

    search_and_add_prots_with_blastp(
        ogs, unique_prots_fasta, ctx.args.tool_args["blastp"]
    )


def search_and_add_prots_with_blastp(
    ogs: list[OrthologGroup], unique_prots_fasta: Fasta, blastp_params: str
) -> None:
    """
    Build a BlastP database from existing ortholog groups, run BlastP and process hits.

    Creates a temporary protein database from the sequences in the
    ortholog groups, searches the unique proteins against it, and adds
    proteins to the groups identified by their best hits.
    """
    ogs_by_subject_id = {}
    with TmpFasta() as blast_db:
        for og in ogs:
            og_seqs = list(og.seqs)
            blast_db.add_seqs(*og_seqs)
            ogs_by_subject_id.update({seq.id: og for seq in og_seqs})

        engines.make_blast_db(blast_db, "prot")
        sorted_hits = engines.blastp_search(
            unique_prots_fasta, blast_db, f"{blastp_params} -max_target_seqs 1"
        )  # only best hit

    seqs_by_seq_id = {seq.id: seq for seq in unique_prots_fasta.seqs}
    added_seq_ids = _add_proteins_to_ogs(sorted_hits, ogs_by_subject_id, seqs_by_seq_id)

    if added_seq_ids:
        unique_prots_fasta.remove_seqs(*added_seq_ids)


def _add_proteins_to_ogs(
    sorted_hits: list[BlastHit],
    ogs_by_subject_id: dict[str, OrthologGroup],
    seqs_by_seq_id: dict[str, Seq],
) -> set[str]:
    """
    Add proteins to ortholog groups based on their best BlastP hits.

    Proteins are added only when their genome is not already represented
    in the corresponding ortholog group.

    Returns:
        The IDs of the proteins that were added to ortholog groups.
    """
    added_seq_ids = set()

    for hit in sorted_hits:
        seq = seqs_by_seq_id[hit.query_id]
        og = ogs_by_subject_id[hit.subject_id]
        if seq.genome_id not in og.genome_ids:
            assert seq.id not in added_seq_ids
            og.add_seqs(seq)
            added_seq_ids.add(seq.id)
            logger.info("Protein %s added to %s", seq.id, og.path)

    return added_seq_ids
