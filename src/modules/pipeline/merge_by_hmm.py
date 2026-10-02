import engines
import utils
from config.context import Context
from models.ortholog_group import OrthologGroup
from models.tmp_fasta import TmpFasta

MIN_SCORE_TO_MERGE_HMMS = 70


def merge_by_hmm(ctx: Context) -> None:
    """
    Merge compatible ortholog groups based on HMM similarity.
    """
    ctx.ui.show("Merging ortholog groups using HMMs...")

    ogs = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)

    merged_ogs = set()
    for i, og1 in enumerate(ctx.ui.progress_bar(ogs)):
        if og1 in merged_ogs:
            continue

        ogs_to_try = _get_ogs_to_try(
            og1, [og for og in ogs[i + 1 :] if og not in merged_ogs]
        )

        for og2 in ogs_to_try:
            if not utils.o_groups_are_compatible(og1, og2):
                continue

            score = engines.hmm_compare_groups(og1, og2)

            if score >= MIN_SCORE_TO_MERGE_HMMS:
                og1.add_seqs(*og2.seqs)
                merged_ogs.add(og2)

    for og in merged_ogs:
        og.delete_fasta()


def _get_ogs_to_try(
    og1: OrthologGroup, available_ogs: list[OrthologGroup]
) -> list[OrthologGroup]:
    """
    Find ortholog groups containing sequences that match the given ortholog group in a BLASTP search.

    Args:
        og1: Ortholog group used as the BLASTP query.
        available_ogs: Ortholog groups available as BLASTP subjects.

    Returns:
        Ortholog groups containing sequences that match `og1` in the BLASTP search.
    """
    if not available_ogs:
        return []

    ogs_by_seq_id = {}

    with TmpFasta() as blastp_db:
        for og in available_ogs:
            ogs_by_seq_id.update({seq.id: og for seq in og.seqs})
            blastp_db.add_seqs(*og.seqs)

        engines.make_blast_db(blastp_db, "prot")
        hits = engines.blastp_search(
            og1, blastp_db, "-evalue 50 -dbsize 100000000 -word_size 2"
        )

    ogs_to_try = []
    for hit in hits:
        og = ogs_by_seq_id[hit.subject_id]
        if og not in ogs_to_try:
            ogs_to_try.append(og)

    return ogs_to_try
