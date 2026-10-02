from collections.abc import Iterator

import engines
import utils
from config.context import Context
from engines.hmmer import HMMHit
from models.fasta import Fasta
from models.ortholog_group import OrthologGroup
from models.seq import Seq
from models.tmp_fasta import TmpFasta


def paralogs_hmm_search(ctx: Context) -> None:
    """Find and assign paralog sequences to ortholog groups using HMM searches."""
    ctx.ui.show("Checking the paralog sequences...")

    ogs = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)
    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)
    paralogs_by_seq_id = {
        seq.id: seq for seq in proteomes_fasta.seqs if seq.id in ctx.runtime.paralog_ids
    }

    sorted_hits = _get_hmm_sorted_hits(ogs, paralogs_by_seq_id)

    added_paralogs = set()
    for hit, og in sorted_hits.items():
        new_paralog = paralogs_by_seq_id[hit.seq_id]

        if new_paralog.id in added_paralogs:
            continue

        # add new paralog to og directly
        if new_paralog.genome_id not in og.genome_ids:
            og.add_seqs(new_paralog)
            ctx.runtime.paralog_ids.remove(new_paralog.id)
            added_paralogs.add(new_paralog.id)
            continue

        # compare the paralog in the og with the new paralog
        old_paralog, other_og_proteins = _filter_og_proteins(
            og.seqs, new_paralog.genome_id
        )
        best_paralog_id = _find_best_paralog(
            new_paralog,
            old_paralog,
            other_og_proteins,
        )

        # remove old paralog from og and add the new one
        if best_paralog_id != old_paralog.id:
            tmp_og = Fasta(og.path)  # tmp possibly invalid og
            tmp_og.remove_seqs(old_paralog.id)
            tmp_og.add_seqs(new_paralog)

            ctx.runtime.paralog_ids.remove(new_paralog.id)
            ctx.runtime.paralog_ids.add(old_paralog.id)
            added_paralogs.add(new_paralog.id)


def _get_hmm_sorted_hits(
    ogs: list[OrthologGroup],
    paralogs_by_seq_id: dict[str, Seq],
) -> dict[HMMHit, OrthologGroup]:
    """
    Find HMM hits for paralog sequences and sort them by e-value.

    Args:
        ogs: Ortholog groups to search against.
        paralogs_by_seq_id: Paralog sequences indexed by sequence ID.

    Returns:
        A dictionary mapping each HMM hit to its corresponding ortholog group,
        sorted by increasing e-value.
    """
    all_hits: dict[HMMHit, OrthologGroup] = {}

    with TmpFasta() as paralogs_fasta:
        paralogs_fasta.add_seqs(*paralogs_by_seq_id.values())
        for og in ogs:
            hits = engines.hmm_search(og, paralogs_fasta)
            all_hits.update({hit: og for hit in hits})

    return dict(sorted(all_hits.items(), key=lambda item: item[0].evalue))


def _filter_og_proteins(
    og_proteins: Iterator[Seq], paralogs_genome_id: str
) -> tuple[Seq, list[Seq]]:
    """
    Separate the existing paralog from the other proteins in an ortholog group.

    Args:
        og_proteins: Proteins belonging to an ortholog group.
        paralogs_genome_id: Genome ID of the paralogs being evaluated.

    Returns:
        A tuple containing the existing paralog and the proteins from other genomes.
    """
    old_paralog = None
    other_og_proteins = []

    for seq in og_proteins:
        if seq.genome_id == paralogs_genome_id:
            assert old_paralog is None
            old_paralog = seq
        else:
            other_og_proteins.append(seq)

    assert old_paralog is not None
    assert other_og_proteins
    return old_paralog, other_og_proteins


def _find_best_paralog(
    new_paralog: Seq, old_paralog: Seq, other_og_proteins: list[Seq]
) -> str:
    """
    Find the best paralog using Blastp.

    Args:
        new_paralog: New paralog candidate being evaluated.
        old_paralog: Existing paralog in the ortholog group.
        other_og_proteins: Proteins from the ortholog group belonging to other
            genomes.
    """
    with TmpFasta() as og_blast_db:
        og_blast_db.add_seqs(*other_og_proteins)
        engines.make_blast_db(og_blast_db, "prot")

        with TmpFasta() as query_fasta:
            query_fasta.add_seqs(new_paralog, old_paralog)

            hits = engines.blastp_search(
                query_fasta,
                og_blast_db,
                "-evalue 50 -dbsize 100000000 -word_size 2",
            )
            assert hits

    return hits[0].query_id
