from collections import defaultdict

import engines
import utils
from config.context import Context
from engines.blast import BlastHit
from models.fasta import Fasta


def find_paralogs(ctx: Context) -> None:
    """
    Find and remove paralog protein sequences from each proteome.

    Reciprocal BLASTP hits are used to identify groups of paralog sequences.
    Within each group, the longest sequence is retained in the proteome and
    the remaining sequences are moved to a separate FASTA file.
    """
    ctx.ui.show("Searching for paralogs...")

    ctx.paths.paralogs_dir.mkdir(exist_ok=True)
    proteomes = utils.get_fastas(ctx.paths.sequences_dir, ".proteome")

    for proteome in ctx.ui.progress_bar(proteomes):
        with engines.BlastDB("prot", proteome) as db:
            blastp_hits = engines.blastp_search(
                proteome, db, ctx.args.tool_args["blastp_paralog_search"]
            )

        # find reciprocal hits
        reciprocal_hits = _find_reciprocal_hits(blastp_hits)
        if not reciprocal_hits:
            continue

        # find groups of paralog seqs
        paralog_groups = _find_paralog_groups(reciprocal_hits)

        # remove paralogs from proteome
        for paralog_group_ids in paralog_groups:
            biggest_prot = max(
                (seq for seq in proteome.seqs if seq.id in paralog_group_ids),
                key=lambda seq: len(seq.seq),
            )
            biggest_prot_id = biggest_prot.id

            prot_ids_to_move = [
                seq_id for seq_id in paralog_group_ids if seq_id != biggest_prot_id
            ]
            seqs_to_move = proteome.get_seqs(*prot_ids_to_move)
            proteome.remove_seqs(*prot_ids_to_move)

            paralogs_fasta = Fasta(ctx.paths.paralogs_dir / f"{biggest_prot_id}.fasta")
            paralogs_fasta.add_seqs(*seqs_to_move)


def _find_reciprocal_hits(blastp_hits: list[BlastHit]) -> set[frozenset[str]]:
    """
    Find reciprocal BLASTP hits between protein sequences.

    A reciprocal hit is a pair of sequences where each sequence is a BLASTP
    hit of the other. Self-hits are excluded, and each reciprocal pair is
    represented as a frozenset to avoid duplicate pairs.

    Args:
        blastp_hits: BLASTP hits to analyze.

    Returns:
        A set of frozensets, where each frozenset contains the IDs of a
        reciprocal hit pair.
    """
    hits_by_query_id = defaultdict(set)
    for hit in blastp_hits:
        if hit.query_id == hit.subject_id:
            continue
        hits_by_query_id[hit.query_id].add(hit.subject_id)

    reciprocal_hits = set()
    for query_id, subject_id_set in hits_by_query_id.items():
        for subject_id in subject_id_set:
            if query_id in hits_by_query_id.get(subject_id, set()):
                reciprocal_hits.add(frozenset((query_id, subject_id)))

    return reciprocal_hits


def _find_paralog_groups(reciprocal_hits: set[frozenset[str]]) -> list[set[str]]:
    """
    Group paralog sequences using a union-find algorithm.

    Args:
        reciprocal_hits: Set of reciprocal BLASTP hits represented as
            unordered pairs of gene IDs.

    Returns:
        A list of sets, where each set contains the gene IDs belonging to
        a group of paralog sequences.
    """
    all_genes = {x for s in reciprocal_hits for x in s}
    parent = {x: x for x in all_genes}

    def find(x: str) -> str:
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def unite(gene1: str, gene2: str) -> None:
        root1 = find(gene1)
        root2 = find(gene2)
        if root1 != root2:
            parent[root2] = root1

    for gene1, gene2 in reciprocal_hits:
        unite(gene1, gene2)

    final_groups: dict[str, set[str]] = {}
    for gene in parent:
        root = find(gene)
        final_groups.setdefault(root, set())
        final_groups[root].add(gene)

    return list(final_groups.values())
