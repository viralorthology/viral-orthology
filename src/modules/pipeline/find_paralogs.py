import logging
from collections import defaultdict

import engines
from config.context import Context
from engines.blast import BlastHit
from models.fasta import Fasta
from models.tmp_fasta import TmpFasta

logger = logging.getLogger(__name__)


def find_paralogs(ctx: Context) -> None:
    """
    Identify paralog proteins within each proteome and select one representative.

    BLASTP reciprocal hits are used to identify groups of paralogous proteins.
    Within each group, the longest protein sequence is selected for further
    analysis, while the remaining paralogs are marked as ignored in the
    runtime context.
    """
    ctx.ui.show("Searching for paralogs...")

    genomes_fasta = Fasta(ctx.paths.genomes_fasta)
    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)

    for genome_id in ctx.ui.progress_bar(genomes_fasta.ids):
        proteome = [seq for seq in proteomes_fasta.seqs if seq.genome_id == genome_id]
        if not proteome:
            logger.info("Genome %s ignored: no annotated proteome found", genome_id)
            continue

        with TmpFasta() as proteome_fasta:
            proteome_fasta.add_seqs(*proteome)
            engines.make_blast_db(proteome_fasta, "prot")

            blastp_hits = engines.blastp_search(
                proteome_fasta,
                proteome_fasta,
                ctx.args.tool_args["blastp_paralog_search"],
            )
            logger.debug("blastp hits for %s: %s", genome_id, blastp_hits)

        # find reciprocal hits
        reciprocal_hits = _find_reciprocal_hits(blastp_hits)
        if not reciprocal_hits:
            logger.info("No paralogs found in genome %s", genome_id)
            continue

        # find groups of paralog seqs
        paralog_groups = _find_paralog_groups(reciprocal_hits)
        logger.info("Paralogs found in genome %s: %s", genome_id, paralog_groups)

        # add ignored paralogs to Runtime
        for paralog_group_ids in paralog_groups:
            selected_paralog = max(
                (seq for seq in proteome if seq.id in paralog_group_ids),
                key=lambda seq: (len(seq.seq), seq.id),  # Ensure reproducibility
            )

            ignored_paralog_ids = [
                seq_id for seq_id in paralog_group_ids if seq_id != selected_paralog.id
            ]

            ctx.runtime.paralog_ids.update(ignored_paralog_ids)

            logger.info(
                "Genome %s proteins ignored: %s, protein %s selected for further analysis",
                genome_id,
                ignored_paralog_ids,
                selected_paralog.id,
            )


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
