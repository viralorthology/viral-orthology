import logging
from collections import defaultdict

import engines
import utils
from config.context import Context
from engines.hmmer import HMMHit
from models.fasta import Fasta
from models.ortholog_group import OrthologGroup

logger = logging.getLogger(__name__)


def clean_ortholog_groups(ctx: Context) -> None:
    """
    Clean ortholog groups by retaining the best HMM hit for each protein.

    Each ortholog group is searched against the available proteomes. If a
    protein is detected by multiple ortholog group HMMs, it is assigned to
    the group with the best HMM score (lowest E-value). Proteins that are
    not assigned to their original ortholog group are removed from the
    group.
    """
    ctx.ui.show("Cleaning ortholog groups...")

    o_groups = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)
    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)

    all_hits: dict[HMMHit, OrthologGroup] = {}
    for og in ctx.ui.progress_bar(o_groups):
        hits = engines.hmm_search(og, proteomes_fasta)
        all_hits.update({hit: og for hit in hits})

    sorted_hits = sorted(all_hits.items(), key=lambda item: item[0].evalue)

    best_hits: defaultdict[OrthologGroup, set[str]] = defaultdict(set)
    assigned_seqs = set()
    for hit, og in sorted_hits:
        if hit.seq_id not in assigned_seqs:
            assigned_seqs.add(hit.seq_id)
            best_hits[og].add(hit.seq_id)

    # remove non-best hit proteins
    for og, seqs_to_keep in best_hits.items():
        seqs_to_remove = set(og.ids) - seqs_to_keep
        if seqs_to_remove:
            og.remove_seqs(*seqs_to_remove, ctx=ctx)
            logger.info(
                "Proteins %s removed from %s during HMM-based cleaning",
                seqs_to_remove,
                og.path.name,
            )
