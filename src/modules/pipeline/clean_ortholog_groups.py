import logging

import engines
import utils
from config.context import Context
from models.fasta import Fasta

logger = logging.getLogger(__name__)


def clean_ortholog_groups(ctx: Context) -> None:
    """
    Remove sequences from ortholog groups that are not detected by HMM search.

    Each ortholog group is searched against the available proteomes, and sequences
    whose IDs are not present among the search hits are removed from the group.
    """
    # TODO remove protein if its og is not the best hit?
    ctx.ui.show("Cleaning ortholog groups...")

    o_groups = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)
    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)
    for og in ctx.ui.progress_bar(o_groups):
        hits = engines.hmm_search(og, proteomes_fasta)
        id_hits = {hit.seq_id for hit in hits}

        ids_to_remove = [id_ for id_ in og.ids if id_ not in id_hits]

        if ids_to_remove:
            og.remove_seqs(*ids_to_remove, ctx=ctx)
            logger.debug("proteins %s removed from %", ids_to_remove, og.path)
