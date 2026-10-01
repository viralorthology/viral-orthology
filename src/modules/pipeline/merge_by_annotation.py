import logging

import utils
from config.context import Context
from models.ortholog_group import OrthologGroup

logger = logging.getLogger(__name__)


def merge_by_annotation(ctx: Context) -> None:
    """
    Merge compatible ortholog groups that share the same annotation.

    Ortholog groups representing hypothetical proteins or already containing
    one sequence per active genome are skipped. Compatible groups sharing an
    annotation are merged into the first group, and the merged group's FASTA
    file is removed.
    """
    ortholog_groups = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)
    ortholog_groups.sort(key=lambda og: og.path.stem)

    active_genomes_count = len(ctx.runtime.active_genome_ids)
    for i, og1 in enumerate(ortholog_groups):
        if _skip_og(og1, active_genomes_count):
            continue

        for og2 in ortholog_groups[i + 1 :]:
            if _skip_og(og2, active_genomes_count):
                continue

            if not _ogs_have_same_annotation(og1.path.stem, og2.path.stem):
                continue

            if not utils.o_groups_are_compatible(og1, og2):
                continue

            if ctx.args.assume_yes or ctx.ui.ask_yes_no(
                ("\n").join(
                    [
                        "The following files can be merged:",
                        f"{og1.path} / {og2.path}",
                        "Want to merge them?",
                    ]
                )
            ):
                og1.add_seqs(*og2.seqs)
                og2.delete_fasta()

                logger.info(
                    "Ortholog groups merged by annotation: %s, %s",
                    og1.path,
                    og2.path,
                )


def _skip_og(og: OrthologGroup, active_genomes_count: int) -> bool:
    """
    Return whether an ortholog group should be skipped during merging.

    Skips groups whose file does not exist, whose annotation is
    hypothetical-protein, or that already contain all active genomes.
    """
    return (
        not og.path.exists()
        or "hypothetical-protein" in og.path.stem
        or og.n_seqs == active_genomes_count
    )


def _ogs_have_same_annotation(og1_stem: str, og2_stem: str) -> bool:
    """
    Return whether two ortholog groups represent the same annotation.

    Ortholog groups with a numeric suffix in their name are considered
    to have the same annotation as the unsuffixed group. For example,
    'pol', 'pol-1', and 'pol-2' are considered the same annotation.
    """

    def _annotation(stem: str) -> str:
        parts = stem.rsplit("-", 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0]
        return stem

    return _annotation(og1_stem) == _annotation(og2_stem)
