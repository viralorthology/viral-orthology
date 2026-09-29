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
        if (
            "hypothetical-protein" in og1.path.stem
            or og1.n_seqs == active_genomes_count
        ):
            continue

        ogs_same_annotation = _get_same_annotation_ogs(
            og1.path.stem,
            ortholog_groups[i + 1 :],
        )

        for og2 in ogs_same_annotation:
            if utils.o_groups_are_compatible(og1, og2) and (
                ctx.args.assume_yes
                or ctx.ui.ask_yes_no(
                    f"The following files can be merged:\n{og1.path} / {og2.path}\nWant to merge them?"
                )
            ):
                og1.add_seqs(*og2.seqs)
                og2.delete_fasta()

                logger.info(
                    "Ortholog groups merged by annotation: %s, %s",
                    og1.path,
                    og2.path,
                )


def _get_same_annotation_ogs(
    og1_stem: str, ogs: list[OrthologGroup]
) -> list[OrthologGroup]:
    """
    Return ortholog groups with the same annotation as ``og1_stem``.

    The annotation is inferred from the ortholog group filename by removing
    the suffix after the last hyphen. Only groups whose FASTA file still
    exists are returned.

    Returns:
        Ortholog groups whose annotation matches ``og1_stem`` and whose
        FASTA files still exist.
    """
    same_annotation_ogs = []

    for og2 in ogs:
        assert og2.path.stem != og1_stem

        if og2.path.exists() and og2.path.stem.rsplit("-", 1)[0] == og1_stem:
            same_annotation_ogs.append(og2)

    return same_annotation_ogs
