import logging

import utils
from config.context import Context
from models.ortholog_group import OrthologGroup
from modules.pipeline.rename_og_fastas_by_annotation import (
    get_og_most_frequent_annotation,
)

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

    og_annotation = {og: get_og_most_frequent_annotation(og) for og in ortholog_groups}
    og_genomes = {og: set(og.genome_ids) for og in ortholog_groups}

    active_genomes_count = len(ctx.runtime.active_genome_ids)
    for i, og1 in enumerate(ortholog_groups):
        if _skip_og(og1, active_genomes_count):
            continue

        for og2 in ortholog_groups[i + 1 :]:
            if _skip_og(og2, active_genomes_count):
                continue

            if og_annotation[og1] != og_annotation[og2]:
                continue

            if og_genomes[og1].intersection(og_genomes[og2]):
                continue

            if ctx.args.assume_yes or ctx.ui.ask_yes_no(
                ("\n").join(
                    [
                        "The following files can be merged:",
                        f"{og1.path.name} / {og2.path.name}",
                        "Want to merge them?",
                    ]
                )
            ):
                og1.add_seqs(*og2.seqs)
                og_annotation[og1] = get_og_most_frequent_annotation(og1)
                og_genomes[og1] = set(og1.genome_ids)
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
