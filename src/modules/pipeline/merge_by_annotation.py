import logging

import utils
from config.context import Context
from modules.pipeline.rename_og_fastas_by_annotation import (
    get_og_most_frequent_annotation,
)

logger = logging.getLogger(__name__)


def merge_by_annotation(ctx: Context) -> None:
    """
    Merge compatible ortholog groups that share the same annotation.

    Ortholog groups annotated as hypothetical protein are skipped. Compatible
    groups sharing annotation are merged into the first group, and the merged
    group's FASTA file is deleted.
    """
    ogs = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)

    _og_genomes_id_list = {og: og.genome_ids for og in ogs}

    og_n_seqs = {og: len(_og_genomes_id_list[og]) for og in ogs}
    og_genomes = {og: set(_og_genomes_id_list[og]) for og in ogs}
    og_annotation = {og: get_og_most_frequent_annotation(og) for og in ogs}
    active_genomes_count = len(ctx.runtime.active_genome_ids)

    merged_ogs = set()
    for i, og1 in enumerate(ogs):
        if og1 in merged_ogs:
            continue
        if _skip_og(og_annotation[og1], og_n_seqs[og1], active_genomes_count):
            continue

        for og2 in ogs[i + 1 :]:
            if og2 in merged_ogs:
                continue
            if _skip_og(og_annotation[og2], og_n_seqs[og2], active_genomes_count):
                continue

            if _ogs_can_be_merged(
                og_annotation[og1], og_annotation[og2], og_genomes[og1], og_genomes[og2]
            ) and (
                ctx.args.assume_yes
                or ctx.ui.ask_yes_no(
                    ("\n").join(
                        [
                            "The following files can be merged:",
                            f"{og1.path.name} / {og2.path.name}",
                            "Want to merge them?",
                        ]
                    )
                )
            ):
                og1.add_seqs(*og2.seqs)
                og_annotation[og1] = get_og_most_frequent_annotation(og1)
                og_n_seqs[og1] = og1.n_seqs
                og_genomes[og1] = set(og1.genome_ids)
                og2.delete_fasta()
                merged_ogs.add(og2)

                logger.info(
                    "Ortholog groups merged by annotation: %s, %s",
                    og1.path.name,
                    og2.path.name,
                )


def _skip_og(og_annotation: str, og_n_seqs: int, active_genomes_count: int) -> bool:
    """
    Determine whether an ortholog group should be skipped.

    Args:
        og_annotation: Annotation assigned to the ortholog group.
        og_n_seqs: Number of sequences in the ortholog group.
        active_genomes_count: Number of active genomes in the dataset.

    Returns:
        True if the ortholog group is annotated as a hypothetical protein or
        contains one sequence per active genome; otherwise, False.
    """
    return og_annotation == "hypothetical-protein" or og_n_seqs == active_genomes_count


def _ogs_can_be_merged(
    og1_annotation: str,
    og2_annotation: str,
    og1_genome_ids: set[str],
    og2_genome_ids: set[str],
) -> bool:
    """
    Determine whether two ortholog groups can be merged.

    Args:
        og1_annotation: Annotation assigned to the first ortholog group.
        og2_annotation: Annotation assigned to the second ortholog group.
        og1_genome_ids: Genome IDs represented in the first ortholog group.
        og2_genome_ids: Genome IDs represented in the second ortholog group.

    Returns:
        True if both ortholog groups have the same annotation and do not
        contain sequences from any of the same genomes; otherwise, False.
    """
    return og1_annotation == og2_annotation and not og1_genome_ids.intersection(
        og2_genome_ids
    )
