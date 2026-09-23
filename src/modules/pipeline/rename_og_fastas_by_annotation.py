import re
from collections import Counter

import utils
from config.context import Context
from models.ortholog_group import OrthologGroup

_ANNOTATION_INVALID_CHARS = re.compile(r"[^a-zA-Z0-9]+")


def rename_og_fastas_by_annotation(ctx: Context) -> None:
    """
    Rename ortholog group FASTA files using their most frequent annotation.

    Temporary names are assigned first to avoid filename collisions during
    the renaming process. Repeated annotations receive a numeric suffix.
    """
    og_fastas = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)

    # give tmp names
    for n, og in enumerate(og_fastas):
        og.rename_fasta(f"{n}_tmp.fasta")

    # rename files
    used_names: set[str] = set()
    for og in og_fastas:
        new_name = _get_og_most_frequent_annotation(og)

        if new_name in used_names:
            new_name = _add_number_to_repeated_name(new_name, used_names)

        og.rename_fasta(f"{new_name}.fasta")
        used_names.add(new_name)


def _get_og_most_frequent_annotation(og: OrthologGroup) -> str:
    """
    Get the most frequent valid annotation in an ortholog group.

    An annotation is used only if it occurs more than once. If no valid
    annotation meets this criterion, "hypothetical-protein" is returned.
    """
    og_annotations = [_get_seq_annotation(seq.description) for seq in og.seqs]
    filtered_annotations = _filter_annotations(og_annotations)

    if filtered_annotations:
        annotation, count = _get_most_frequent_annotation(filtered_annotations)

        if count > 1 and annotation:
            return annotation

    return "hypothetical-protein"


def _filter_annotations(og_annotations: list[str | None]) -> list[str]:
    """Filter out missing, unknown, and hypothetical annotations."""
    filtered_annotations = []

    for seq_annotation in og_annotations:
        if (
            seq_annotation is None
            or "unknown" in seq_annotation
            or "hypothetical" in seq_annotation
        ):
            continue
        filtered_annotations.append(seq_annotation)

    return filtered_annotations


def _get_seq_annotation(seq_description: str) -> str | None:
    """Extract and sanitize a protein or gene annotation from a sequence description."""
    for tag in ("[protein=", "[gene="):
        if tag in seq_description:
            annotation = seq_description.split(tag)[1].split("]")[0].lower()
            return _sanitize_annotation(annotation)

    return None


def _sanitize_annotation(seq_annotation: str) -> str:
    """Replace non-alphanumeric characters with hyphens and strip leading and trailing hyphens."""
    return _ANNOTATION_INVALID_CHARS.sub("-", seq_annotation).strip("-")


def _get_most_frequent_annotation(annotations: list[str]) -> tuple[str, int]:
    """
    Get the most frequent annotation and its count.

    Ties are broken alphabetically to ensure deterministic results.
    """
    annotation, count = min(
        Counter(annotations).items(),
        key=lambda item: (
            -item[1],
            item[0],
        ),  # Sort by frequency, then alphabetically (ensure reproducibility)
    )

    return annotation, count


def _add_number_to_repeated_name(annotation: str, used_names: set[str]) -> str:
    """Generate a unique name by appending an incrementing numeric suffix."""
    n = 0
    while True:
        n += 1
        possible_name = f"{annotation}-{n}"

        if possible_name not in used_names:
            return possible_name
