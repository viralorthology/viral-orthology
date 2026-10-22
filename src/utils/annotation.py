import re
from collections import Counter

from models.ortholog_group import OrthologGroup

_ANNOTATION_INVALID_CHARS = re.compile(r"[^a-zA-Z0-9]+")


def get_og_most_frequent_annotation(og: OrthologGroup) -> str:
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
