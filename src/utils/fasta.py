from models.ortholog_group import OrthologGroup


def o_groups_are_compatible(*ortholog_groups: OrthologGroup) -> bool:
    """
    Return whether ortholog groups contain distinct genome IDs.

    Args:
        *ortholog_groups: Ortholog groups to check.

    Raises:
        ValueError: If fewer than two ortholog groups are provided.
    """
    if len(ortholog_groups) < 2:
        raise ValueError("At least two ortholog groups must be provided")

    genome_ids: set[str] = set()

    for og in ortholog_groups:
        og_genome_ids = og.genome_ids
        if genome_ids.intersection(og_genome_ids):
            return False

        genome_ids.update(og_genome_ids)

    return True
