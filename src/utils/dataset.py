import hashlib


def get_dataset_hash(genome_ids: list[str]) -> str:
    """
    Return a short deterministic hash for a collection of genome IDs.

    The genome IDs are sorted before hashing, so the resulting hash is independent of their input order

    Args:
        genome_ids: A non-empty collection of genome IDs.
    """
    if not genome_ids:
        raise ValueError("At least one sequence ID must be provided")
    if len(genome_ids) != len(set(genome_ids)):
        raise ValueError("Sequence IDs must be unique")

    return hashlib.sha256(("\n").join(sorted(genome_ids)).encode()).hexdigest()[:8]
