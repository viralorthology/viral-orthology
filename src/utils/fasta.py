from io import StringIO
from pathlib import Path

from Bio import SeqIO

from models.fasta import Fasta
from models.ortholog_group import OrthologGroup
from models.seq import Seq, get_seq_from_seqrecord


def o_groups_are_compatible(*ortholog_groups: OrthologGroup) -> bool:  # TODO delete
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


def get_fastas(dir_path: Path) -> list[Fasta]:
    """
    Get FASTA files from a directory.

    Args:
        dir_path: Directory containing the FASTA files
    """

    fastas = [
        Fasta(p) for p in dir_path.iterdir() if p.is_file() and p.suffix == ".fasta"
    ]

    assert fastas

    return fastas


def get_ortholog_groups(dir_path: Path) -> list[OrthologGroup]:
    """
    Load all OrthologGroup FASTA files from a directory.

    Returns:
        Ortholog groups sorted by number of sequences in descending order,
        with filename used as the secondary sorting criterion.
    """

    ogs = [
        OrthologGroup(p)
        for p in dir_path.iterdir()
        if p.is_file() and p.suffix == ".fasta"
    ]

    assert ogs

    return sorted(
        ogs, key=lambda og: (-og.n_seqs, og.path.stem)
    )  # stem provides a deterministic tie-breaker.


def get_seqs_from_fasta_str(
    fasta_str: str,
) -> list[Seq]:
    """
    Parse a FASTA string into a list of Seq objects.

    Raises:
        ValueError: If fasta_str is empty.
    """
    if not fasta_str:
        raise ValueError("fasta_str cannot be empty")

    return [
        get_seq_from_seqrecord(record)
        for record in SeqIO.parse(StringIO(fasta_str), "fasta")
    ]
