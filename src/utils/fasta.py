from io import StringIO
from pathlib import Path

from Bio import SeqIO

from models.fasta import Fasta
from models.ortholog_group import OrthologGroup
from models.seq import Seq, get_seq_from_seqrecord


def get_fastas(dir_path: Path) -> list[Fasta]:
    """Get FASTA files from a directory."""
    fastas = [
        Fasta(p) for p in dir_path.iterdir() if p.is_file() and p.suffix == ".fasta"
    ]

    if not fastas:
        raise FileNotFoundError(f"No fastas found in {dir_path}")

    return sorted(
        fastas, key=lambda fasta: (-fasta.n_seqs, fasta.path.stem)
    )  # stem provides a deterministic tie-breaker.


def get_ortholog_groups(dir_path: Path) -> list[OrthologGroup]:
    """
    Get all OrthologGroup FASTA files from a directory.

    Returns:
        Ortholog groups sorted by number of sequences in descending order,
        with filename used as the secondary sorting criterion.
    """
    ogs = [
        OrthologGroup(p)
        for p in dir_path.iterdir()
        if p.is_file() and p.suffix == ".fasta"
    ]

    if not ogs:
        raise FileNotFoundError(f"No ortholog groups found in {dir_path}")

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
