from io import StringIO
from pathlib import Path

from Bio import SeqIO

from models.fasta import Fasta
from models.fasta_type import FastaType
from models.seq import Seq


def get_fastas(
    dir_path: Path, fasta_type: FastaType, file_extension: str
) -> list[Fasta]:
    """
    Get FASTA files from a directory.

    Args:
        dir_path: Directory containing the FASTA files
        fasta_type: Type of the FASTA files
        file_extension: File extension to filter by (it must start with .)

    Raises:
        ValueError: If the file extension is not valid
    """
    if not file_extension.startswith("."):
        raise ValueError(f"{file_extension} is not a valid file extension")

    fastas = [
        Fasta(p, fasta_type)
        for p in dir_path.iterdir()
        if p.is_file() and p.suffix == file_extension
    ]

    assert fastas

    return fastas


def ensure_files_exist(*file_paths: Path) -> None:
    """
    Ensure that all given files exist.

    Args:
        *file_paths: Paths to files that must exist.

    Raises:
        ValueError: If no file path is provided.
        FileNotFoundError: If any of the files do not exist.
    """
    if not file_paths:
        raise ValueError("At least one file path must be provided.")

    for file_path in file_paths:
        if not file_path.is_file():
            raise FileNotFoundError(f"{file_path} does not exist")


def ensure_files_do_not_exist(*file_paths: Path) -> None:
    """
    Ensure that none of the given files exist.

    Args:
        *file_paths: Paths to files that must not exist.

    Raises:
        ValueError: If no file paths are provided.
        FileExistsError: If any of the files already exist.
    """
    if not file_paths:
        raise ValueError("At least one file path must be provided")

    for file_path in file_paths:
        if file_path.exists():
            raise FileExistsError(f"{file_path} already exists")


def ensure_files_have_content(*file_paths: Path) -> None:
    """
    Ensure that all given files exist and have content.

    Args:
        *file_paths: Paths to files that must exist and contain content.

    Raises:
        FileNotFoundError: If no file path is provided or any of the files do not exist
        ValueError: If any of the files is empty
    """
    ensure_files_exist(*file_paths)
    for file_path in file_paths:
        if not file_path.stat().st_size > 0:
            raise ValueError(f"{file_path} is empty")


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
        Seq(
            seq=record.seq,
            seq_id=record.id,
            seq_description=record.description,
        )
        for record in SeqIO.parse(StringIO(fasta_str), "fasta")
    ]
