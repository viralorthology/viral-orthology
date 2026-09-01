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
        dir_path: directory containing the FASTA files
        fasta_type: type of the FASTA files
        file_extension: file extension to filter by
    """
    assert isinstance(dir_path, Path)
    assert dir_path.is_dir()
    assert file_extension.startswith(".") and len(file_extension) > 1

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
        FileNotFoundError: If any of the files does not exist.
    """
    assert file_paths

    for file_path in file_paths:
        if not file_path.is_file():
            raise FileNotFoundError(f"{file_path} does not exist")


def ensure_files_do_not_exist(*file_paths: Path) -> None:
    """
    Ensure that none of the given files exist.

    Args:
        *file_paths: Paths to files that must not exist.

    Raises:
        FileExistsError: If any of the files already exists.
    """
    assert file_paths

    for file_path in file_paths:
        if file_path.is_file():
            raise FileExistsError(f"{file_path} already exists")


def ensure_files_have_content(*file_paths: Path) -> None:
    """
    Ensure that all given files exist and have content.

    Raises:
        FileNotFoundError: if any of the files does not exist
        ValueError: if any of the files is empty
    """
    ensure_files_exist(*file_paths)
    for file_path in file_paths:
        if not file_path.stat().st_size > 0:
            raise ValueError(f"{file_path} is empty")


def get_seqs_from_fasta_str(
    fasta_str: str,
) -> list[Seq]:
    """Parse a FASTA string into a list of Seq objects."""
    assert fasta_str

    return [
        Seq(
            seq=record.seq,
            seq_id=record.id,
            seq_description=record.description,
        )
        for record in SeqIO.parse(StringIO(fasta_str), "fasta")
    ]
