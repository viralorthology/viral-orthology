import shutil
from collections.abc import Iterator
from pathlib import Path

from Bio import SeqIO

from models.fasta_type import FastaType
from models.seq import Seq, get_seq_from_seqrecord, get_seqrecord_from_seq


class Fasta:
    """
    Represents a FASTA file and provides common operations.

    An instance of Fasta can be created without the file existing.

    Attributes:
        path: Path of the file
        fasta_type: type of sequences in the file
    """

    def __init__(self, path: Path, fasta_type: FastaType):
        self.path = path
        self.fasta_type = fasta_type

    @property
    def seqs(self) -> Iterator[Seq]:
        """
        Return an iterator over the sequences in the FASTA file.

        Raises:
            FileNotFoundError: if the file does not exist
            ValueError: if the file is empty, or if it contains a malformed record or duplicate ID
        """
        if not self.path.is_file():
            raise FileNotFoundError(f"{self.path} does not exist")
        if self.path.stat().st_size == 0:
            raise ValueError(f"{self.path} is empty")

        ids = set()
        for seq in SeqIO.parse(self.path, "fasta"):
            if seq.id in ids:
                raise ValueError(
                    f"{self.path} contains duplicate sequence ID: {seq.id}"
                )

            ids.add(seq.id)
            yield get_seq_from_seqrecord(seq, self.fasta_type, self.path)

    @property
    def n_seqs(self) -> int:
        return sum(1 for _ in self.seqs)

    @property
    def ids(self) -> list[str]:
        """
        Return the IDs of all sequences in the FASTA file.
        """
        return [seq.id for seq in self.seqs]

    @property
    def genome_ids(self) -> list[str]:
        """
        Return the genome IDs of all protein sequences.

        Raises:
            ValueError: if used on a non-protein FASTA
        """
        if self.fasta_type != FastaType.PROTEIN:
            raise ValueError(
                f"genome_ids is only available for protein FASTAs: {self.path}"
            )

        return [seq.genome_id for seq in self.seqs]

    def get_seqs(self, *ids: str) -> list[Seq]:
        """
        Return the sequences in the same order as the provided IDs.

        Raises:
            ValueError: If no IDs are provided, if the IDs are not unique, or if
                any of the requested sequences is not found in the FASTA file.
        """
        if not ids:
            raise ValueError("At least one sequence ID must be provided")
        if len(ids) != len(set(ids)):
            raise ValueError("Sequence IDs must be unique")

        seqs = [seq for seq in self.seqs if seq.id in ids]

        if len(seqs) != len(ids):
            raise ValueError(f"Some sequence was not found in {self.path}")

        positions = {id_: i for i, id_ in enumerate(ids)}
        return sorted(seqs, key=lambda seq: positions[seq.id])

    def add_seqs(self, *seqs: Seq) -> None:
        """
        Append the given sequences to the FASTA file.

        Raises:
            ValueError: If no sequences are provided.
        """
        if not seqs:
            raise ValueError("At least one sequence must be provided")

        seq_records = [get_seqrecord_from_seq(seq) for seq in seqs]

        with self.path.open("a", encoding="utf-8") as fh:
            SeqIO.write(seq_records, fh, "fasta")

    def remove_seqs(self, *ids: str) -> None:
        """
        Remove the sequences with the given IDs from the FASTA file.

        Raises:
            ValueError: If no IDs are provided, if the IDs are not unique, or if
                any of the requested sequences is not found in the FASTA file.
        """
        if not ids:
            raise ValueError("At least one sequence ID must be provided")
        if len(ids) != len(set(ids)):
            raise ValueError("Sequence IDs must be unique")

        ids_to_remove = set(ids)
        missing_ids = ids_to_remove - set(self.ids)
        if missing_ids:
            raise ValueError(f"Some IDs were not found in {self.path}: {missing_ids}")
        seqs_to_keep = [seq for seq in self.seqs if seq.id not in ids_to_remove]

        self.delete_fasta()
        if seqs_to_keep:
            self.add_seqs(*seqs_to_keep)

    def rename_fasta(self, new_filename: str) -> None:
        """
        Rename the FASTA file.

        Raises:
            FileNotFoundError: If the FASTA file does not exist.
            ValueError: If no filename is provided, if the FASTA file is empty,
                or if the new filename has no suffix.
            FileExistsError: If a file with the new name already exists.
        """
        if not new_filename:
            raise ValueError("A new filename must be provided")
        if not self.path.is_file():
            raise FileNotFoundError(f"{self.path} does not exist")
        if self.path.stat().st_size == 0:
            raise ValueError(f"{self.path} is empty")

        new_path = self.path.with_name(new_filename)
        if not new_path.suffix:
            raise ValueError("new filename has to have a suffix")
        if new_path.exists():
            raise FileExistsError(f"{new_path} already exists")

        self.path = self.path.rename(new_path)

    def move_fasta(self, directory_path: Path) -> None:
        """
        Move the file from its current directory to the given directory.

        Raises:
            FileNotFoundError: if the file does not exist or the given directory does not exist
            ValueError: if the file is empty
            FileExistsError: if the destination path already exists
        """
        if not self.path.is_file():
            raise FileNotFoundError(f"{self.path} does not exist")
        if self.path.stat().st_size == 0:
            raise ValueError(f"{self.path} is empty")

        new_path = directory_path / self.path.name
        if new_path.exists():
            raise FileExistsError(f"{new_path} already exists")

        shutil.move(self.path, new_path)
        self.path = new_path

    def delete_fasta(self) -> None:
        """
        Delete the file.

        Raises:
            FileNotFoundError: if the file does not exist
        """
        self.path.unlink()
