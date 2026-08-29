from collections.abc import Iterator
from pathlib import Path

from Bio import SeqIO

from models.fasta_type import FastaType
from models.file import File
from models.seq import Seq


class Fasta(File):
    """
    Represents a FASTA file and provides common operations.

    An instance of Fasta can be created without the file existing.

    Attributes:
        path: Path of the file
        fasta_type: type of sequences in the file
    """

    def __init__(self, path: Path, fasta_type: FastaType):
        super().__init__(path)
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
        if not self.has_content:
            raise ValueError(f"{self.path} is empty")

        ids = set()
        for seq in SeqIO.parse(self.path, "fasta"):
            if not seq.id:
                raise ValueError(f"{self.path} contains a malformed FASTA record")
            if seq.id in ids:
                raise ValueError(
                    f"{self.path} contains duplicate sequence ID: {seq.id}"
                )

            ids.add(seq.id)
            yield Seq(seq, self.fasta_type)

    @property
    def n_seqs(self) -> int:
        return len(list(self.seqs))

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
        """
        assert self.fasta_type == FastaType.PROTEIN
        return [seq.genome_id for seq in self.seqs]

    def get_seqs(self, *ids: str) -> list[Seq]:
        """
        Return the sequences with the given IDs.

        Raises:
            ValueError: if any of the requested sequences was not found
        """
        assert ids
        assert len(ids) == len(set(ids))

        seqs = [seq for seq in self.seqs if seq.id in ids]

        if len(seqs) != len(ids):
            raise ValueError(f"Some sequence was not found in {self.path}")

        positions = {id_: i for i, id_ in enumerate(ids)}
        return sorted(seqs, key=lambda seq: positions[seq.id])

    def add_seqs(self, *seqs: Seq) -> None:
        """
        Append the given sequences to the FASTA file.
        """
        assert seqs

        with self.path.open("a", encoding="utf-8") as fh:
            SeqIO.write(seqs, fh, "fasta")

    def remove_seqs(self, *ids: str) -> None:
        """
        Remove the sequences with the given IDs from the FASTA file.

        Raises:
            ValueError: if any of the sequences was not found
        """
        assert ids
        assert len(ids) == len(set(ids))

        ids_to_remove = set(ids)
        missing_ids = ids_to_remove - set(self.ids)
        if missing_ids:
            raise ValueError(f"Some IDs were not found in {self.path}: {missing_ids}")
        seqs_to_keep = [seq for seq in self.seqs if seq.id not in ids_to_remove]

        self.delete_file()
        if seqs_to_keep:
            self.add_seqs(*seqs_to_keep)
