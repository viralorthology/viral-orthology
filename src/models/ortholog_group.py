from pathlib import Path

from models.fasta import Fasta
from models.fasta_type import FastaType
from models.seq import Seq


class OrthologGroup(Fasta):
    """
    Represents a FASTA file containing an ortholog group.

    An instance can only be created if the FASTA file contains a valid ortholog group. See _validate() for the validation rules.
    """

    def __init__(
        self, path: Path, annotated_prot_db_path: Path, orffinder_prot_db_path: Path
    ):
        super().__init__(path, FastaType.PROTEIN)
        self.annotated_prot_db_path = annotated_prot_db_path
        self.orffinder_prot_db_path = orffinder_prot_db_path
        self._validate()

    @property
    def alignment_path(self) -> Path:
        return self.path.with_suffix(".muscle")

    @property
    def hmm_hmmer_path(self) -> Path:
        return self.path.with_suffix(".hmm")

    @property
    def hmm_hhsuite_path(self) -> Path:
        return self.path.with_suffix(".hhm")

    @property
    def a2m_path(self) -> Path:
        return self.path.with_suffix(".a2m")

    @property
    def associated_files(self) -> list[Path]:
        return [
            self.alignment_path,
            self.hmm_hmmer_path,
            self.hmm_hhsuite_path,
            self.a2m_path,
        ]

    def add_seqs(self, *seqs: Seq) -> None:
        """Add sequences to the ortholog group and revalidate it."""
        super().add_seqs(*seqs)
        self._delete_associated_files()
        self._validate()

    def remove_seqs(self, *ids: str) -> None:
        """Remove sequences from the ortholog group.

        If only one sequence remains, it is moved to the corresponding protein db and the ortholog group FASTA file is deleted.
        """
        self._delete_associated_files()
        remaining_seq_ids = set(self.ids) - set(ids)

        if len(remaining_seq_ids) == 1:
            remaining_seq = self.get_seqs(next(iter(remaining_seq_ids)))[0]
            prot_db_path = (
                self.orffinder_prot_db_path
                if remaining_seq.id.startswith("ORFFINDER")
                else self.annotated_prot_db_path
            )
            prot_db = Fasta(prot_db_path, FastaType.GENERIC)
            prot_db.add_seqs(remaining_seq)
            super().remove_seqs(*ids, remaining_seq.id)  # check ids
            assert not self.exists
        else:
            super().remove_seqs(*ids)
            if self.exists:
                self._validate()

    def delete_file(self) -> None:
        """Delete the ortholog group FASTA file and its associated files."""
        self._delete_associated_files()
        super().delete_file()

    def rename_file(self, new_filename: str) -> None:
        """Rename the ortholog group FASTA file and delete its associated files."""
        self._delete_associated_files()  # associated files depend upon self.path, so delete them before renaming the file
        super().rename_file(new_filename)

    def _validate(self) -> None:
        """
        Check if the ortholog group is valid.

        It has to have at least two sequences, and no two sequences from the same genome.
        """
        genome_ids = self.genome_ids
        n_seqs = len(genome_ids)

        if n_seqs < 2:
            raise ValueError(f"{self.path} ortholog group has fewer than two sequences")

        if len(genome_ids) != len(set(genome_ids)):
            raise ValueError(
                f"{self.path} ortholog group contains multiple sequences from the same genome"
            )

    def _delete_associated_files(self) -> None:
        for path in self.associated_files:
            if path.is_file():
                path.unlink()
