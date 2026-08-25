from pathlib import Path

from classes.fasta import Fasta
from classes.fasta_type import FastaType
from classes.seq import Seq


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
        super().add_seqs(*seqs)
        self._delete_associated_files()
        self._validate()

    def remove_seqs(self, *ids: str) -> None:
        ids_to_remove = set(ids)
        fasta_ids = set(self.ids)
        remaining_seq_ids = fasta_ids - ids_to_remove

        if (
            len(remaining_seq_ids) == 1
        ):  # move remaining seq to the corresponding prot db
            assert len(ids) == len(ids_to_remove)
            missing_ids = ids_to_remove - fasta_ids
            if missing_ids:
                raise ValueError(
                    f"Some IDs were not found in {self.path}: {missing_ids}"
                )

            remaining_seq = self.get_seqs(next(iter(remaining_seq_ids)))[0]
            prot_db_path = (
                self.orffinder_prot_db_path
                if remaining_seq.id.startswith("ORFFINDER")
                else self.annotated_prot_db_path
            )
            prot_db = Fasta(prot_db_path, FastaType.GENERIC)
            prot_db.add_seqs(remaining_seq)
            self.delete_file()
        else:
            super().remove_seqs(*ids)
            self._delete_associated_files()
            if self.exists:
                self._validate()

    def delete_file(self) -> None:
        self._delete_associated_files()
        super().delete_file()

    def rename_file(self, new_filename: str) -> None:
        self._delete_associated_files()  # associated files depend upon self.path, so delete them before renaming the file
        super().rename_file(new_filename)

    def _validate(self) -> None:
        """
        Check if the ortholog group is valid.

        It has to have at least 2 sequences, and no 2 sequences from the same genome.
        """
        genome_ids = self.genome_ids
        n_seqs = len(genome_ids)

        if n_seqs < 2:
            raise ValueError(f"{self.path} ortholog group has less than 2 sequences")

        if len(genome_ids) != len(set(genome_ids)):
            raise ValueError(
                f"{self.path} ortholog group contains multiple sequences from the same genome"
            )

    def _delete_associated_files(self) -> None:
        for path in self.associated_files:
            if path.is_file():
                path.unlink()
