from collections.abc import Iterator
from pathlib import Path

from config.context import Context
from models.fasta import Fasta
from models.seq import Seq


class OrthologGroup:
    """
    Represents a FASTA file containing an ortholog group.

    An instance can only be created if the FASTA file contains a valid ortholog group. See _validate() for the validation rules.
    """

    def __init__(self, path: Path):
        self._fasta = Fasta(path)
        self._validate_ortholog_group()

    @property
    def path(self) -> Path:
        return self._fasta.path

    @property
    def ids(self) -> list[str]:
        return self._fasta.ids

    @property
    def genome_ids(self) -> list[str]:
        return self._fasta.genome_ids

    @property
    def seqs(self) -> Iterator[Seq]:
        return self._fasta.seqs

    @property
    def n_seqs(self) -> int:
        return self._fasta.n_seqs

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

    def get_seqs(self, *ids: str) -> list[Seq]:
        return self._fasta.get_seqs(*ids)

    def add_seqs(self, *seqs: Seq) -> None:
        """Add sequences to the ortholog group and revalidate it."""
        self._fasta.add_seqs(*seqs)
        self._delete_associated_files()
        self._validate_ortholog_group()

    def remove_seqs(self, *ids: str, ctx: Context) -> None:
        """
        Remove sequences from the ortholog group.

        If only one sequence remains, it is moved to the corresponding protein db and the ortholog group FASTA file is deleted.

        Raises:
            ValueError: If ctx is not provided.
        """
        self._delete_associated_files()
        remaining_seq_ids = set(self.ids) - set(ids)

        if len(remaining_seq_ids) == 1:
            remaining_seq = self.get_seqs(next(iter(remaining_seq_ids)))[0]
            prot_db_path = (
                ctx.paths.predicted_unique_prots_fasta
                if remaining_seq.is_predicted
                else ctx.paths.annotated_unique_prots_fasta
            )
            prot_db = Fasta(prot_db_path)
            prot_db.add_seqs(remaining_seq)
            self._fasta.remove_seqs(
                *ids, remaining_seq.id
            )  # Fasta.remove_seqs will check the ids and delete the file
            assert not self.path.is_file()
        else:
            self._fasta.remove_seqs(*ids)
            if self.path.exists():
                self._validate_ortholog_group()

    def delete_fasta(self) -> None:
        """Delete the ortholog group FASTA file and its associated files."""
        self._delete_associated_files()
        self._fasta.delete_fasta()

    def rename_fasta(self, new_filename: str) -> None:
        """Rename the ortholog group FASTA file and delete its associated files."""
        self._delete_associated_files()  # associated files depend upon self.path, so delete them before renaming the file
        self._fasta.rename_fasta(new_filename)

    def _validate_ortholog_group(self) -> None:
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
            if path.exists():
                path.unlink()
