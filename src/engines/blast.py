import multiprocessing
import shutil
import tempfile
from pathlib import Path
from types import TracebackType
from typing import Literal

from typing_extensions import Self

import utils
from models.fasta import Fasta


class BlastDB:
    """
    Context manager for creating a temporary BLAST database.

    If a single subject FASTA is provided, it is copied directly to the
    temporary directory. If multiple subject FASTAs are provided, they
    are combined into a single FASTA file. The BLAST database is created
    from the resulting FASTA file and is automatically removed when
    leaving the context manager.

    Args:
        db_type: BLAST database type, e.g. "prot" or "nucl".
        *subject_fastas: One or more FASTA objects to use as BLAST
            database subjects.

    Example:
        with BlastDB("prot", fasta1, fasta2) as blast_db:
            run_blast(query, blast_db)
    """

    def __init__(self, db_type: str, *subject_fastas: Fasta):
        self.db_type = db_type
        self.subject_fastas = subject_fastas
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp_dir.name) / "blast_db.fasta"

    def __enter__(self) -> Self:
        if len(self.subject_fastas) == 1:
            shutil.copy(self.subject_fastas[0].path, self.path)
        else:
            utils.get_combined_fasta(self.path, *self.subject_fastas)

        db_fasta = Fasta(self.path)
        make_blast_db(db_fasta, self.db_type)

        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> Literal[False]:
        self.tmp_dir.cleanup()
        return False


class BlastHit:
    """Represents a single BLAST hit."""

    def __init__(self, results_string: str):
        results_split = results_string.split("@")
        assert len(results_split) == 5
        self.query_id = results_split[0]
        self.subject_id = results_split[1]
        self.qcov = int(results_split[2])
        self.ident = float(results_split[3])
        self.evalue = float(results_split[4])


def make_blast_db(fasta: Fasta, db_type: str) -> None:
    """
    Create a BLAST database from a FASTA file.

    Args:
        fasta: FASTA file containing the sequences to index.
        db_type: Type of BLAST database to create (e.g. "prot" or "nucl").
    """
    utils.run_cmd(f"makeblastdb -dbtype {db_type} -in {fasta.path}")


def blastp_search(query_fasta: Fasta, blast_db: BlastDB, params: str) -> list[BlastHit]:
    """
    Run BLASTP against a protein database and return the hits sorted by
    E-value in ascending order (from best to worst).

    Args:
        query_fasta: FASTA file containing the query protein sequences.
        blast_db: BLAST database to search against.
        params: Additional command-line parameters to pass to BLASTP.

    Returns:
        A list of BLAST hits sorted by E-value. Returns an empty list if
        BLAST produces no output.
    """
    num_threads = max(1, multiprocessing.cpu_count() - 1)

    output = utils.run_cmd(
        f"blastp -query {query_fasta.path} -db {blast_db.path} {params} -max_hsps 1 -num_threads {num_threads} -outfmt '6 delim=@ qseqid sseqid qcovs pident evalue'"
    )

    if not output.strip():
        return []

    return sorted(
        [BlastHit(line) for line in output.split("\n") if line.strip()],
        key=lambda hit: hit.evalue,
    )


def run_blastn(query_fasta: Fasta, db_fasta: Fasta, params: str) -> list[BlastHit]:
    """
    Run BLASTN against a nucleotide database and return the hits sorted by evalue.

    Args:
        query_fasta: FASTA file containing the query nucleotide sequences.
        db_fasta: FASTA file containing the nucleotide BLAST database.
        params: Additional command-line parameters to pass to BLASTN.
    """
    num_threads = max(1, multiprocessing.cpu_count() - 1)

    output = utils.run_cmd(
        f"blastn -query {query_fasta.path} -db {db_fasta.path} {params} -num_threads {num_threads} -outfmt '6 delim=@ qseqid sseqid qcovs pident evalue'"
    )

    if not output.strip():
        return []

    return sorted(
        [BlastHit(line) for line in output.split("\n") if line.strip()],
        key=lambda hit: hit.evalue,
    )
