import multiprocessing
import tempfile
from pathlib import Path

import utils
from models.fasta import Fasta


class BlastHit:
    """Represents a single BLASTP hit."""

    def __init__(self, results_string: str):
        results_split = results_string.split("@")
        assert len(results_split) == 5
        self.query_id = results_split[0]
        self.subject_id = results_split[1]
        self.qcov = int(results_split[2])
        self.ident = float(results_split[3])
        self.evalue = float(results_split[4])


def search(query_fasta: Fasta, params: str, *subject_fastas: Fasta) -> list[BlastHit]:
    """
    Search query proteins against one or more subject FASTA files and
    get the hits sorted by evalue.

    The subject FASTA files are combined into a temporary protein BLAST
    database. The database and all associated files are automatically
    removed after the search completes.

    Args:
        query_fasta: FASTA file containing the query protein sequences.
        params: Additional command-line parameters to pass to BLASTP.
        *subject_fastas: FASTA files containing the subject protein sequences.

    Raises:
        ValueError: If no subject FASTA files are provided.
    """

    if not subject_fastas:
        raise ValueError("At least one subject fasta must be provided")

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir_path = Path(tmp_dir)

        blastp_db_fasta = utils.get_combined_fasta(
            tmp_dir_path / "blastp_db.fasta", *subject_fastas
        )
        make_blastp_db(blastp_db_fasta)

        blastp_hits = run_blastp(
            query_fasta, blastp_db_fasta, f"{params} -max_target_seqs 1"
        )  # only best hit per query prot

        hits_per_query = [hit.query_id for hit in blastp_hits]
        assert len(hits_per_query) == len(set(hits_per_query))

        return blastp_hits


def make_blastp_db(fasta: Fasta) -> None:
    """
    Create a protein BLASTP database from a FASTA file.

    Args:
        fasta: FASTA file containing the protein sequences to index.
    """
    utils.run_cmd(f"makeblastdb -dbtype prot -in {fasta.path}")


def run_blastp(query_fasta: Fasta, db_fasta: Fasta, params: str) -> list[BlastHit]:
    """
    Run BLASTP against a protein database and return the hits sorted by evalue.

    Args:
        query_fasta: FASTA file containing the query protein sequences.
        db_fasta: FASTA file containing the protein BLAST database.
        params: Additional command-line parameters to pass to BLASTP.
    """
    num_threads = max(1, multiprocessing.cpu_count() - 1)

    output = utils.run_cmd(
        f"blastp -query {query_fasta.path} -db {db_fasta.path} {params} -max_hsps 1 -num_threads {num_threads} -outfmt '6 delim=@ qseqid sseqid qcovs pident evalue'"
    )

    if not output.strip():
        return []

    return sorted(
        [BlastHit(line) for line in output.split("\n") if line.strip()],
        key=lambda hit: hit.evalue,
    )
