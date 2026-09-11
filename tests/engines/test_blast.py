from unittest.mock import MagicMock, patch

import pytest

from engines.blast import BlastDB, BlastHit, blastp_search, make_blast_db
from models.fasta import Fasta


def test_blast_result_parses_result():
    result = BlastHit("query1@subject1@95@87.5@1e-20")

    assert result.query_id == "query1"
    assert result.subject_id == "subject1"
    assert result.qcov == 95
    assert result.ident == 87.5
    assert result.evalue == 1e-20


def test_blast_result_rejects_invalid_result():
    with pytest.raises(AssertionError):
        BlastHit("query1@subject1@95")


@patch("engines.blast.utils.run_cmd")
def test_make_blastp_db(run_cmd, tmp_path):
    fasta_path = tmp_path / "fasta.fasta"
    fasta_path.write_text(">seq1 description\nMKMMKMKMKMMMKMKMK\n")
    fasta = Fasta(fasta_path)
    make_blast_db(fasta, "prot")

    run_cmd.assert_called_once_with(f"makeblastdb -dbtype prot -in {fasta.path}")


@patch("engines.blast.utils.run_cmd")
def test_run_blastp_returns_sorted_hits(run_cmd, tmp_path):
    query = Fasta(tmp_path / "query.fasta")
    blast_db = MagicMock()
    blast_db.path = tmp_path / "blast_db"

    run_cmd.return_value = (
        "query2@subject2@80@90.0@1e-5\nquery1@subject1@95@95.0@1e-20\n"
    )

    hits = blastp_search(query, blast_db, "-evalue 1e-5")

    assert [hit.query_id for hit in hits] == ["query1", "query2"]


@patch("engines.blast.utils.run_cmd")
def test_run_blastp_returns_empty_list(run_cmd, tmp_path):
    query = Fasta(tmp_path / "query.fasta")
    blast_db = MagicMock()
    blast_db.path = tmp_path / "blast_db"

    run_cmd.return_value = " "

    assert blastp_search(query, blast_db, "") == []


@patch("engines.blast.make_blast_db")
def test_blast_db_creates_database(make_blast_db, tmp_path):
    fasta = Fasta(tmp_path / "db.fasta")
    fasta.path.touch()

    with BlastDB("prot", fasta) as blast_db:
        assert blast_db.path.exists()
        make_blast_db.assert_called_once()


@patch("engines.blast.make_blast_db")
def test_blast_db_cleans_up_after_context(make_blast_db, tmp_path):
    fasta = Fasta(tmp_path / "db.fasta")
    fasta.path.touch()

    with BlastDB("prot", fasta) as blast_db:
        db_path = blast_db.path
        assert db_path.exists()

    assert not db_path.exists()
