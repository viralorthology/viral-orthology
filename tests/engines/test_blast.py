from unittest.mock import patch

import pytest

from engines.blast import BlastHit, blastp_search, make_blast_db, run_blastp
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
def test_make_blastp_db(run_cmd, fasta):
    fasta = Fasta(fasta)
    make_blast_db(fasta, "prot")

    run_cmd.assert_called_once_with(f"makeblastdb -dbtype prot -in {fasta.path}")


@patch("engines.blast.utils.run_cmd")
def test_run_blastp_returns_sorted_hits(run_cmd, tmp_path):
    query = Fasta(tmp_path / "query.fasta")
    db = Fasta(tmp_path / "db.fasta")

    run_cmd.return_value = (
        "query2@subject2@80@90.0@1e-5\nquery1@subject1@95@95.0@1e-20\n"
    )

    hits = run_blastp(query, db, "-evalue 1e-5")

    assert [hit.query_id for hit in hits] == ["query1", "query2"]


@patch("engines.blast.utils.run_cmd")
def test_run_blastp_returns_empty_list(run_cmd, tmp_path):
    query = Fasta(tmp_path / "query.fasta")
    db = Fasta(tmp_path / "db.fasta")

    run_cmd.return_value = " "

    assert run_blastp(query, db, "") == []


def test_search_requires_subject_fasta(tmp_path):
    query = Fasta(tmp_path / "query.fasta")

    with pytest.raises(ValueError):
        blastp_search(query, "")
