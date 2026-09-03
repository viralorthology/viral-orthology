from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from engines.hmmer import (
    HMMHit,
    _parse_hmmsearch_output,
    _search_hmm_hmmer,
    hmm_search,
)


@pytest.fixture
def query_fasta():
    query = MagicMock()
    query.alignment_path = Path("/tmp/query.alignment.fasta")
    query.hmm_hmmer_path = Path("/tmp/query.hmm")
    return query


@pytest.fixture
def db_fasta():
    db = MagicMock()
    db.path = Path("/tmp/database.fasta")
    return db


def test_search_hmm_hmmer(query_fasta, db_fasta):
    with (
        patch(
            "engines.hmmer.utils.run_cmd",
        ),
        patch(
            "engines.hmmer._parse_hmmsearch_output",
            return_value=[
                HMMHit(
                    evalue=1e-20,
                    seq_id="protein1",
                    genome_id="genome1",
                )
            ],
        ),
    ):
        hits = _search_hmm_hmmer(query_fasta, db_fasta, "")

    assert len(hits) == 1
    assert hits[0].evalue == 1e-20
    assert hits[0].seq_id == "protein1"
    assert hits[0].genome_id == "genome1"


def test_parse_hmmsearch_output():
    output = """
    Some HMMER header information

    E-value  score  bias
    1e-10    100.0  0.1  0  0  0  0 0 proteinA genomeA
    1e-50    150.0  0.2  0  0  0  0 0 proteinB genomeB

    Domain annotation for each sequence
    """

    hits = _parse_hmmsearch_output(output)

    assert hits == [
        HMMHit(
            evalue=1e-50,
            seq_id="proteinB",
            genome_id="genomeB",
        ),
        HMMHit(
            evalue=1e-10,
            seq_id="proteinA",
            genome_id="genomeA",
        ),
    ]


def test_parse_hmmsearch_output_sorts_by_evalue():
    output = """
    random stuff
    E-value  score  bias
    1e-5     100.0  0.1  0  0  0  0 0 proteinA genomeA
    1e-20    150.0  0.2  0  0  0  0 0 proteinB genomeB
    1e-10    120.0  0.3  0  0  0  0 0 proteinC genomeC

    Domain annotation for each sequence
    """

    hits = _parse_hmmsearch_output(output)

    assert [hit.evalue for hit in hits] == [
        1e-20,
        1e-10,
        1e-5,
    ]


def test_parse_hmmsearch_output_stops_at_inclusion_threshold():
    output = """
    E-value  score  bias
    1e-10    100.0  0.1  0  0  0  0 0 proteinA genomeA

    ------ inclusion threshold ------
    1e-50    200.0  0.1  0  0  0  0  proteinB genomeB
    """

    hits = _parse_hmmsearch_output(output)

    assert hits == [
        HMMHit(
            evalue=1e-10,
            seq_id="proteinA",
            genome_id="genomeA",
        )
    ]


def test_parse_hmmsearch_output_returns_empty_list():
    output = """
    Some HMMER output

    E-value  score  bias

    Domain annotation for each sequence
    """

    hits = _parse_hmmsearch_output(output)

    assert hits == []


def test_hmm_search_returns_search_results(
    query_fasta,
    db_fasta,
):
    expected_hits = [
        HMMHit(
            evalue=1e-30,
            seq_id="proteinA",
            genome_id="genomeA",
        ),
        HMMHit(
            evalue=1e-10,
            seq_id="proteinB",
            genome_id="genomeB",
        ),
    ]

    with (
        patch.object(Path, "is_file", return_value=True),
        patch(
            "engines.hmmer._search_hmm_hmmer",
            return_value=expected_hits,
        ) as search_hmm,
    ):
        hits = hmm_search(
            query_fasta,
            db_fasta,
            "-test 1e-5",
        )

    assert hits == expected_hits

    search_hmm.assert_called_once_with(
        query_fasta,
        db_fasta,
        "-test 1e-5",
    )
