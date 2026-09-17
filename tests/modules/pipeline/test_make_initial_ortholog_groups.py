from unittest.mock import Mock, patch

import pytest

from modules.pipeline.make_initial_ortholog_groups import (
    _evaluate_paralogs,
    _get_genome_ids_with_paralogs,
    _remove_paralogs,
)


def make_seq(seq_id, genome_id):
    seq = Mock()
    seq.id = seq_id
    seq.genome_id = genome_id
    return seq


def test_remove_paralogs_all_seqs_are_paralogs(tmp_path):
    og = Mock()
    og.genome_ids = ["A", "A"]

    deleted = _remove_paralogs(og, tmp_path)

    assert deleted is True
    og.delete_fasta.assert_called_once()


@patch("modules.pipeline.make_initial_ortholog_groups.Fasta")
@patch("modules.pipeline.make_initial_ortholog_groups._evaluate_paralogs")
@patch("modules.pipeline.make_initial_ortholog_groups._filter_seqs_for_evaluation")
@patch("modules.pipeline.make_initial_ortholog_groups._get_genome_ids_with_paralogs")
def test_remove_paralogs(
    mock_get_genome_ids_with_paralogs,
    mock_filter_seqs_for_evaluation,
    mock_evaluate_paralogs,
    mock_fasta,
    tmp_path,
):
    seq1 = make_seq("seq1", "genome1")
    seq2 = make_seq("seq2", "genome1")

    og = Mock()
    og.genome_ids = ["A", "B", "A"]
    mock_get_genome_ids_with_paralogs.return_value = ["A"]
    mock_filter_seqs_for_evaluation.return_value = [1, 2]
    mock_evaluate_paralogs.return_value = (seq1, [seq2])

    deleted = _remove_paralogs(og, tmp_path)

    assert deleted is False
    og.delete_fasta.assert_not_called()
    og.remove_seqs.assert_called_once_with(seq2.id)
    mock_fasta.assert_called_once_with(tmp_path / f"{seq1.id}.fasta")
    mock_fasta.return_value.add_seqs.assert_called_once_with(seq2)


@pytest.mark.parametrize(
    ("genome_ids", "expected"),
    [
        (["A", "B", "A"], ["A"]),
        (["C", "A", "B", "A", "C", "C", "B", "D"], ["A", "B", "C"]),
        (["C", "A", "B"], []),
        (["A", "B", "A", "A", "B"], ["A", "B"]),
    ],
)
def test_get_genome_ids_with_paralogs(genome_ids, expected):
    assert _get_genome_ids_with_paralogs(genome_ids) == expected


@pytest.mark.parametrize(
    "genome_ids",
    [
        [],
        ["A"],
        ["A", "A"],
    ],
)
def test_get_genome_ids_with_invalid_input(genome_ids):
    with pytest.raises(AssertionError):
        _get_genome_ids_with_paralogs(genome_ids)


@patch("modules.pipeline.make_initial_ortholog_groups.engines.blastp_search")
@patch("modules.pipeline.make_initial_ortholog_groups.engines.make_blast_db")
@patch("modules.pipeline.make_initial_ortholog_groups.Fasta")
def test_evaluate_paralogs_first_seq_is_the_best(
    _mock_fasta, _mock_make_blast_db, mock_blastp_search
):
    seq1 = make_seq("seq1", "genome1")
    seq2 = make_seq("seq2", "genome1")
    seq3 = make_seq("seq3", "genome2")
    mock_blastp_search.return_value = [Mock(query_id="seq1"), Mock(query_id="seq2")]

    best_seq, other_seqs = _evaluate_paralogs([seq1, seq2], [seq3])

    assert best_seq.id == "seq1"
    assert other_seqs == [seq2]


@patch("modules.pipeline.make_initial_ortholog_groups.engines.blastp_search")
@patch("modules.pipeline.make_initial_ortholog_groups.engines.make_blast_db")
@patch("modules.pipeline.make_initial_ortholog_groups.Fasta")
def test_evaluate_paralogs_second_seq_is_the_best(
    _mock_fasta, _mock_make_blast_db, mock_blastp_search
):
    seq1 = make_seq("seq1", "genome1")
    seq2 = make_seq("seq2", "genome1")
    seq3 = make_seq("seq3", "genome2")
    mock_blastp_search.return_value = [Mock(query_id="seq2"), Mock(query_id="seq1")]

    best_seq, other_seqs = _evaluate_paralogs([seq1, seq2], [seq3])

    assert best_seq.id == "seq2"
    assert other_seqs == [seq1]
