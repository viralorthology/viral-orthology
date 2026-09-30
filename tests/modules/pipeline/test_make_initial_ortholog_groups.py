from types import SimpleNamespace
from unittest.mock import Mock, call, patch

import pytest

from modules.pipeline.make_initial_ortholog_groups import (
    _get_genome_ids_with_paralogs,
    _remove_paralogs,
)


def test_remove_paralogs_all_seqs_are_paralogs(tmp_path):
    og = Mock()
    og.genome_ids = ["A", "A"]

    deleted = _remove_paralogs(og, tmp_path)

    assert deleted is True
    og.delete_fasta.assert_called_once()


def test_remove_paralogs_keeps_longest_sequence():
    seq1 = SimpleNamespace(
        id="gene1",
        genome_id="genome1",
        seq="A" * 100,
    )
    seq2 = SimpleNamespace(
        id="gene2",
        genome_id="genome1",
        seq="A" * 150,
    )
    seq3 = SimpleNamespace(
        id="gene3",
        genome_id="genome2",
        seq="A" * 120,
    )

    og = Mock()
    og.genome_ids = ["genome1", "genome1", "genome2"]
    og.seqs = [seq1, seq2, seq3]
    og.path.exists.return_value = True

    paralog_ids = set()

    with patch(
        "modules.pipeline.make_initial_ortholog_groups._get_genome_ids_with_paralogs",
        return_value=["genome1"],
    ):
        result = _remove_paralogs(og, paralog_ids)

    assert result is False
    og.remove_seqs.assert_called_once_with("gene1")
    assert paralog_ids == {"gene1"}


def test_remove_paralogs_keeps_longest_sequence_per_genome():
    seqs = [
        SimpleNamespace(id="gene1", genome_id="genome1", seq="A" * 100),
        SimpleNamespace(id="gene2", genome_id="genome1", seq="A" * 150),
        SimpleNamespace(id="gene3", genome_id="genome2", seq="A" * 200),
        SimpleNamespace(id="gene4", genome_id="genome2", seq="A" * 120),
    ]

    og = Mock()
    og.genome_ids = [
        "genome1",
        "genome1",
        "genome2",
        "genome2",
    ]
    og.seqs = seqs
    og.path.exists.return_value = True

    paralog_ids = set()

    with patch(
        "modules.pipeline.make_initial_ortholog_groups._get_genome_ids_with_paralogs",
        return_value=["genome1", "genome2"],
    ):
        result = _remove_paralogs(og, paralog_ids)

    assert result is False

    assert og.remove_seqs.call_args_list == [
        call("gene1"),
        call("gene4"),
    ]

    assert paralog_ids == {"gene1", "gene4"}


def test_remove_paralogs_does_nothing_when_there_are_no_paralogs():
    seq1 = SimpleNamespace(id="gene1", genome_id="genome1")
    seq2 = SimpleNamespace(id="gene2", genome_id="genome2")

    og = Mock()
    og.genome_ids = ["genome1", "genome2"]
    og.seqs = [seq1, seq2]

    paralog_ids = set()

    with patch(
        "modules.pipeline.make_initial_ortholog_groups._get_genome_ids_with_paralogs",
        return_value=[],
    ):
        result = _remove_paralogs(og, paralog_ids)

    assert result is False
    og.remove_seqs.assert_not_called()
    og.delete_fasta.assert_not_called()
    assert paralog_ids == set()


def test_remove_paralogs_uses_sequence_id_as_tiebreaker():
    seq1 = SimpleNamespace(
        id="gene1",
        genome_id="genome1",
        seq="A" * 150,
    )
    seq2 = SimpleNamespace(
        id="gene2",
        genome_id="genome1",
        seq="A" * 150,
    )
    seq3 = SimpleNamespace(
        id="gene3",
        genome_id="genome2",
        seq="A" * 100,
    )

    og = Mock()
    og.genome_ids = ["genome1", "genome1", "genome2"]
    og.seqs = [seq1, seq2, seq3]
    og.path.exists.return_value = True

    paralog_ids = set()

    with patch(
        "modules.pipeline.make_initial_ortholog_groups._get_genome_ids_with_paralogs",
        return_value=["genome1"],
    ):
        result = _remove_paralogs(og, paralog_ids)

    assert result is False
    og.remove_seqs.assert_called_once_with("gene1")
    assert paralog_ids == {"gene1"}


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
