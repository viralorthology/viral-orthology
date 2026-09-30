from types import SimpleNamespace
from unittest.mock import Mock, patch

from engines.hmmer import HMMHit
from modules.pipeline.hmm_search import (
    _add_proteins_to_ogs,
    search_and_add_prots_with_hmm,
)


@patch("modules.pipeline.hmm_search.engines.hmm_search")
def test_search_and_add_prots_with_hmm_searches_incomplete_ogs(
    mock_hmm_search,
):
    og_complete = Mock()
    og_complete.n_seqs = 3

    og_incomplete = Mock()
    og_incomplete.n_seqs = 2

    ui = Mock()
    ui.progress_bar.return_value = [
        og_complete,
        og_incomplete,
    ]

    unique_prots_fasta = Mock()

    mock_hmm_search.return_value = []

    result = search_and_add_prots_with_hmm(
        ui,
        [og_complete, og_incomplete],
        unique_prots_fasta,
        "hmmsearch params",
        3,
    )

    mock_hmm_search.assert_called_once_with(
        og_incomplete,
        unique_prots_fasta,
        "hmmsearch params",
    )

    assert result is False


@patch("modules.pipeline.hmm_search.engines.hmm_search")
def test_search_and_add_prots_with_hmm_returns_false_when_nothing_added(
    mock_hmm_search,
):
    ui = Mock()
    ui.progress_bar.return_value = []

    unique_prots_fasta = Mock()
    mock_hmm_search.return_value = []

    result = search_and_add_prots_with_hmm(
        ui,
        [],
        unique_prots_fasta,
        "params",
        3,
    )

    assert result is False
    unique_prots_fasta.remove_seqs.assert_not_called()


@patch("modules.pipeline.hmm_search.engines.hmm_search")
def test_search_and_add_prots_with_hmm_removes_added_proteins(
    mock_hmm_search,
):
    hit = HMMHit(seq_id="gene1", genome_id="genome1", evalue=1e-50)

    og = Mock()
    og.n_seqs = 1
    og.genome_ids = set()

    seq = Mock()
    seq.id = "gene1"

    ui = Mock()
    ui.progress_bar.return_value = [og]

    unique_prots_fasta = Mock()
    unique_prots_fasta.get_seqs.return_value = [seq]

    mock_hmm_search.return_value = [hit]

    result = search_and_add_prots_with_hmm(
        ui,
        [og],
        unique_prots_fasta,
        "params",
        3,
    )

    assert result is True
    og.add_seqs.assert_called_once_with(seq)
    unique_prots_fasta.remove_seqs.assert_called_once_with("gene1")


def test_add_proteins_to_ogs_adds_protein_to_og():
    hit = HMMHit(
        seq_id="gene1",
        genome_id="genome1",
        evalue=1e-20,
    )

    og = Mock()
    og.genome_ids = set()

    seq = SimpleNamespace(id="gene1")

    fasta = Mock()
    fasta.get_seqs.return_value = [seq]

    added = _add_proteins_to_ogs(
        {hit: og},
        fasta,
    )

    og.add_seqs.assert_called_once_with(seq)
    assert added == {"gene1"}


def test_add_proteins_to_ogs_skips_protein_from_same_genome():
    hit = HMMHit(
        seq_id="gene1",
        genome_id="genome1",
        evalue=1e-20,
    )

    og = Mock()
    og.genome_ids = {"genome1"}

    unique_prots_fasta = Mock()

    added = _add_proteins_to_ogs(
        {hit: og},
        unique_prots_fasta,
    )

    assert added == set()
    og.add_seqs.assert_not_called()
    unique_prots_fasta.get_seqs.assert_not_called()


def test_add_proteins_to_ogs_uses_best_hit():
    hit1 = HMMHit(
        seq_id="gene1",
        genome_id="genome1",
        evalue=1e-10,
    )

    hit2 = HMMHit(
        seq_id="gene1",
        genome_id="genome1",
        evalue=1e-50,
    )

    og1 = Mock()
    og1.genome_ids = set()

    og2 = Mock()
    og2.genome_ids = set()

    seq = SimpleNamespace(id="gene1")

    unique_prots_fasta = Mock()
    unique_prots_fasta.get_seqs.return_value = [seq]

    sorted_hits = {
        hit2: og2,
        hit1: og1,
    }

    added = _add_proteins_to_ogs(
        sorted_hits,
        unique_prots_fasta,
    )

    og2.add_seqs.assert_called_once_with(seq)
    og1.add_seqs.assert_not_called()

    assert added == {"gene1"}


def test_add_proteins_to_ogs_does_not_add_same_protein_twice():
    hit1 = HMMHit(
        seq_id="gene1",
        genome_id="genome1",
        evalue=1e-50,
    )

    hit2 = HMMHit(
        seq_id="gene1",
        genome_id="genome1",
        evalue=1e-20,
    )

    og1 = Mock()
    og1.genome_ids = set()

    og2 = Mock()
    og2.genome_ids = set()

    seq = SimpleNamespace(id="gene1")

    unique_prots_fasta = Mock()
    unique_prots_fasta.get_seqs.return_value = [seq]

    added = _add_proteins_to_ogs(
        {hit1: og1, hit2: og2},
        unique_prots_fasta,
    )

    assert added == {"gene1"}
    og1.add_seqs.assert_called_once_with(seq)

    assert og1.add_seqs.call_count + og2.add_seqs.call_count == 1


def test_add_proteins_to_ogs_adds_multiple_proteins():
    hit1 = HMMHit(
        seq_id="gene1",
        genome_id="genome1",
        evalue=1e-50,
    )

    hit2 = HMMHit(
        seq_id="gene2",
        genome_id="genome2",
        evalue=1e-40,
    )

    og1 = Mock()
    og1.genome_ids = set()

    og2 = Mock()
    og2.genome_ids = set()

    seq1 = SimpleNamespace(id="gene1")
    seq2 = SimpleNamespace(id="gene2")

    unique_prots_fasta = Mock()
    unique_prots_fasta.get_seqs.side_effect = [[seq1], [seq2]]

    added = _add_proteins_to_ogs(
        {hit1: og1, hit2: og2},
        unique_prots_fasta,
    )

    og1.add_seqs.assert_called_once_with(seq1)
    og2.add_seqs.assert_called_once_with(seq2)

    assert added == {"gene1", "gene2"}
