from unittest.mock import call, patch

import pytest
from Bio.Seq import Seq as BioSeq

from models.seq import Seq
from modules.download_seqs.run import (
    _analyze_genomes,
    _download_fasta,
    _get_genome_ids_to_download,
    _get_seq_id_and_description_protein_seqs,
)


def test_get_genome_ids_to_download():
    ids = "ABC123.1\nABC123.1\nABC123.2\nBCD123.1\n\n"
    assert _get_genome_ids_to_download(ids) == {"ABC123", "BCD123"}


def test_analyze_genomes():
    genome_seqs = [
        Seq(BioSeq("GGCCATNN"), "genome_1", ""),
        Seq(BioSeq("GGCCATNN"), "genome_2", ""),
        Seq(BioSeq("ATATATAT"), "genome_3", ""),
    ]

    genome_lens, genome_gc_perc, genome_n_counts, identical_genomes = _analyze_genomes(
        genome_seqs
    )

    assert genome_lens == {
        "genome_1": 8,
        "genome_2": 8,
        "genome_3": 8,
    }
    assert genome_gc_perc == {
        "genome_1": 66.67,
        "genome_2": 66.67,
        "genome_3": 0.0,
    }
    assert genome_n_counts == {
        "genome_1": 2,
        "genome_2": 2,
        "genome_3": 0,
    }
    assert identical_genomes == [
        ("genome_1", "genome_2"),
    ]


def test_analyze_genomes_duplicate_genomes():
    genome_seqs = [
        Seq(BioSeq("GGCCATNN"), "genome1", "genome1"),
        Seq(BioSeq("GGCCATNN"), "genome1", "genome1"),
    ]
    with pytest.raises(AssertionError):
        _ = _analyze_genomes(genome_seqs)


@pytest.mark.parametrize(
    "old_description, expected_seq_id, expected_description",
    [
        (
            "protein|GCF_000001405.40_prot_1 DNA-directed RNA polymerase subunit beta [protein_id=NP_000001.2]",
            "NP_000001.2",
            "GCF_000001405.40 DNA-directed RNA polymerase subunit beta [protein_id=NP_000001.2]",
        ),
        (
            "protein|GCF_000002315.6_prot_42 hypothetical protein [protein_id=WP_012345678.1]",
            "WP_012345678.1",
            "GCF_000002315.6 hypothetical protein [protein_id=WP_012345678.1]",
        ),
        (
            "protein|GCF_000005845.3_prot_108 ATP synthase subunit beta [protein_id=YP_009724390.1]",
            "YP_009724390.1",
            "GCF_000005845.3 ATP synthase subunit beta [protein_id=YP_009724390.1]",
        ),
        (
            "protein|GCF_000006765.12_prot_256 putative membrane protein [protein_id=WP_098765432.2]",
            "WP_098765432.2",
            "GCF_000006765.12 putative membrane protein [protein_id=WP_098765432.2]",
        ),
        (
            "protein|GCF_000009999.1_prot_731 elongation factor Tu [protein_id=NP_414543.1]",
            "NP_414543.1",
            "GCF_000009999.1 elongation factor Tu [protein_id=NP_414543.1]",
        ),
    ],
)
def test_get_seq_id_and_description(
    old_description,
    expected_seq_id,
    expected_description,
):
    seq_id, description = _get_seq_id_and_description_protein_seqs(old_description)

    assert seq_id == expected_seq_id
    assert description == expected_description


@pytest.mark.parametrize(
    "old_description",
    [
        "protein|genome_123 hypothetical protein [protein_id=ABC123]",
        "protein|genome_123_prot_456 hypothetical protein",
        "protein|genome_123 hypothetical protein",
        "",
    ],
)
def test_get_seq_id_and_description_invalid_description(old_description):
    with pytest.raises(AssertionError):
        _get_seq_id_and_description_protein_seqs(old_description)


def test_download_fasta_success():
    with patch(
        "modules.download_seqs.run.utils.run_cmd",
        return_value=">seq1\nATGC\n",
    ) as mock_run_cmd:
        result = _download_fasta("efetch command")

    assert result == ">seq1\nATGC\n"
    mock_run_cmd.assert_called_once_with("efetch command")


def test_download_fasta_retry_then_success():
    with (
        patch(
            "modules.download_seqs.run.utils.run_cmd",
            side_effect=["error", ">seq1\nATGC\n"],
        ) as mock_run_cmd,
        patch("modules.download_seqs.run.time.sleep") as mock_sleep,
    ):
        result = _download_fasta("efetch command")

    assert result == ">seq1\nATGC\n"
    assert mock_run_cmd.call_count == 2
    mock_sleep.assert_called_once_with(3)


def test_download_fasta_fails_three_times():
    with (
        patch(
            "modules.download_seqs.run.utils.run_cmd",
            side_effect=["error", "error", "error"],
        ) as mock_run_cmd,
        patch("modules.download_seqs.run.time.sleep") as mock_sleep,
    ):
        result = _download_fasta("efetch command")

    assert result is None
    assert mock_run_cmd.call_count == 3
    assert mock_sleep.call_count == 2
    assert mock_sleep.call_args_list == [
        call(3),
        call(6),
    ]
