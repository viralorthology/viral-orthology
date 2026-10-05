import pytest

from modules.pipeline.merge_by_annotation import _ogs_can_be_merged, _skip_og


@pytest.mark.parametrize(
    ("annotation", "n_seqs", "active_genomes", "expected"),
    [
        ("hypothetical-protein", 1, 5, True),
        ("ABC transporter", 5, 5, True),
        ("ABC transporter", 3, 5, False),
        ("hypothetical-protein", 5, 5, True),
    ],
)
def test_skip_og(annotation, n_seqs, active_genomes, expected):
    assert _skip_og(annotation, n_seqs, active_genomes) is expected


@pytest.mark.parametrize(
    (
        "og1_annotation",
        "og2_annotation",
        "og1_genome_ids",
        "og2_genome_ids",
        "expected",
    ),
    [
        (
            "ABC transporter",
            "ABC transporter",
            {"genome1", "genome2"},
            {"genome3", "genome4"},
            True,
        ),
        (
            "ABC transporter",
            "DNA-binding protein",
            {"genome1", "genome2"},
            {"genome3", "genome4"},
            False,
        ),
        (
            "ABC transporter",
            "ABC transporter",
            {"genome1", "genome2"},
            {"genome2", "genome3"},
            False,
        ),
        (
            "ABC transporter",
            "DNA-binding protein",
            {"genome1", "genome2"},
            {"genome2", "genome3"},
            False,
        ),
    ],
)
def test_ogs_can_be_merged(
    og1_annotation,
    og2_annotation,
    og1_genome_ids,
    og2_genome_ids,
    expected,
):
    assert (
        _ogs_can_be_merged(
            og1_annotation,
            og2_annotation,
            og1_genome_ids,
            og2_genome_ids,
        )
        is expected
    )
