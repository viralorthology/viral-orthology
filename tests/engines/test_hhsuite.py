from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from engines.hhsuite import (
    _hhalign,
    hmm_compare_groups,
)


@pytest.fixture
def group1():
    group = MagicMock()
    group.alignment_path = Path("/tmp/group1.aln.fasta")
    group.a2m_path = Path("/tmp/group1.a2m")
    group.hmm_hhsuite_path = Path("/tmp/group1.hmm")
    return group


@pytest.fixture
def group2():
    group = MagicMock()
    group.alignment_path = Path("/tmp/group2.aln.fasta")
    group.a2m_path = Path("/tmp/group2.a2m")
    group.hmm_hhsuite_path = Path("/tmp/group2.hmm")
    return group


def test_hmm_compare_groups(group1, group2):
    with (
        patch("engines.hhsuite.muscle_align"),
        patch("engines.hhsuite._convert_aligned_fasta_to_a2m"),
        patch("engines.hhsuite._build_hmm_hhsuite"),
        patch(
            "engines.hhsuite._hhalign",
            return_value=97.5,
        ),
    ):
        score = hmm_compare_groups(group1, group2)

    assert score == pytest.approx(97.5)


def test_hhalign_returns_probability_score(group1, group2):
    output = """
    No 1
    Probab=98.73 E-value=1.2e-20 Score=123.4
    """

    with patch(
        "engines.hhsuite.utils.run_cmd",
        return_value=output,
    ):
        score = _hhalign(group1, group2)

    assert score == pytest.approx(98.73)


def test_hhalign_returns_zero_when_no_hit(group1, group2):
    with patch(
        "engines.hhsuite.utils.run_cmd",
        return_value="No Hit",
    ):
        score = _hhalign(group1, group2)

    assert score == 0.0
