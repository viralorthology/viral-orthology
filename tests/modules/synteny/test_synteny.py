import statistics

import pytest

from modules.synteny.synteny import _conserved_sw_is_compatible, _sws_are_compatible
from modules.synteny.synteny_windows import (
    SyntenyWindow,
    _classify_sw,
    _get_pairs_of_contiguos_proteins,
    _get_seq_mid_position,
    _get_seq_strand,
)


@pytest.mark.parametrize(
    ("sw1", "sw2", "expected"),
    [
        (
            SyntenyWindow("1", "2", True, False),
            SyntenyWindow("1", "2", True, False),
            True,
        ),
        (
            SyntenyWindow("1", "2", True, True),
            SyntenyWindow("2", "1", False, False),
            True,
        ),
        (
            SyntenyWindow("1", "2", True, False),
            SyntenyWindow("2", "1", False, False),
            False,
        ),
        (
            SyntenyWindow("1", "2", True, False),
            SyntenyWindow("1", "3", True, False),
            False,
        ),
        (
            SyntenyWindow("1", "2", True, False),
            SyntenyWindow("1", "2", False, True),
            False,
        ),
        (
            SyntenyWindow("1", "2", True, True),
            SyntenyWindow("3", "4", True, True),
            False,
        ),
        (
            SyntenyWindow("1", "2", True, True),
            SyntenyWindow("3", "4", False, False),
            False,
        ),
    ],
)
def test_sws_are_compatible(sw1, sw2, expected):
    assert (
        _sws_are_compatible(
            sw1,
            sw2,
        )
        is expected
    )


@pytest.mark.parametrize(
    (
        "incomplete_sw",
        "conserved_sw",
        "expected",
    ),
    [
        (
            SyntenyWindow("1", "2", True, True),
            SyntenyWindow("1", "3", True, True),
            True,
        ),
        (
            SyntenyWindow("1", "2", True, True),
            SyntenyWindow("3", "1", False, False),
            True,
        ),
        (
            SyntenyWindow("1", "2", True, False),
            SyntenyWindow("3", "1", True, False),
            True,
        ),
        (
            SyntenyWindow("1", "2", True, False),
            SyntenyWindow("1", "3", True, False),
            True,
        ),
        (
            SyntenyWindow("1", "2", True, True),
            SyntenyWindow("3", "4", True, True),
            False,
        ),
        (
            SyntenyWindow("1", "2", True, False),
            SyntenyWindow("1", "3", False, True),
            False,
        ),
    ],
)
def test_conserved_sw_is_compatible(
    incomplete_sw,
    conserved_sw,
    expected,
):
    assert (
        _conserved_sw_is_compatible(
            incomplete_sw,
            conserved_sw,
        )
        is expected
    )


@pytest.mark.parametrize(
    ("prot1_id", "prot2_id", "seqs_in_og_ids", "expected"),
    [
        ("1", "2", {"1", "2"}, "complete"),
        ("1", "3", {"1", "2"}, "incomplete"),
        ("3", "4", {"1", "2"}, "discarded"),
    ],
)
def test_classify_sw(prot1_id, prot2_id, seqs_in_og_ids, expected):
    assert _classify_sw(prot1_id, prot2_id, seqs_in_og_ids) == expected


def test_get_pairs_of_contigouos_proteins():
    ids = ["1", "2", "3", "4", "5"]

    assert _get_pairs_of_contiguos_proteins(ids) == [
        ("1", "2"),
        ("2", "3"),
        ("3", "4"),
        ("4", "5"),
    ]


@pytest.mark.parametrize(
    ("seq_description_str", "expected"),
    [
        ("complement(123...123)", False),
        ("(123...123)", True),
    ],
)
def test_get_seq_strand(seq_description_str, expected):
    assert _get_seq_strand(seq_description_str) is expected


@pytest.mark.parametrize(
    ("seq_description_str", "expected"),
    [
        ("complement(123...234)", [123, 234]),
        ("123...234", [123, 234]),
        ("3132...13232", [3132, 13232]),
        ("3132...12312...23233...232323...13232", [3132, 13232]),
        ("complement(3132...12312...23233...232323...13232)", [3132, 13232]),
    ],
)
def test_get_seq_mid_position(seq_description_str, expected):
    assert _get_seq_mid_position(seq_description_str) == statistics.mean(expected)
