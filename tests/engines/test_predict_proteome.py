import pytest

from engines.orffinder import (
    _get_predicted_prot_location,
    _predicted_prot_is_annotated,
)


@pytest.mark.parametrize(
    "annotated, predicted, expected",
    [
        ("MKKLL", "MKKLL", True),
        ("MKKLLAA", "KKLL", True),
        ("KKLL", "MKKLLAA", True),
        ("MKKLL", "AAAA", False),
    ],
)
def test_predicted_prot_is_annotated(annotated, predicted, expected):
    assert (
        _predicted_prot_is_annotated(
            annotated,
            predicted,
        )
        == expected
    )


@pytest.mark.parametrize(
    "seq_id, expected",
    [
        ("genome:100:200", "[location=100..200]"),
        ("genome:200:100", "[location=complement(100..200)]"),
    ],
)
def test_get_predicted_prot_location(seq_id, expected):
    assert _get_predicted_prot_location(seq_id) == expected
