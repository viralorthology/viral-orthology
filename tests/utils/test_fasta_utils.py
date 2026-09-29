from unittest.mock import Mock

import pytest

from utils import o_groups_are_compatible


def make_og(*genome_ids: str) -> Mock:
    og = Mock()
    og.genome_ids = list(genome_ids)
    return og


@pytest.mark.parametrize(
    "ortholog_groups",
    [
        (),
        (make_og("genome1"),),
    ],
)
def test_raises_if_fewer_than_two_ortholog_groups_are_provided(
    ortholog_groups,
):
    with pytest.raises(
        ValueError,
        match="At least two ortholog groups must be provided",
    ):
        o_groups_are_compatible(*ortholog_groups)


@pytest.mark.parametrize(
    "ortholog_groups, expected",
    [
        (
            (make_og("genome1"), make_og("genome2")),
            True,
        ),
        (
            (make_og("genome1", "genome2"), make_og("genome3")),
            True,
        ),
        (
            (make_og("genome1", "genome2"), make_og("genome2", "genome3")),
            False,
        ),
        (
            (
                make_og("genome1"),
                make_og("genome2"),
                make_og("genome3"),
            ),
            True,
        ),
        (
            (
                make_og("genome1"),
                make_og("genome2"),
                make_og("genome1", "genome3"),
            ),
            False,
        ),
    ],
)
def test_o_groups_are_compatible(ortholog_groups, expected):
    assert o_groups_are_compatible(*ortholog_groups) is expected
