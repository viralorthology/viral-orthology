from unittest.mock import Mock

import pytest

from modules.pipeline.merge_by_annotation import _skip_og


@pytest.mark.parametrize(
    "stem, n_seqs, active_genomes_count",
    [
        ("abc", 3, 3),
        ("abc", 2, 2),
        ("abc-hypothetical-protein", 1, 3),
    ],
)
def test_skip_og(stem, n_seqs, active_genomes_count):
    og = Mock()
    og.path.exists.return_value = True
    og.path.stem = stem
    og.n_seqs = n_seqs

    assert _skip_og(og, active_genomes_count)


def test_skip_og_when_path_does_not_exist():
    og = Mock()
    og.path.exists.return_value = False

    assert _skip_og(og, active_genomes_count=3)


@pytest.mark.parametrize(
    "stem, n_seqs, active_genomes_count",
    [
        ("abc", 2, 3),
        ("abc", 1, 3),
    ],
)
def test_do_not_skip_og(stem, n_seqs, active_genomes_count):
    og = Mock()
    og.path.exists.return_value = True
    og.path.stem = stem
    og.n_seqs = n_seqs

    assert not _skip_og(og, active_genomes_count)
