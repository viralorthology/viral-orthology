from unittest.mock import Mock

import pytest

from modules.pipeline.merge_by_annotation import _ogs_have_same_annotation, _skip_og


@pytest.mark.parametrize(
    "stem1, stem2",
    [
        ("pol", "pol"),
        ("pol", "pol-1"),
        ("pol", "pol-2"),
        ("pol-1", "pol-2"),
        ("pol-2", "pol-3"),
        ("pol-2", "pol"),
    ],
)
def test_ogs_have_same_annotation(stem1, stem2):
    assert _ogs_have_same_annotation(stem1, stem2)


@pytest.mark.parametrize(
    "stem1, stem2",
    [
        ("pol", "env"),
        ("pol-1", "env-1"),
        ("pol", "pol-alpha"),
        ("pol", "pol-hypothetical-protein"),
    ],
)
def test_ogs_do_not_have_same_annotation(stem1, stem2):
    assert not _ogs_have_same_annotation(stem1, stem2)


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
