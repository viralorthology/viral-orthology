from pathlib import Path
from unittest.mock import Mock

from modules.pipeline.merge_by_annotation import _get_same_annotation_ogs


def make_og(stem: str, exists: bool = True) -> Mock:
    og = Mock()
    og.path = Mock(spec=Path)
    og.path.stem = stem
    og.path.exists.return_value = exists
    return og


def test_returns_ogs_with_same_annotation():
    og1 = make_og("annotation-1")
    og2 = make_og("annotation-2")
    og3 = make_og("other-1")

    result = _get_same_annotation_ogs(
        "annotation",
        [og1, og2, og3],
    )

    assert result == [og1, og2]


def test_ignores_ogs_whose_fasta_does_not_exist():
    og1 = make_og("annotation-1")
    og2 = make_og("annotation-2", exists=False)

    result = _get_same_annotation_ogs(
        "annotation",
        [og1, og2],
    )

    assert result == [og1]


def test_returns_empty_list_when_no_og_has_same_annotation():
    og1 = make_og("other-1")
    og2 = make_og("another-1")

    result = _get_same_annotation_ogs(
        "annotation",
        [og1, og2],
    )

    assert result == []


def test_uses_last_hyphen_to_determine_annotation():
    og1 = make_og("annotation-with-hyphen-1")
    og2 = make_og("annotation-with-hyphen-2")
    og3 = make_og("annotation-with-hyphen-other-1")

    result = _get_same_annotation_ogs(
        "annotation-with-hyphen",
        [og1, og2, og3],
    )

    assert result == [og1, og2]
