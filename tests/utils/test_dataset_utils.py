import pytest

from utils import get_dataset_hash


def test_get_dataset_hash_is_deterministic():
    assert get_dataset_hash(["genome_a", "genome_b"]) == get_dataset_hash(
        ["genome_a", "genome_b"]
    )


def test_get_dataset_hash_is_independent_of_input_order():
    assert get_dataset_hash(["genome_a", "genome_b"]) == get_dataset_hash(
        ["genome_b", "genome_a"]
    )


def test_get_dataset_hash_empty_input():
    with pytest.raises(ValueError):
        get_dataset_hash([])


def test_get_dataset_hash_duplicate_ids():
    with pytest.raises(ValueError):
        get_dataset_hash(["genome_a", "genome_a"])


def test_get_dataset_hash_returns_short_hash():
    result = get_dataset_hash(["genome_a", "genome_b"])

    assert len(result) == 8
    assert all(c in "0123456789abcdef" for c in result)
