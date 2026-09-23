from modules.pipeline.rename_og_fastas_by_annotation import (
    _add_number_to_repeated_name,
    _filter_annotations,
    _get_most_frequent_annotation,
    _get_seq_annotation,
    _sanitize_annotation,
)


def test_filter_annotations():
    annotations = [
        "kinase",
        None,
        "unknown protein",
        "hypothetical protein",
        "transporter",
    ]

    assert _filter_annotations(annotations) == [
        "kinase",
        "transporter",
    ]


def test_filter_annotations_empty():
    assert _filter_annotations([None, "unknown", "hypothetical"]) == []


def test_get_seq_annotation_from_protein_tag():
    description = "some protein [protein=DNA polymerase III] [gene=dnaE]"

    assert _get_seq_annotation(description) == "dna-polymerase-iii"


def test_get_seq_annotation_from_gene_tag():
    description = "some protein [gene=dnaE]"

    assert _get_seq_annotation(description) == "dnae"


def test_get_seq_annotation_returns_none_when_no_annotation():
    description = "some protein without annotation"

    assert _get_seq_annotation(description) is None


def test_get_seq_annotation_from_protein():
    description = "[protein=DNA polymerase III]"

    assert _get_seq_annotation(description) == "dna-polymerase-iii"


def test_sanitize_annotation():
    assert _sanitize_annotation("DNA polymerase III") == "DNA-polymerase-III"


def test_sanitize_annotation_removes_leading_and_trailing_hyphens():
    assert _sanitize_annotation("---DNA polymerase III---") == "DNA-polymerase-III"


def test_sanitize_annotation_replaces_special_characters():
    assert _sanitize_annotation("ABC/DEF (protein)") == "ABC-DEF-protein"


def test_get_most_frequent_annotation():
    annotations = [
        "kinase",
        "transporter",
        "kinase",
        "receptor",
        "kinase",
    ]

    assert _get_most_frequent_annotation(annotations) == ("kinase", 3)


def test_get_most_frequent_annotation_breaks_ties_alphabetically():
    annotations = [
        "transporter",
        "kinase",
        "transporter",
        "kinase",
    ]

    assert _get_most_frequent_annotation(annotations) == ("kinase", 2)


def test_get_most_frequent_annotation_single_annotation():
    assert _get_most_frequent_annotation(["kinase"]) == ("kinase", 1)


def test_add_number_to_repeated_name():
    used_names = {"kinase"}

    assert _add_number_to_repeated_name("kinase", used_names) == "kinase-1"


def test_add_number_to_repeated_name_skips_existing_suffixes():
    used_names = {
        "kinase",
        "kinase-1",
        "kinase-2",
    }

    assert _add_number_to_repeated_name("kinase", used_names) == "kinase-3"


def test_add_number_to_repeated_name_returns_first_available_suffix():
    used_names = {
        "kinase",
        "kinase-1",
        "kinase-3",
    }

    assert _add_number_to_repeated_name("kinase", used_names) == "kinase-2"
