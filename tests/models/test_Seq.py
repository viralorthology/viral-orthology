import pytest
from Bio.Seq import Seq as BioSeq
from Bio.SeqRecord import SeqRecord

from models.constants import PREDICTED_PROTEINS_PREFIX
from models.seq import (
    Seq,
    SeqType,
    get_seq_from_seqrecord,
    get_seqrecord_from_seq,
)


def test_seq_has_expected_attributes():
    seq = Seq(
        BioSeq("ATGC"),
        "seq1",
        "seq1 description",
    )

    assert seq.seq == BioSeq("ATGC")
    assert seq.id == "seq1"
    assert seq.description == "seq1 description"
    assert seq.seq_type == SeqType.GENERIC


def test_seq_detects_protein():
    seq = Seq(
        BioSeq("MKT"),
        "seq1",
        "seq1 genome123 [protein_id=ABC123]",
    )

    assert seq.seq_type == SeqType.PROTEIN


def test_protein_seq_adds_metadata():
    seq = Seq(
        BioSeq("MKT"),
        "seq1",
        "seq1 genome123 [protein_id=ABC123]",
    )

    assert seq.genome_id == "genome123"
    assert seq.is_predicted is False


def test_generic_seq_does_not_add_metadata():
    seq = Seq(
        BioSeq("MKT"),
        "seq1",
        "seq1 genome123",
    )

    with pytest.raises(AttributeError):
        assert seq.genome_id


def test_predicted_protein_is_detected():
    seq = Seq(
        BioSeq("MKT"),
        f"{PREDICTED_PROTEINS_PREFIX}_1",
        f"{PREDICTED_PROTEINS_PREFIX}_1 genome123 [protein_id=ABC123]",
    )

    assert seq.seq_type == SeqType.PROTEIN
    assert seq.is_predicted is True


def test_non_predicted_protein_is_detected():
    seq = Seq(
        BioSeq("MKT"),
        "protein001",
        "protein001 genome123 [protein_id=ABC123]",
    )

    assert seq.seq_type == SeqType.PROTEIN
    assert seq.is_predicted is False


def test_get_seq_from_seqrecord():
    record = SeqRecord(
        BioSeq("ATGC"),
        id="seq1",
        description="seq1 description",
    )

    seq = get_seq_from_seqrecord(record)

    assert isinstance(seq, Seq)
    assert seq.seq == BioSeq("ATGC")
    assert seq.id == "seq1"
    assert seq.description == "seq1 description"
    assert seq.seq_type == SeqType.GENERIC


def test_get_seq_from_seqrecord_detects_protein():
    record = SeqRecord(
        BioSeq("MKT"),
        id="protein1",
        description="protein1 genome123 [protein_id=ABC123]",
    )

    seq = get_seq_from_seqrecord(record)

    assert seq.seq_type == SeqType.PROTEIN
    assert seq.genome_id == "genome123"


def test_get_seq_from_seqrecord_rejects_empty_sequence():
    record = SeqRecord(
        BioSeq(""),
        id="seq1",
        description="seq1 description",
    )

    with pytest.raises(ValueError):
        get_seq_from_seqrecord(record)


def test_get_seq_from_seqrecord_rejects_missing_id():
    record = SeqRecord(
        BioSeq("ATGC"),
        id="",
        description="seq1 description",
    )

    with pytest.raises(ValueError):
        get_seq_from_seqrecord(record)


def test_get_seqrecord_from_seq():
    seq = Seq(
        BioSeq("ATGC"),
        "seq1",
        "seq1 description",
    )

    record = get_seqrecord_from_seq(seq)

    assert isinstance(record, SeqRecord)
    assert record.seq == BioSeq("ATGC")
    assert record.id == "seq1"
    assert record.description == "seq1 description"


def test_seqrecord_conversion_is_reversible():
    original = Seq(
        BioSeq("ATGC"),
        "seq1",
        "seq1 description",
    )

    record = get_seqrecord_from_seq(original)
    converted = get_seq_from_seqrecord(record)

    assert converted.seq == original.seq
    assert converted.id == original.id
    assert converted.description == original.description
    assert converted.seq_type == original.seq_type
