import pytest
from Bio.Seq import Seq as BioSeq

from models.seq import Seq


def test_generic_seq():
    _ = Seq(
        BioSeq("ATGC"),
        "seq1",
        "seq1",
    )


def test_protein_seq():
    _ = Seq(
        BioSeq("ATGC"),
        "seq1",
        "seq1 genome1",
    )


def test_protein_seq_without_genome_id():
    with pytest.raises(AttributeError):
        seq = Seq(
            BioSeq("ATGC"),
            "seq1",
            "seq1",
        )
        print(seq.genome_id)
