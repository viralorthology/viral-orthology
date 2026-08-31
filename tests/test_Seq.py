import pytest
from Bio.Seq import Seq as BioSeq

from models.fasta_type import FastaType
from models.seq import Seq


def test_generic_seq():
    _ = Seq(
        BioSeq("ATGC"),
        "seq1",
        "seq1",
        FastaType.GENERIC,
    )


def test_protein_seq():
    _ = Seq(
        BioSeq("ATGC"),
        "seq1",
        "seq1 genome1",
        FastaType.PROTEIN,
    )


def test_protein_seq_without_genome_id():
    with pytest.raises(ValueError):
        _ = Seq(
            BioSeq("ATGC"),
            "seq1",
            "seq1",
            FastaType.PROTEIN,
        )
