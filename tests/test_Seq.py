import pytest

from classes.fasta_type import FastaType
from classes.seq import Seq


def test_generic_seq(seq_without_genome_id):
    _ = Seq(seq_without_genome_id, FastaType.GENERIC)


def test_protein_seq(seq_with_genome_id):
    _ = Seq(seq_with_genome_id, FastaType.PROTEIN)


def test_protein_seq_without_genome_id(seq_without_genome_id):
    with pytest.raises(ValueError):
        _ = Seq(seq_without_genome_id, FastaType.PROTEIN)
