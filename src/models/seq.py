from Bio.Seq import Seq as BioSeq
from Bio.SeqRecord import SeqRecord

from models.constants import PREDICTED_PROTEINS_PREFIX


class Seq:
    genome_id: str
    """
    Represents a biological sequence and its associated metadata.

    Args:
            seq: The biological sequence as a Biopython Seq object.
            seq_id: The sequence identifier.
            seq_description: The FASTA record description.
    """

    def __init__(self, seq: BioSeq, seq_id: str, seq_description: str):
        self.seq = seq
        self.id = seq_id
        self.description = seq_description
        self.is_predicted = self.id.startswith(PREDICTED_PROTEINS_PREFIX)


def get_seq_from_seqrecord(seqrecord: SeqRecord, seq_id: str) -> Seq:
    """Convert a Biopython SeqRecord into a Seq model."""

    seq = Seq(
        seq=BioSeq(seqrecord.seq),
        seq_id=seq_id,
        seq_description=seqrecord.description,
    )

    if is_protein_seq(seq.description):
        seq.genome_id = seq.description.split()[1]

    return seq


def get_seqrecord_from_seq(seq: Seq) -> SeqRecord:
    return SeqRecord(
        seq=seq.seq,
        id=seq.id,
        description=seq.description,
    )


def is_protein_seq(seq_description: str) -> bool:
    return "[protein_id=" in seq_description
