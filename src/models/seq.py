from enum import Enum, auto

from Bio.Seq import Seq as BioSeq
from Bio.SeqRecord import SeqRecord

from models.constants import PREDICTED_PROTEINS_PREFIX


class SeqType(Enum):
    GENERIC = auto()
    PROTEIN = auto()


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
        self.seq_type = self._get_seq_type()
        self._add_metadata()

    def _get_seq_type(self) -> SeqType:
        if "[protein_id=" in self.description:
            return SeqType.PROTEIN
        return SeqType.GENERIC

    def _add_metadata(self) -> None:
        if self.seq_type == SeqType.PROTEIN:
            self.genome_id = self.description.split()[1]
            self.is_predicted = self.id.startswith(PREDICTED_PROTEINS_PREFIX)


def get_seq_from_seqrecord(seqrecord: SeqRecord) -> Seq:
    """Convert a Biopython SeqRecord into a Seq model."""
    if not seqrecord.seq:
        raise ValueError("Empty sequence found.")
    if not seqrecord.id:
        raise ValueError("Malformed record found.")

    return Seq(
        seq=BioSeq(seqrecord.seq),
        seq_id=seqrecord.id,
        seq_description=seqrecord.description,
    )


def get_seqrecord_from_seq(seq: Seq) -> SeqRecord:
    return SeqRecord(
        seq=seq.seq,
        id=seq.id,
        description=seq.description,
    )
