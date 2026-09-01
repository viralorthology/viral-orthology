from pathlib import Path

from Bio.Seq import Seq as BioSeq
from Bio.SeqRecord import SeqRecord

from models.constants import PREDICTED_PROTEINS_PREFIX
from models.fasta_type import FastaType


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


def get_seq_from_seqrecord(
    seqrecord: SeqRecord, fasta_type: FastaType, fasta_path: Path
) -> Seq:
    """
    Convert a Biopython SeqRecord into a Seq model.

    Args:
        seqrecord: The Biopython record containing the sequence and metadata.
        fasta_type: The type of FASTA record being processed.
        fasta_path: Path to the FASTA file, used to provide context in error
            messages.

    Raises:
        ValueError: If the sequence is empty, the record ID is missing, or
            the protein sequence description is malformed.
    """
    if not seqrecord.seq:
        raise ValueError(f"{fasta_path} contains and empty sequence")
    if not seqrecord.id:
        raise ValueError(f"{fasta_path} contains a malformed FASTA record")

    seq = Seq(
        seq=BioSeq(seqrecord.seq),
        seq_id=seqrecord.id,
        seq_description=seqrecord.description,
    )

    match fasta_type:
        case FastaType.GENERIC:
            pass
        case FastaType.PROTEIN:
            description_elements = seq.description.split()
            if len(description_elements) < 2:
                raise ValueError(
                    f"Invalid protein sequence description: {seq.description}"
                )

            seq.genome_id = description_elements[1]

    return seq


def get_seqrecord_from_seq(seq: Seq) -> SeqRecord:
    return SeqRecord(
        seq=seq.seq,
        id=seq.id,
        description=seq.description,
    )
