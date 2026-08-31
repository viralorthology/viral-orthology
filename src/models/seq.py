from Bio.Seq import Seq as BioSeq

from models.constants import PREDICTED_PROTEINS_PREFIX
from models.fasta_type import FastaType


class Seq:
    genome_id: str

    def __init__(
        self, seq: BioSeq, seq_id: str, seq_description: str, fasta_type: FastaType
    ):
        self.seq = seq
        self.id = seq_id
        self.description = seq_description
        self.is_predicted = self.id.startswith(PREDICTED_PROTEINS_PREFIX)
        self._add_data(fasta_type)

    def _add_data(self, fasta_type: FastaType) -> None:
        match fasta_type:
            case FastaType.GENERIC:
                pass
            case FastaType.PROTEIN:
                self._add_protein_data()

    def _add_protein_data(self) -> None:
        description_elements = self.description.split()
        if len(description_elements) < 2:
            raise ValueError(
                f"Invalid protein sequence description: {self.description}"
            )

        self.genome_id = description_elements[1]
