from Bio.SeqRecord import SeqRecord

from models.fasta_type import FastaType


class Seq(SeqRecord):
    id: str
    genome_id: str
    description: str

    def __init__(self, seq_record: SeqRecord, fasta_type: FastaType):
        super().__init__(
            seq=seq_record.seq,
            id=seq_record.id,
            description=seq_record.description,
        )
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
