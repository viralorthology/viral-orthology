import utils
from models.constants import PREDICTED_PROTEINS_PREFIX
from models.fasta import Fasta
from models.seq import Seq


def run(
    genome_fasta: Fasta,
    params: str,
    predicted_prots_n: int,
    annotated_proteome_fasta: Fasta | None = None,
) -> list[Seq]:
    """
    Predict genome proteins and return the formatted sequences.

    Predicted proteins matching annotated proteins are removed when an
    annotated proteome is provided. The remaining sequences are assigned
    predicted protein IDs and genomic location descriptions.

    Args:
        genome_fasta: Genome FASTA file used for protein prediction.
        params: Command-line parameters passed to ORFfinder.
        predicted_prots_n: Number of previously predicted proteins, used to
            generate consecutive protein IDs.
        annotated_proteome_fasta: Optional annotated proteome used to remove
            already annotated proteins.

    Raises:
        ValueError: If the genome FASTA does not contain exactly one sequence.

    Returns:
        A list of predicted protein sequences with formatted IDs and
        descriptions.
    """
    if genome_fasta.n_seqs != 1:
        raise ValueError(
            f"{genome_fasta.path} must contain exactly one sequence for protein prediction"
        )

    # predict proteins
    predicted_prots_fasta_str = _run_orffinder(genome_fasta, params)
    predicted_prot_seqs = utils.get_seqs_from_fasta_str(predicted_prots_fasta_str)

    # remove annotated prots from predicted proteome
    if annotated_proteome_fasta is not None:
        ids_prots_to_remove = []

        for predicted_prot in predicted_prot_seqs:
            predicted_prot_str = str(predicted_prot.seq)

            for annotated_prot in annotated_proteome_fasta.seqs:
                if _predicted_prot_is_annotated(
                    str(annotated_prot.seq), predicted_prot_str
                ):
                    ids_prots_to_remove.append(predicted_prot.id)
                    break

        if ids_prots_to_remove:
            predicted_prot_seqs = [
                seq for seq in predicted_prot_seqs if seq.id not in ids_prots_to_remove
            ]

    # format seq descriptions

    genome_id = next(iter(genome_fasta.ids))
    seqs_with_description_format = []
    for n, seq in enumerate(predicted_prot_seqs):
        seq_location_str = _get_predicted_prot_location(seq.id)
        seq.id = f"{PREDICTED_PROTEINS_PREFIX}_{predicted_prots_n + n + 1}"
        seq.description = f"{seq.id} {genome_id} {seq_location_str}"
        seqs_with_description_format.append(seq)

    return seqs_with_description_format


def _run_orffinder(genome_fasta: Fasta, params: str) -> str:
    """
    Run ORFfinder on a genome FASTA file and return its output.

    Raises:
        ValueError: If ORFfinder does not produce a valid FASTA output.
    """
    cmd_stdout = utils.run_cmd(f"ORFfinder -in {genome_fasta.path} {params}")
    if not cmd_stdout.startswith(">"):
        raise ValueError(f"Cannot predict proteome on {genome_fasta.path}")

    return cmd_stdout


def _predicted_prot_is_annotated(
    annotated_prot_seq: str, predicted_prot_seq: str
) -> bool:
    """Check whether either protein sequence is contained in the other."""

    if len(predicted_prot_seq) <= len(annotated_prot_seq):
        return predicted_prot_seq in annotated_prot_seq

    return annotated_prot_seq in predicted_prot_seq


def _get_predicted_prot_location(seq_id: str) -> str:
    """Extract the genomic location from a predicted sequence ID."""
    start, end = map(int, seq_id.split(":")[-2:])
    assert start != end
    position_str = ("..").join(map(str, sorted((start, end))))

    if start > end:
        position_str = f"complement({position_str})"

    return f"[location={position_str}]"
