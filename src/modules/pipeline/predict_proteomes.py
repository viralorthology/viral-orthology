import utils
from config.context import Context
from engines import predict_proteome
from models.fasta import Fasta
from models.fasta_type import FastaType


def predict_proteomes(ctx: Context) -> None:
    """Predict proteomes from genome sequences and save the resulting proteomes."""
    genome_fastas = utils.get_fastas(
        ctx.paths.sequences_dir, FastaType.GENERIC, ".genome"
    )
    n_predicted_proteins = 0  # TODO if -download_seqs predicts proteome, change this
    for genome_fasta in ctx.ui.progress_bar(genome_fastas):
        genome_id = genome_fasta.ids[0]
        proteome_fasta = Fasta(
            ctx.paths.sequences_dir / f"{genome_id}.proteome", FastaType.PROTEIN
        )
        predicted_proteins = predict_proteome(
            genome_fasta,
            ctx.args.tool_args["orffinder"],
            n_predicted_proteins,
            proteome_fasta,
        )
        predicted_fasta = Fasta(
            ctx.paths.sequences_dir / f"{genome_id}.predicted", FastaType.GENERIC
        )
        predicted_fasta.add_seqs(*predicted_proteins)
        n_predicted_proteins += len(predicted_proteins)
