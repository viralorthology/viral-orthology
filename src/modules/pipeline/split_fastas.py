from collections import defaultdict

from config.context import Context
from models.fasta import Fasta
from models.fasta_type import FastaType


def split_fastas(ctx: Context) -> None:
    """Split genomes and proteomes into one file per genome."""
    ctx.paths.sequences_dir.mkdir()

    # make a fasta for every genome
    genomes_fasta = Fasta(ctx.paths.genomes_fasta, FastaType.GENERIC)
    for genome_seq in genomes_fasta.seqs:
        genome_fasta = Fasta(
            ctx.paths.sequences_dir / f"{genome_seq.id}.genome", FastaType.GENERIC
        )
        genome_fasta.add_seqs(genome_seq)

    # make a fasta for every proteome
    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta, FastaType.PROTEIN)
    proteins_by_genome_id = defaultdict(list)

    for prot_seq in proteomes_fasta.seqs:
        proteins_by_genome_id[prot_seq.genome_id].append(prot_seq)

    for genome_id, prot_list in proteins_by_genome_id.items():
        proteome_fasta = Fasta(
            ctx.paths.sequences_dir / f"{genome_id}.proteome", FastaType.GENERIC
        )
        proteome_fasta.add_seqs(*prot_list)
