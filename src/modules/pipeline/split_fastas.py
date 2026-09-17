from collections import defaultdict

from config.context import Context
from models.fasta import Fasta


def split_fastas(ctx: Context) -> None:
    """Split proteomes into one file per genome."""
    ctx.paths.proteomes_dir.mkdir()

    # make a fasta for every proteome
    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)
    proteins_by_genome_id = defaultdict(list)

    for prot_seq in proteomes_fasta.seqs:
        proteins_by_genome_id[prot_seq.genome_id].append(prot_seq)

    for genome_id, prot_list in proteins_by_genome_id.items():
        proteome_fasta = Fasta(ctx.paths.proteomes_dir / f"{genome_id}.fasta")
        proteome_fasta.add_seqs(*prot_list)
