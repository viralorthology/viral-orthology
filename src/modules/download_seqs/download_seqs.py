import time

from Bio.SeqUtils import gc_fraction

import utils
from cli.ui import UI
from config.context import Context
from models.fasta import Fasta
from models.seq import Seq


def run(ctx: Context) -> None:
    ##### VALIDATION #####
    utils.check_dependencies("efetch")
    utils.ensure_files_have_content(ctx.paths.ids_txt)
    utils.ensure_paths_do_not_exist(ctx.paths.genomes_fasta, ctx.paths.proteomes_fasta)
    ctx.paths.output_dir.mkdir(parents=True, exist_ok=True)

    ##### RUN #####
    genome_ids_to_download = _get_genome_ids_to_download(
        ctx.paths.ids_txt.read_text(encoding="utf-8")
    )

    # download seqs
    ctx.ui.show("Downloading sequences from GenBank...")

    genomes_fasta = Fasta(ctx.paths.genomes_fasta)
    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)

    genome_ids_failed, _ = (
        _download_sequences(  # TODO predict proteins if no proteome could be downloaded
            ctx.ui, genome_ids_to_download, genomes_fasta, proteomes_fasta
        )
    )

    for genome_id in genome_ids_failed:
        ctx.ui.show_error(f"{genome_id} could not be downloaded")
    if len(genome_ids_failed) == len(genome_ids_to_download):
        raise ValueError("No genome could be downloaded")

    # format protein sequence descriptions
    all_protein_seqs = []
    for seq in proteomes_fasta.seqs:
        seq.id, seq.description = _get_seq_id_description_from_gb_description(
            seq.description
        )
        all_protein_seqs.append(seq)

    proteomes_fasta.delete_fasta()
    proteomes_fasta.add_seqs(*all_protein_seqs)

    ##### MAKE REPORT #####
    ctx.ui.show("Running sequence analysis...")

    genome_seqs = list(genomes_fasta.seqs)
    genome_lens, genome_gc_perc, genome_n_counts, identical_genomes = _analyze_genomes(
        genome_seqs
    )

    for genome1_id, genome2_id in identical_genomes:
        ctx.ui.show(f"{genome1_id} and {genome2_id} have an identical sequence.")

    genome_ids = [genome.id for genome in genome_seqs]
    dataset_hash = utils.get_dataset_hash(genome_ids)
    genome_report_path = ctx.paths.output_dir / f"dataset_{dataset_hash}.csv"

    report_content = ["genome_id,genome_length,gc_perc,Ns_in_genome"]
    for genome_id in sorted(genome_ids):
        length = genome_lens[genome_id]
        gc_perc = genome_gc_perc[genome_id]
        ns = genome_n_counts[genome_id]
        report_content.append(f"{genome_id},{length},{gc_perc},{ns}")

    genome_report_path.write_text(("\n").join(report_content))


def _download_sequences(
    ui: UI, genome_ids: set[str], genomes_fasta: Fasta, proteomes_fasta: Fasta
) -> tuple[list[str], list[str]]:
    """
    Download genome and proteome sequences for a set of genome IDs.

    Genome and proteome sequences are downloaded from GenBank using
    ``efetch``. Genome download failures are recorded and the corresponding
    genome is skipped. Proteome download failures are also recorded, but do
    not prevent the successfully downloaded genome sequence from being saved.

    Returns:
        - A list of genome ids for genomes that could not be
            downloaded.
        - A list of genome IDs for which no proteome could be downloaded.
    """
    assert genome_ids

    genome_ids_failed = []
    genome_ids_without_proteome = []

    genome_fastas_to_write = []
    proteome_fastas_to_write = []

    for genome_id in ui.progress_bar(genome_ids):
        genome_fasta_str = _download_fasta(
            f'efetch -db nuccore -id "{genome_id}" -format fasta'
        )
        if genome_fasta_str is None:
            genome_ids_failed.append(genome_id)
            continue

        genome_fastas_to_write.append(genome_fasta_str.rstrip("\n") + "\n")

        proteome_fasta_str = _download_fasta(
            f'efetch -db nuccore -id "{genome_id}" -format fasta_cds_aa'
        )

        if proteome_fasta_str is None:
            genome_ids_without_proteome.append(genome_id)
            continue

        proteome_fastas_to_write.append(proteome_fasta_str.rstrip("\n") + "\n")

    if genome_fastas_to_write:
        with genomes_fasta.path.open("a", encoding="utf-8") as fh:
            fh.write(("").join(genome_fastas_to_write))

    if proteome_fastas_to_write:
        with proteomes_fasta.path.open("a", encoding="utf-8") as fh:
            fh.write(("").join(proteome_fastas_to_write))

    return genome_ids_failed, genome_ids_without_proteome


def _analyze_genomes(
    genome_seqs: list[Seq],
) -> tuple[dict[str, int], dict[str, float], dict[str, int], list[tuple[str, str]]]:
    """
    Calculate statistics and identify identical genome sequences.

    For each genome, calculates its sequence length, GC content percentage,
    and number of ambiguous N nucleotides. Also identifies pairs of
    genomes with identical sequences.

    Args:
        genome_seqs: List of genome sequences to analyze.

    Returns:
        - A mapping of genome IDs to sequence lengths.
        - A mapping of genome IDs to GC content percentages.
        - A mapping of genome IDs to the number of ``N`` nucleotides.
        - A list of pairs of genome IDs with identical sequences.
    """
    assert len([seq.id for seq in genome_seqs]) == len({seq.id for seq in genome_seqs})
    genome_n_counts = {}
    genome_lens = {}
    genome_gc_perc = {}
    identical_genomes = []

    for n, genome1 in enumerate(genome_seqs):
        genome_lens[genome1.id] = len(str(genome1.seq))
        genome_gc_perc[genome1.id] = round(gc_fraction(genome1.seq) * 100, 2)
        genome_n_counts[genome1.id] = str(genome1.seq).count("N")

        for genome2 in genome_seqs[n + 1 :]:
            if genome1.seq == genome2.seq:
                identical_genomes.append((genome1.id, genome2.id))

    return genome_lens, genome_gc_perc, genome_n_counts, identical_genomes


def _get_seq_id_description_from_gb_description(
    old_description: str,
) -> tuple[str, str]:
    """
    Extract the protein sequence ID and reconstruct its description.

    The input description is expected to contain a [protein_id=...]
    field and a genome identifier in the format <genome_id>_prot_...

    Args:
        old_description: Original GenBank sequence description.

    Returns:
        - The protein sequence ID
        - The reconstructed description prefixed with the genome ID.
    """
    if "_prot_" not in old_description or "[protein_id=" not in old_description:
        raise ValueError(
            f"Protein description does not contain the required data: {old_description}"
        )

    seq_id = old_description.split("[protein_id=")[1].split("]")[0]
    genome_id = old_description.split("|")[1].split("_prot_")[0]
    remaining_description = (" ").join(old_description.split()[1:])
    return seq_id, f"{seq_id} {genome_id} {remaining_description}"


def _download_fasta(cmd_str: str) -> str | None:
    """
    Run a command to download sequences from GenBank, retrying on failure.

    The command is attempted up to three times. A successful response is
    identified by content starting with >. Failed attempts are retried
    with an increasing delay between attempts.

    Args:
        cmd_str: Command used to download the FASTA.

    Returns:
        The downloaded FASTA content if successful, otherwise None.
    """
    assert cmd_str

    sleep_time = 3

    for n_attempt in range(3):
        fasta_str = utils.run_cmd(cmd_str)

        if fasta_str.startswith(">"):
            return fasta_str

        if n_attempt < 2:
            time.sleep((n_attempt + 1) * sleep_time)

    return None


def _get_genome_ids_to_download(ids_file_content: str) -> set[str]:
    """
    Extract unique genome IDs from file content.

    Version numbers are removed from genome IDs before duplicates are
    discarded.

    Args:
        ids_file_content: Contents of the IDs file.

    Returns:
        A set of unique genome IDs without version numbers.
    """
    assert ids_file_content

    genome_ids_to_download = set()
    for line in ids_file_content.splitlines():
        if not line.strip():
            continue
        genome_id = line.split(".")[0]
        genome_ids_to_download.add(genome_id)

    return genome_ids_to_download
