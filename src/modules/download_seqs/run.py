import time

from Bio.SeqUtils import gc_fraction

import utils
from cli.ui import UI
from config.context import Context
from models.fasta import Fasta
from models.fasta_type import FastaType
from models.seq import Seq


def run(ctx: Context) -> None:
    ##### VALIDATION #####
    utils.check_dependencies("efetch")
    utils.ensure_files_have_content(ctx.paths.ids_txt)
    utils.ensure_files_do_not_exist(ctx.paths.genomes_fasta, ctx.paths.proteomes_fasta)
    ctx.paths.output_dir.mkdir(parents=True, exist_ok=True)

    ##### RUN #####
    genome_ids_to_download = _get_genome_ids_to_download(
        ctx.paths.ids_txt.read_text(encoding="utf-8")
    )

    # download seqs
    ctx.ui.show("Downloading sequences from GenBank...")

    genomes_fasta = Fasta(ctx.paths.genomes_fasta, FastaType.GENERIC)
    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta, FastaType.GENERIC)

    errors = _download_sequences(
        ctx.ui, genome_ids_to_download, genomes_fasta, proteomes_fasta
    )

    for error_msg in errors:
        ctx.ui.show_error(error_msg)
    if len(errors) == len(genome_ids_to_download):
        raise ValueError("No genome could be downloaded")

    # format sequence descriptions
    all_protein_seqs = []
    for seq in proteomes_fasta.seqs:
        seq.id, seq.description = _get_seq_id_and_description(seq.description)
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


def _get_seq_id_and_description(old_description: str) -> tuple[str, str]:
    """
    Extract the protein sequence ID and reconstruct its description.

    The input description is expected to contain a [protein_id=...]
    field and a genome identifier in the format <genome_id>_prot_...

    Args:
        old_description: Original GenBank sequence description.

    Returns:
        A tuple containing the protein sequence ID and the reconstructed
        description prefixed with the genome ID.
    """
    assert "_prot_" in old_description
    assert "[protein_id=" in old_description

    seq_id = old_description.split("[protein_id=")[1].split("]")[0]
    genome_id = old_description.split("|")[1].split("_prot_")[0]
    remaining_description = (" ").join(old_description.split()[1:])
    return seq_id, f"{genome_id} {remaining_description}"


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
    genome_ids_to_download = set()
    for line in ids_file_content.splitlines():
        if not line.strip():
            continue
        genome_id = line.split(".")[0]
        genome_ids_to_download.add(genome_id)

    return genome_ids_to_download


def _download_sequences(
    ui: UI, genome_ids: set[str], genomes_fasta: Fasta, proteomes_fasta: Fasta
) -> list[str]:
    """
    Download genome and proteome sequences for a set of genome IDs.

    Genome and proteome sequences are downloaded from GenBank using efetch.
    Genome download failures are recorded as errors, while proteome download
    failures are skipped.

    Args:
        ui: UI instance used to display download progress.
        genome_ids: Set of genome IDs to download.
        genomes_fasta: Fasta object used to write downloaded genome sequences.
        proteomes_fasta: Fasta object used to write downloaded protein sequences.

    Returns:
        A list of error messages for genomes that could not be downloaded.
    """
    errors = []
    for genome_id in ui.progress_bar(genome_ids):
        genome_fasta_str = _download_fasta(
            f'efetch -db nuccore -id "{genome_id}" -format fasta'
        )
        if genome_fasta_str is None:
            errors.append(f"genome {genome_id} could not be downloaded")
            continue

        with genomes_fasta.path.open("a", encoding="utf-8") as fh:
            fh.write(genome_fasta_str.rstrip("\n") + "\n")

        proteome_fasta_str = _download_fasta(
            f'efetch -db nuccore -id "{genome_id}" -format fasta_cds_aa'
        )

        # TODO if proteome did not download, annotate with orffinder

        if proteome_fasta_str is not None:
            with proteomes_fasta.path.open("a", encoding="utf-8") as fh:
                fh.write(proteome_fasta_str.rstrip("\n") + "\n")

    return errors


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
