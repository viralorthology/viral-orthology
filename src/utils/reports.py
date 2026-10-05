import statistics
from collections.abc import Iterator

import utils
from config.context import Context
from models.seq import Seq


def make_og_report(ctx: Context) -> None:
    """Generate a CSV report with statistics for each ortholog group."""
    lines = ["name,n_seqs,min_prot_len,max_prot_len,av_prot_len,st_dev_prot_len"]

    ogs = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)
    for og in ogs:
        min_prot_len, max_prot_len, av_prot_len, stdev = _get_prot_len_statistics(
            og.seqs
        )
        lines.append(
            f"{og.path.name},{og.n_seqs},{min_prot_len},{max_prot_len},{av_prot_len},{stdev}"
        )

    report_file = (
        ctx.paths.output_dir / "ortholog_groups_report.csv"
    )  # TODO include dataset hash in file name
    report_file.write_text(("\n").join(lines), encoding="utf-8")


def _get_prot_len_statistics(seqs: Iterator[Seq]) -> tuple[int, int, float, float]:
    """
    Calculate summary statistics for protein sequence lengths.

    Args:
        seqs: Iterator containing protein sequences.

    Returns:
        A tuple containing the minimum length, maximum length, mean length,
        and standard deviation of the protein sequences.
    """
    prot_lens = [len(seq.seq) for seq in seqs]

    return (
        min(prot_lens),
        max(prot_lens),
        round(statistics.mean(prot_lens), 2),
        round(statistics.stdev(prot_lens), 2),
    )
