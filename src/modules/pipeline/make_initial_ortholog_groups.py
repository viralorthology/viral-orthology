import logging
import tempfile
from collections import Counter
from pathlib import Path

import utils
from config.context import Context
from models.fasta import Fasta

logger = logging.getLogger(__name__)


def make_initial_ortholog_groups(ctx: Context) -> None:
    """
    Create initial ortholog groups from the non-redundant proteomes.

    Proteomes are compared with ProteinOrtho to generate ortholog groups.
    Groups containing paralogs are filtered by retaining the longest protein
    from each genome. Groups containing proteins from only one genome are
    deleted.
    """
    ctx.ui.show("Creating initial ortholog groups...")

    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)

    with tempfile.TemporaryDirectory() as tmp_path:
        tmp_dir = Path(tmp_path)
        proteome_fastas = _get_proteome_fastas(
            tmp_dir,
            proteomes_fasta,
            ctx.runtime.active_genome_ids,
            ctx.runtime.paralog_ids,
        )

        # run proteinortho
        og_tmp_dir = tmp_dir / "ortholog_groups"
        og_tmp_dir.mkdir()
        _run_proteinortho(
            proteome_fastas, ctx.args.tool_args["proteinortho"], og_tmp_dir
        )

        # remove paralogs from ortholog groups
        ortholog_groups = utils.get_fastas(og_tmp_dir)
        assert ortholog_groups
        deleted_ogs = set()
        for og in ortholog_groups:
            og_was_deleted = _remove_paralogs(og, ctx.runtime.paralog_ids)
            if og_was_deleted:
                deleted_ogs.add(og)

        # move ortholog groups
        ctx.paths.ortholog_groups_dir.mkdir()
        for og in ortholog_groups:
            if og in deleted_ogs:
                continue
            og.move_fasta(ctx.paths.ortholog_groups_dir)


def _get_proteome_fastas(
    tmp_dir: Path,
    proteomes_fasta: Fasta,
    active_genome_ids: set[str],
    paralog_ids: set[str],
) -> list[Fasta]:
    """
    Create temporary FASTA files containing the proteome of each active genome.

    Proteins identified as paralogs are excluded from the generated proteomes.
    Genomes without proteins are skipped.
    """
    tmp_proteome_fastas = []
    for genome_id in active_genome_ids:
        proteome_seqs = [
            seq
            for seq in proteomes_fasta.seqs
            if seq.genome_id == genome_id and seq.id not in paralog_ids
        ]
        if not proteome_seqs:
            continue

        proteome_fasta = Fasta(tmp_dir / f"{genome_id}.fasta")
        proteome_fasta.add_seqs(*proteome_seqs)
        tmp_proteome_fastas.append(proteome_fasta)

    return tmp_proteome_fastas


def _run_proteinortho(proteome_fastas: list[Fasta], params: str, cwd: Path) -> None:
    """
    Run ProteinOrtho and extract the resulting ortholog group proteins.

    Args:
        proteome_fastas: Proteome FASTA files to compare.
        params: Additional parameters to pass to ProteinOrtho.
        cwd: Working directory for ProteinOrtho and its output files.
    """
    fasta_paths = (" ").join(str(fasta.path) for fasta in proteome_fastas)

    utils.run_cmd(
        f"proteinortho {fasta_paths} {params}",
        cwd=cwd,
    )
    utils.run_cmd(
        f"proteinortho_grab_proteins.pl -tofiles myproject.proteinortho.tsv -exact {fasta_paths}",
        cwd=cwd,
    )


def _remove_paralogs(og: Fasta, paralog_ids: set[str]) -> bool:
    """
    Remove paralogous sequences from an ortholog group.

    If all sequences in the ortholog group belong to the same genome, the
    ortholog group is deleted. Otherwise, for each genome containing multiple
    sequences, the longest sequence is retained and the remaining sequences
    are removed and added to the set of paralog IDs. Sequence IDs are used
    as a deterministic tie-breaker when sequences have the same length.

    Args:
        og: FASTA file representing an ortholog group.
        paralog_ids: Set of sequence IDs identified as paralogs.

    Returns:
        True if the ortholog group was deleted because all sequences belonged
        to the same genome, otherwise False.
    """
    if len(set(og.genome_ids)) == 1:  # all seqs from the same genome
        og.delete_fasta()
        logger.debug(
            "Ortholog group %s deleted: all seqs are from the same genome", og.path
        )
        return True

    for genome_id in _get_genome_ids_with_paralogs(og.genome_ids):
        paralog_seqs = [seq for seq in og.seqs if seq.genome_id == genome_id]
        assert len(paralog_seqs) >= 2

        selected_paralog = max(
            paralog_seqs,
            key=lambda seq: (len(seq.seq), seq.id),  # Ensure reproducibility
        )

        ignored_paralog_ids = [
            seq.id for seq in paralog_seqs if seq.id != selected_paralog.id
        ]

        logger.info(
            "%s seqs from %s were detected as paralogs. %s kept in ortholog group",
            (", ").join([seq_id for seq_id in ignored_paralog_ids]),
            genome_id,
            selected_paralog.id,
        )

        # remove paralogs from og
        og.remove_seqs(*ignored_paralog_ids)
        assert og.path.exists()

        # add ignored paralogs to Runtime
        paralog_ids.update(ignored_paralog_ids)

    return False


def _get_genome_ids_with_paralogs(genome_ids: list[str]) -> list[str]:
    """
    Find genome IDs that occur more than once.

    Returns:
        Sorted list of genome IDs that have multiple sequences.
    """
    assert len(set(genome_ids)) > 1
    counts = Counter(genome_ids)

    return sorted(genome_id for genome_id, count in counts.items() if count > 1)
