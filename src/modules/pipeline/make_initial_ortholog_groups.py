import shutil
import tempfile
from pathlib import Path

import utils
from config.context import Context
from models.fasta import Fasta


def make_initial_ortholog_groups(ctx: Context) -> None:
    """
    Create initial ortholog groups from the non-redundant proteomes.

    Temporary files are used to run ProteinOrtho, and the resulting
    ortholog group FASTA files are copied to the ortholog groups directory.
    """
    ctx.ui.show("Creating initial ortholog groups...")

    proteome_fastas = utils.get_fastas(ctx.paths.sequences_dir, ".proteome")
    non_redundant_proteome_fastas = [
        fasta
        for fasta in proteome_fastas
        if fasta.path.stem not in ctx.runtime.redundant_genome_ids
    ]

    with tempfile.TemporaryDirectory() as tmp_path:
        tmp_dir = Path(tmp_path)
        tmp_proteomes = []
        for fasta in non_redundant_proteome_fastas:
            shutil.copy(fasta.path, tmp_dir / fasta.path.name)
            tmp_proteomes.append(Fasta(tmp_dir / fasta.path.name))

        _run_proteinortho(tmp_proteomes, ctx.args.tool_args["proteinortho"], tmp_dir)

        ctx.paths.ortholog_groups_dir.mkdir()
        ortholog_groups = utils.get_fastas(tmp_dir, ".fasta")
        for og in ortholog_groups:
            shutil.copy(og.path, ctx.paths.ortholog_groups_dir / og.path.name)


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
