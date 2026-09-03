import utils
from engines.muscle import muscle_align
from models.ortholog_group import OrthologGroup


def hmm_compare_groups(
    group1_fasta: OrthologGroup, group2_fasta: OrthologGroup
) -> float:
    """
    Compare two ortholog groups using HHsuite profile-profile alignment.

    Each ortholog group is aligned with MUSCLE if an alignment does not
    already exist. The alignments are then converted to A2M format and
    used to build HHsuite HMM profiles. The two profiles are finally
    compared with HHalign.

    Returns:
        The HHalign probability score between the two HMM profiles.
    """
    if not group1_fasta.alignment_path.is_file():
        muscle_align(group1_fasta)
    if not group2_fasta.alignment_path.is_file():
        muscle_align(group2_fasta)

    _convert_aligned_fasta_to_a2m(group1_fasta)
    _convert_aligned_fasta_to_a2m(group2_fasta)
    _build_hmm_hhsuite(group1_fasta)
    _build_hmm_hhsuite(group2_fasta)

    score = _hhalign(group1_fasta, group2_fasta)

    return score


def _convert_aligned_fasta_to_a2m(fasta: OrthologGroup) -> None:
    """
    Convert a FASTA alignment to A2M format using HHsuite.
    """
    utils.run_cmd(f"reformat_msa fas a2m {fasta.alignment_path} {fasta.a2m_path}")


def _build_hmm_hhsuite(fasta: OrthologGroup) -> None:
    """
    Build an HHsuite HMM profile from an A2M alignment.

    Args:
        fasta: Ortholog group containing the alignment used to build the HMM.
    """
    utils.run_cmd(f"hhmake -i {fasta.a2m_path} -o {fasta.hmm_hhsuite_path}")


def _hhalign(o_group1: OrthologGroup, o_group2: OrthologGroup) -> float:
    """
    Compare two HHsuite HMM profiles using HHalign.

    Returns:
        The HHalign probability score between the two HMM profiles.
    """
    output = utils.run_cmd(
        f"hhalign -i {o_group1.hmm_hhsuite_path} -t {o_group2.hmm_hhsuite_path} -o stdout -glob"
    )

    score = 0.0
    for line in output.splitlines():
        line = line.strip()
        if line.startswith("Probab="):
            score = float(line.split()[0].split("=")[1])
            break

    return score
