import pytest

from cli.args import get_args


def test_unexpected_argument():
    with pytest.raises(ValueError):
        get_args("-pipeline", ["-pipeline", "test"])


def test_module_no_hyphen():
    with pytest.raises(AssertionError):
        get_args("pipeline", ["-pipeline"])


def test_default_params_no_args():
    args = get_args("-pipeline", ["-pipeline"])
    assert args.tool_args["blastp"] == [
        "-word_size",
        "2",
        "-evalue",
        "0.001",
        "-qcov_hsp_perc",
        "40",
    ]


def test_default_params_with_args():
    args = get_args("-pipeline", ["-pipeline", "-blastp", "-word_size", "10"])
    assert args.tool_args["blastp"] == [
        "-word_size",
        "10",
        "-evalue",
        "0.001",
        "-qcov_hsp_perc",
        "40",
    ]


def test_global_flags_default_to_false():
    args = get_args("-pipeline", ["-pipeline"])

    assert args.debug is False
    assert args.assume_yes is False


def test_debug_flag():
    args = get_args("-pipeline", ["-pipeline", "--debug"])

    assert args.debug is True


def test_assume_yes_flag():
    args = get_args("-pipeline", ["-pipeline", "--assume_yes"])

    assert args.assume_yes is True


def test_both_global_flags():
    args = get_args(
        "-pipeline",
        ["-pipeline", "-debug", "-assume_yes"],
    )

    assert args.debug is True
    assert args.assume_yes is True


def test_multiple_tools():
    args = get_args(
        "-pipeline",
        [
            "-pipeline",
            "-blastp",
            "-evalue",
            "10",
            "-hmmsearch",
            "--nobias",
        ],
    )

    assert args.tool_args["blastp"] == [
        "-evalue",
        "10",
        "-word_size",
        "2",
        "-qcov_hsp_perc",
        "40",
    ]
    assert args.tool_args["hmmsearch"] == ["--nobias"]
