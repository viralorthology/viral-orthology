from dataclasses import dataclass


@dataclass(frozen=True)
class Args:
    debug: bool
    assume_yes: bool
    tool_args: dict[str, list[str]]


AVAILABLE_TOOLS_BY_PUBLIC_MODULE = {
    "pipeline": {
        "blastp",
        "proteinortho",
        "hmmsearch",
        "blastp_paralog_search",
        "orffinder",
    },
    "blastp_module": {"blastp"},
    "hmmsearch_module": {"hmmsearch"},
    "synteny": {"blastp"},
}

DEFAULT_TOOL_PARAMS: dict[str, dict[str, dict[str, str | None]]] = {
    "pipeline": {
        "blastp": {"-word_size": "2", "-evalue": "0.001", "-qcov_hsp_perc": "40"},
        "orffinder": {"-ml": "90", "-s": "0"},
        "blastp_paralog_search": {"-evalue": "0.00001"},
    },
    "blastp_module": {"blastp": {"-word_size": "2", "-evalue": "0.1"}},
    "hmmsearch_module": {
        "hmmsearch": {"--nobias": None}
    },  # None indicates a flag without an argument
}

GLOBAL_BOOL_FLAGS = {"assume_yes", "debug"}


def get_args(selected_public_module_flag: str, argv: list[str]) -> Args:
    """
    Parse command-line arguments and apply default parameters.

    Assumes exactly one module is selected to run.
    """
    assert argv.count(selected_public_module_flag) == 1
    argv = [arg for arg in argv if arg != selected_public_module_flag]
    assert "-" in selected_public_module_flag
    selected_public_module_flag = selected_public_module_flag.lstrip("-")
    assert selected_public_module_flag in AVAILABLE_TOOLS_BY_PUBLIC_MODULE

    # parse args
    global_bool_flags = {flag: False for flag in GLOBAL_BOOL_FLAGS}
    current_tool = None
    tool_args: dict[str, list[str]] = {}
    for arg in argv:
        clean_arg = arg.lstrip("-")

        if clean_arg in AVAILABLE_TOOLS_BY_PUBLIC_MODULE[selected_public_module_flag]:
            current_tool = clean_arg
            tool_args[clean_arg] = []
            continue
        if clean_arg in GLOBAL_BOOL_FLAGS:
            global_bool_flags[clean_arg] = True
            continue

        if current_tool is None:
            raise ValueError(f"Unexpected argument: {arg}")

        tool_args[current_tool].append(arg)

    # TODO add arg and param validation

    # apply default params
    module_default_params = DEFAULT_TOOL_PARAMS.get(selected_public_module_flag, {})
    for tool, default_params in module_default_params.items():
        tool_params = tool_args.setdefault(tool, [])

        for param, value in default_params.items():
            if param not in tool_params:
                if value is None:
                    tool_params.append(param)
                else:
                    tool_params.extend([param, value])

    return Args(**global_bool_flags, tool_args=tool_args)
