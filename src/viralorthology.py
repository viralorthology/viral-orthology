import logging
import sys
from pathlib import Path

import utils
from cli.args import get_args
from cli.ui import CLI, UI, DebugCLI
from config.context import Context
from config.paths import Paths
from config.runtime import Runtime
from modules import (
    aa_composition,
    blastp_search,
    download_seqs,
    hmm_search,
    kimura,
    pipeline,
    protein_domain_search,
    secondary_structure,
    synteny,
    tertiary_structure,
    zscore,
)

logger = logging.getLogger(__name__)


def get_selected_module_flag() -> str:
    """
    Get the selected module flag.

    Exit if none or more than one module was selected
    """
    selected_modules = [flag for flag in MODULES if flag in sys.argv]

    if not selected_modules:
        sys.exit("No modules were selected. For help run viralorthology -help")

    if len(selected_modules) > 1:
        sys.exit("Two or more modules were selected. For help run viralorthology -help")

    return selected_modules[0]


def print_help(ui: UI) -> None:
    # TODO
    ui.show("help")


MODULES = {
    "-pipeline": pipeline,
    "-download-seqs": download_seqs,
    # enrichment modules
    "-composition": aa_composition,
    "-blastp": blastp_search,
    "-hmm": hmm_search,
    "-secondary-structure": secondary_structure,
    "-tertiary-structure": tertiary_structure,
    "-synteny": synteny,
    # analysis modules
    "-protein-domain": protein_domain_search,
    "-zscore": zscore,
    "-kimura": kimura,
}


def main() -> None:
    ui = DebugCLI() if "-debug" in sys.argv else CLI()

    if "-help" in sys.argv or "-h" in sys.argv:
        print_help(ui)
        return

    selected_module_flag = get_selected_module_flag()
    module = MODULES[selected_module_flag]
    paths = Paths(Path.cwd())

    try:
        args = get_args(selected_module_flag, sys.argv[1:])
        utils.configure_logging(paths.base, args.debug)
        ctx = Context(args, paths, ui, Runtime())
        module.run(ctx)

    except Exception as e:
        if "-debug" in sys.argv:
            logger.exception("EXCEPTION")
        else:
            logger.error("ERROR: %s", e)
            ui.show_error(str(e))  # TODO can show error two times when debug=False

        sys.exit(1)


if __name__ == "__main__":
    main()
