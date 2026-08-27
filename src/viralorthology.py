import sys
from pathlib import Path

from cli.args import get_args
from cli.ui import CLI, UI
from config.context import Context
from config.paths import Paths
from modules import pipeline


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


MODULES = {"-pipeline": pipeline}


def main() -> None:
    ui = CLI()

    if "-help" in sys.argv:
        print_help(ui)
        return

    selected_module_flag = get_selected_module_flag()
    module = MODULES[selected_module_flag]
    paths = Paths(Path.cwd())

    try:
        args = get_args(selected_module_flag, sys.argv[1:])
        ctx = Context(args, paths, ui)
        module.validate(ctx)
        module.run(ctx)

    except Exception as e:
        if "-debug" in sys.argv:
            raise
        ui.show_error(str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
