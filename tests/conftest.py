from collections.abc import Iterable

import pytest

from cli.args import Args
from cli.ui import UI
from config.context import Context
from config.paths import Paths


@pytest.fixture
def ui_always_yes():
    class UIAlwaysYes(UI):
        def show(self, text: str) -> None:
            pass

        def show_error(self, text: str) -> None:
            pass

        def progress_bar[T](self, iterable: Iterable[T]) -> Iterable[T]:
            return iterable

        def ask_yes_no(self, question: str) -> bool:
            return True

    return UIAlwaysYes()


@pytest.fixture
def context_always_yes(tmp_path, ui_always_yes):
    return Context(
        Args(debug=False, assume_yes=False, tool_args={}),
        Paths(tmp_path),
        ui_always_yes,
    )
