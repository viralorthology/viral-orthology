from abc import ABC, abstractmethod
from collections.abc import Iterable

from tqdm import tqdm  # type: ignore[import-untyped]


class UI(ABC):
    """Interface for interacting with the user."""

    @abstractmethod
    def show(self, text: str) -> None: ...

    @abstractmethod
    def progress_bar[T](self, iterable: Iterable[T]) -> Iterable[T]:
        """Provide progress feedback while iterating."""
        ...

    @abstractmethod
    def ask_yes_no(self, question: str) -> bool:
        """Prompt the user for a yes/no answer."""
        ...


class CLI(UI):
    def show(self, text: str) -> None:
        print(text)

    def progress_bar[T](self, iterable: Iterable[T]) -> Iterable[T]:
        return tqdm(iterable)  # type: ignore[no-any-return]

    def ask_yes_no(self, question: str) -> bool:
        """Prompt the user until a valid yes/no answer is provided."""
        while True:
            yes_no = input(f"{question} [y/n]: ").strip().lower()

            if yes_no in {"y", "n"}:
                return yes_no == "y"

            print("Wrong input, try again.")
