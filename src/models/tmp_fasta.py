import tempfile
from pathlib import Path
from types import TracebackType
from typing import Literal

from typing_extensions import Self

from models.fasta import Fasta


class TmpFasta(Fasta):
    """
    A temporary FASTA file that is automatically cleaned up on exit.

    The FASTA file is created inside a temporary directory and can be used
    as a context manager. The temporary directory and its contents are
    removed when leaving the context manager.
    """

    def __init__(self) -> None:
        self.tmp_dir = tempfile.TemporaryDirectory()
        super().__init__(Path(self.tmp_dir.name) / "tmp.fasta")

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> Literal[False]:
        self.tmp_dir.cleanup()
        return False
