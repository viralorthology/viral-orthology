import shutil
from pathlib import Path


class File:
    """
    Represents a file and provides common file operations.

    An instance of File can be created without the file existing.

    Attributes:
        path: Path of a file
    """

    def __init__(self, path: Path) -> None:
        self.path = path

    @property
    def exists(self) -> bool:
        return self.path.is_file()

    @property
    def has_content(self) -> bool:
        return self.exists and self.path.stat().st_size > 0

    def read_file(self) -> str:
        """
        Get the contents of an existing file.

        Raises:
            FileNotFoundError: if the file does not exist or is empty
        """
        if not self.has_content:
            raise FileNotFoundError(f"{self.path} does not exist or is empty")

        return self.path.read_text(encoding="utf-8")

    def write_to_file(self, text: str) -> None:
        """
        Append text to the file.

        Args:
            text: text to append to the file
        """
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(text)

    def rename_file(self, new_filename: str) -> None:
        """
        Rename the file.

        Args:
            new_filename: new name for the file.

        Raises:
            FileNotFoundError: if the file does not exist
            FileExistsError: if a file with the new name already exists
        """
        new_path = self.path.with_name(new_filename)

        if new_path.exists():
            raise FileExistsError(f"{new_path} already exists")

        self.path = self.path.rename(new_path)

    def move_file(self, directory_path: Path) -> None:
        """
        Move the file from its current directory to the given directory.

        Args:
            directory_path: path to the directory where the file will be moved

        Raises:
            FileNotFoundError: if the file does not exist or the given directory does not exist
            FileExistsError: if a file with the same name already exists in the destination directory
        """
        new_path = directory_path / self.path.name

        if new_path.exists():
            raise FileExistsError(f"{new_path} already exists")

        shutil.move(self.path, new_path)
        self.path = new_path

    def delete_file(self) -> None:
        """
        Delete the file.

        Raises:
            FileNotFoundError: if the file does not exist
        """
        self.path.unlink()
