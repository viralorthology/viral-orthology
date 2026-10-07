from pathlib import Path


def ensure_files_exist(*file_paths: Path) -> None:
    """
    Ensure that all given files exist.

    Args:
        *file_paths: Paths to files that must exist.

    Raises:
        ValueError: If no file path is provided.
        FileNotFoundError: If any of the files do not exist.
    """
    if not file_paths:
        raise ValueError("At least one file path must be provided.")

    for file_path in file_paths:
        if not file_path.is_file():
            raise FileNotFoundError(f"{file_path} does not exist")


def ensure_paths_do_not_exist(*paths: Path) -> None:
    """
    Ensure that none of the given files exist.

    Args:
        *file_paths: Paths to files that must not exist.

    Raises:
        ValueError: If no file paths are provided.
        FileExistsError: If any of the files already exist.
    """
    if not paths:
        raise ValueError("At least one file path must be provided")

    for path in paths:
        if path.exists():
            raise FileExistsError(f"{path} already exists")


def ensure_files_have_content(*file_paths: Path) -> None:
    """
    Ensure that all given files exist and have content.

    Args:
        *file_paths: Paths to files that must exist and contain content.

    Raises:
        FileNotFoundError: If no file path is provided or any of the files do not exist
        ValueError: If any of the files is empty
    """
    ensure_files_exist(*file_paths)
    for file_path in file_paths:
        if not file_path.stat().st_size > 0:
            raise ValueError(f"{file_path} is empty")
