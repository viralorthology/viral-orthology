import pytest

from classes.file import File


@pytest.fixture
def nonexistent_file(tmp_path):
    path = tmp_path / "this_file_does_not_exist.txt"
    return File(path)


@pytest.fixture
def empty_file(tmp_path):
    path = tmp_path / "empty_file.txt"
    path.write_text("", encoding="utf-8")
    return File(path)


@pytest.fixture
def text_file(tmp_path):
    path = tmp_path / "text_file.txt"
    path.write_text("sometext", encoding="utf-8")
    return File(path)
