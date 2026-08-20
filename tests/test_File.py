import pytest


@pytest.mark.parametrize(
    "file_fixture,expected",
    [("nonexistent_file", False), ("empty_file", True), ("text_file", True)],
)
def test_exists(request, file_fixture, expected):
    file = request.getfixturevalue(file_fixture)
    assert file.exists is expected


@pytest.mark.parametrize(
    "file_fixture,expected",
    [("nonexistent_file", False), ("empty_file", False), ("text_file", True)],
)
def test_has_content(request, file_fixture, expected):
    file = request.getfixturevalue(file_fixture)
    assert file.has_content is expected


def test_move_nonexistent_file(nonexistent_file, tmp_path):
    with pytest.raises(FileNotFoundError):
        nonexistent_file.move_file(tmp_path)


def test_move_file(text_file, tmp_path):
    dir_to = tmp_path / "dir_to"
    dir_to.mkdir()
    original_path = text_file.path
    text_file.move_file(dir_to)
    assert text_file.path == dir_to / original_path.name
    assert text_file.path.exists()
    assert not original_path.exists()


def test_move_file_to_nonexistent_dir(text_file, tmp_path):
    dir_to = tmp_path / "dir_to"
    with pytest.raises(FileNotFoundError):
        text_file.move_file(dir_to)


def test_move_file_to_existing_file_path(text_file, tmp_path):
    with pytest.raises(FileExistsError):
        text_file.move_file(tmp_path)


def test_rename_nonexistent_file(nonexistent_file):
    with pytest.raises(FileNotFoundError):
        nonexistent_file.rename_file("testfile.txt")


def test_rename_file(text_file):
    text_file.rename_file("testfile.txt")
    assert text_file.path.name == "testfile.txt"


def test_rename_file_to_existent_file_path(empty_file, text_file):
    other_file_name = text_file.path.name
    with pytest.raises(FileExistsError):
        empty_file.rename_file(other_file_name)


def test_delete_nonexistent_file(nonexistent_file):
    with pytest.raises(FileNotFoundError):
        nonexistent_file.delete_file()


def test_delete_file(text_file):
    text_file.delete_file()
    assert text_file.exists is False


def test_read_nonexistent_file(nonexistent_file):
    with pytest.raises(FileNotFoundError):
        nonexistent_file.read_file()


def test_read_empty_file(empty_file):
    with pytest.raises(FileNotFoundError):
        empty_file.read_file()


def test_read_file(text_file):
    file_content = text_file.read_file()
    assert file_content == "sometext"


def test_write_to_nonexistent_file(nonexistent_file):
    assert nonexistent_file.exists is False
    nonexistent_file.write_to_file("test")
    assert nonexistent_file.exists is True
    assert nonexistent_file.has_content is True
    assert nonexistent_file.read_file() == "test"
