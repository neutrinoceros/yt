from yt.data_objects.time_series import get_filenames_from_glob_pattern
import pytest
from pathlib import Path

def test_pattern_expansion(tmpdir):
    file_list = ["fake_data_file_{}".format(str(i).zfill(4)) for i in range(10)]
    for file in file_list:
        tmpdir.join(file).write("")
    pattern = tmpdir.join("fake_data_file_*")
    found = get_filenames_from_glob_pattern(pattern)
    assert found == [tmpdir.join(file) for file in file_list]

    found2 = get_filenames_from_glob_pattern(Path(pattern))
    assert found2 == found

def test_no_match_pattern(tmpdir):
    pattern = tmpdir.join("fake_data_file_*")
    with pytest.raises(FileNotFoundError):
        get_filenames_from_glob_pattern(pattern)
