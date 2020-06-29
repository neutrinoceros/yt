from yt.data_objects.time_series import get_filenames_from_glob_pattern

from pathlib import Path
import os
import tempfile
from yt.testing import assert_raises

def test_pattern_expansion():
    file_list = ["fake_data_file_{}".format(str(i).zfill(4)) for i in range(10)]

    with tempfile.TemporaryDirectory() as tmpdir:
        for file in file_list:
            (Path(tmpdir) / file).touch()

        pattern = os.path.join(tmpdir, "fake_data_file_*")
        found = get_filenames_from_glob_pattern(pattern)
        assert found == [os.path.join(tmpdir, file) for file in file_list]

        found2 = get_filenames_from_glob_pattern(Path(pattern))
        assert found2 == found

# future: equivalent implementation using pytest fixtures
"""
import pytest
def test_pattern_expansion(tmpdir):
    file_list = ["fake_data_file_{}".format(str(i).zfill(4)) for i in range(10)]
    for file in file_list:
        tmpdir.join(file).write("")
    pattern = tmpdir.join("fake_data_file_*")
    found = get_filenames_from_glob_pattern(pattern)
    assert found == [tmpdir.join(file) for file in file_list]

    found2 = get_filenames_from_glob_pattern(Path(pattern))
    assert found2 == found
"""

def test_no_match_pattern():
    with tempfile.TemporaryDirectory() as tmpdir:
        pattern = os.path.join(tmpdir, "fake_data_file_*")
        assert_raises(IOError, get_filenames_from_glob_pattern, pattern)

# future: equivalent implementation using pytest fixtures
"""
def test_no_match_pattern(tmpdir):
    pattern = tmpdir.join("fake_data_file_*")
    with pytest.raises(IOError):
        get_filenames_from_glob_pattern(pattern)
"""
