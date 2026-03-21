"""Tests for the utils subpackage.

Covers: seq_utils (rev_comp, calc_gc), file_paths, required_files,
and input_dispatch.
"""

import os
import tempfile

import pytest

from oligominer.utils.seq_utils import rev_comp, calc_gc
from oligominer.utils.file_paths import get_abs_path, get_dir_name
from oligominer.utils.required_files import (
    check_input_exists,
    check_output_exists,
    check_dir_exists,
)
from oligominer.utils.input_dispatch import require_one_of
from oligominer.utils.exceptions import (
    MissingInputFile,
    MissingOutputFile,
    MissingDirectory,
    InvalidInputError,
)


# ---------------------------------------------------------------------------
# seq_utils
# ---------------------------------------------------------------------------

class TestRevComp:

    def test_basic(self):
        assert rev_comp("ATCG") == "CGAT"

    def test_palindrome(self):
        assert rev_comp("AATT") == "AATT"

    def test_single_base(self):
        assert rev_comp("A") == "T"
        assert rev_comp("C") == "G"

    def test_lowercase(self):
        assert rev_comp("atcg") == "cgat"

    def test_empty(self):
        assert rev_comp("") == ""

    def test_round_trip(self):
        seq = "ACGTACGT"
        assert rev_comp(rev_comp(seq)) == seq


class TestCalcGc:

    def test_all_gc(self):
        assert calc_gc("GCGC") == pytest.approx(1.0)

    def test_all_at(self):
        assert calc_gc("ATAT") == pytest.approx(0.0)

    def test_half_gc(self):
        assert calc_gc("ATGC") == pytest.approx(0.5)

    def test_as_percent(self):
        assert calc_gc("ATGC", as_percent=True) == pytest.approx(50.0)

    def test_lowercase(self):
        assert calc_gc("atgc") == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# file_paths
# ---------------------------------------------------------------------------

class TestFilePaths:

    def test_get_abs_path(self):
        result = get_abs_path("test.txt")
        assert os.path.isabs(result)

    def test_get_dir_name(self, tmp_path):
        path = str(tmp_path / "subdir" / "file.txt")
        dirname = get_dir_name(path)
        assert dirname.endswith("subdir")


# ---------------------------------------------------------------------------
# required_files
# ---------------------------------------------------------------------------

class TestRequiredFiles:

    def test_check_input_exists(self, example_fasta_path):
        assert check_input_exists(example_fasta_path) is True

    def test_check_input_missing(self):
        with pytest.raises(MissingInputFile):
            check_input_exists("/nonexistent/file.txt")

    def test_check_output_missing(self):
        with pytest.raises(MissingOutputFile):
            check_output_exists("/nonexistent/file.txt")

    def test_check_dir_exists(self, tmp_path):
        assert check_dir_exists(str(tmp_path)) is True

    def test_check_dir_missing(self):
        with pytest.raises(MissingDirectory):
            check_dir_exists("/nonexistent/directory")

    def test_check_dir_create(self, tmp_path):
        new_dir = str(tmp_path / "new_subdir")
        assert not os.path.exists(new_dir)
        check_dir_exists(new_dir, create=True)
        assert os.path.exists(new_dir)

    def test_check_dir_parent(self, tmp_path):
        file_path = str(tmp_path / "existing" / "file.txt")
        check_dir_exists(file_path, parent_dir=True, create=True)
        assert os.path.exists(str(tmp_path / "existing"))


# ---------------------------------------------------------------------------
# input_dispatch
# ---------------------------------------------------------------------------

class TestInputDispatch:

    def test_require_one_of_first(self):
        # should not raise when exactly one is provided
        require_one_of("hello", None, "a", "b")

    def test_require_one_of_second(self):
        require_one_of(None, "world", "a", "b")

    def test_require_one_of_neither(self):
        with pytest.raises(InvalidInputError):
            require_one_of(None, None, "a", "b")

    def test_require_one_of_both(self):
        with pytest.raises(InvalidInputError):
            require_one_of("hello", "world", "a", "b")
