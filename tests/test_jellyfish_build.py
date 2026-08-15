"""
Tests for sizing and building a Jellyfish index.

Jellyfish sizes its hash by the number of distinct k-mers it expects. A fixed
default makes a small genome allocate for a large one, which is what
om2_issues #20 records.
"""

import os
import shutil

import pytest

from oligominer.specificity.kmers.jellyfish_build import (
    MIN_HASH_SIZE, estimate_hash_size, jellyfish_build,
)


@pytest.fixture
def small_genome(tmp_path):
    # above the floor, so the estimate is the file's own size
    path = tmp_path / 'small.fa'
    path.write_text('>chr1\n' + 'ACGTACGTAC' * 150_000 + '\n')
    return str(path)


@pytest.fixture
def tiny_genome(tmp_path):
    path = tmp_path / 'tiny.fa'
    path.write_text('>chr1\nACGTACGTACGTACGTACGTACGT\n')
    return str(path)


class TestEstimateHashSize:

    def test_it_scales_with_the_input(self, small_genome, tiny_genome):
        assert estimate_hash_size(small_genome) > estimate_hash_size(
            tiny_genome)

    def test_a_tiny_input_still_gets_a_usable_hash(self, tiny_genome):
        # a 24-base genome must not ask jellyfish for a 24-slot table
        assert estimate_hash_size(tiny_genome) == MIN_HASH_SIZE

    def test_a_large_input_is_sized_from_its_bytes(self, small_genome):
        assert estimate_hash_size(small_genome) == os.path.getsize(
            small_genome)

    def test_the_floor_is_configurable(self, tiny_genome):
        assert estimate_hash_size(tiny_genome, floor=50) == 50


@pytest.mark.skipif(shutil.which('jellyfish') is None,
                    reason='jellyfish is not installed')
class TestBuild:

    def test_an_index_is_written(self, small_genome, tmp_path):
        out = tmp_path / 'index.jf'
        path = jellyfish_build(small_genome, str(out), k=12)

        assert os.path.exists(path)
        assert os.path.getsize(path) > 0

    def test_it_builds_without_being_told_a_size(self, tiny_genome, tmp_path):
        # the whole point of the estimate: a caller need not know the genome
        out = tmp_path / 'tiny.jf'
        assert os.path.exists(jellyfish_build(tiny_genome, str(out), k=12))

    def test_an_explicit_size_is_still_honoured(self, small_genome, tmp_path):
        out = tmp_path / 'explicit.jf'
        assert os.path.exists(
            jellyfish_build(small_genome, str(out), k=12, size='10M'))
