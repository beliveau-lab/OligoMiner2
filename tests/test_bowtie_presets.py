"""
Tests for the Bowtie2 parameter presets.

The presets are transcriptions of Bowtie2's own named modes. A transcription
error changes alignment sensitivity without any error, and would show up only
as a different set of off-target hits, so the values are checked against the
installed Bowtie2 rather than against a second copy of the same numbers.
"""

import re
import shutil
import subprocess

import pytest

from oligominer.specificity.alignment import bowtie_presets as presets

# preset name as Bowtie2 spells it -> the module constant holding it
PRESET_NAMES = {
    'very-fast': 'VERY_FAST',
    'fast': 'FAST',
    'sensitive': 'SENSITIVE',
    'very-sensitive': 'VERY_SENSITIVE',
    'very-fast-local': 'VERY_FAST_LOCAL',
    'fast-local': 'FAST_LOCAL',
    'sensitive-local': 'SENSITIVE_LOCAL',
    'very-sensitive-local': 'VERY_SENSITIVE_LOCAL',
}

needs_bowtie2 = pytest.mark.skipif(shutil.which('bowtie2') is None,
                                   reason='bowtie2 is not installed')


def bowtie2_presets():
    """
    Parse Bowtie2's own preset definitions out of its help text.

    Returns:
        parsed (dict): preset name mapped to its D, R, N, L and i values.
    """
    out = subprocess.run(['bowtie2', '--help'], capture_output=True, text=True)
    text = out.stdout + out.stderr

    pattern = re.compile(
        r'--([a-z-]+)\s+-D\s+(\d+)\s+-R\s+(\d+)\s+-N\s+(\d+)\s+-L\s+(\d+)'
        r'\s+-i\s+(\S+)')

    parsed = {}
    for name, d, r, n, ell, i in pattern.findall(text):
        if name in PRESET_NAMES:
            parsed[name] = {'D': int(d), 'R': int(r), 'N': int(n),
                            'L': int(ell), 'i': i}

    # success
    return parsed


class TestPresetsMatchBowtie2:

    @needs_bowtie2
    def test_every_preset_is_found_in_the_help_text(self):
        found = bowtie2_presets()
        assert set(found) == set(PRESET_NAMES), 'bowtie2 help parse changed'

    @needs_bowtie2
    @pytest.mark.parametrize('name', sorted(PRESET_NAMES))
    def test_the_transcription_is_faithful(self, name):
        expected = bowtie2_presets()[name]
        preset = getattr(presets, PRESET_NAMES[name])

        for key, value in expected.items():
            assert preset[key] == value, f'{name}: {key}'

    @needs_bowtie2
    def test_local_presets_are_flagged_local(self):
        for name, constant in PRESET_NAMES.items():
            preset = getattr(presets, constant)
            assert preset['local'] is name.endswith('-local'), name


class TestPresetShape:
    """These hold without bowtie2 installed."""

    def test_every_preset_carries_the_same_keys(self):
        keys = {'D', 'R', 'N', 'L', 'i', 'local'}
        for constant in PRESET_NAMES.values():
            assert set(getattr(presets, constant)) == keys, constant

    def test_sensitivity_increases_the_seed_extension_effort(self):
        # D is the number of consecutive failed extensions before giving up
        assert (presets.VERY_FAST['D'] < presets.FAST['D']
                < presets.SENSITIVE['D'] < presets.VERY_SENSITIVE['D'])

    def test_the_index_extensions_are_the_six_bowtie2_writes(self):
        assert presets.BT2_INDEX_EXTENSIONS == [
            '.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
