"""Tests for core-count resolution.

The failure this guards against is oversubscription: resolving to the machine's
core count, or to an affinity mask wider than the scheduler's grant, and then
running that many workers inside a smaller slot allocation.
"""

import os

import pytest

from oligominer.utils import cores


@pytest.fixture
def clean_env(monkeypatch):
    """Remove every scheduler and override variable so a test starts from nothing."""
    for var in cores._GRANT_VARS:
        monkeypatch.delenv(var, raising=False)
    monkeypatch.delenv(cores.ENV_VAR, raising=False)


class TestResolutionOrder:

    def test_explicit_wins_over_everything(self, clean_env, monkeypatch):
        monkeypatch.setenv("NSLOTS", "8")
        monkeypatch.setenv(cores.ENV_VAR, "4")
        assert cores.resolve_cores(2) == 2

    def test_env_var_wins_over_scheduler(self, clean_env, monkeypatch):
        monkeypatch.setenv("NSLOTS", "8")
        monkeypatch.setenv(cores.ENV_VAR, "4")
        assert cores.resolve_cores() == 4

    def test_scheduler_grant_wins_over_affinity(self, clean_env, monkeypatch):
        # the measured trap: an 8-slot SGE job whose affinity mask reports 16
        monkeypatch.setenv("NSLOTS", "8")
        monkeypatch.setattr(cores, "_affinity", lambda: 16)
        assert cores.resolve_cores() == 8

    def test_affinity_used_when_no_scheduler(self, clean_env, monkeypatch):
        monkeypatch.setattr(cores, "_affinity", lambda: 6)
        assert cores.resolve_cores() == 6

    def test_scheduler_var_preference_order(self, clean_env, monkeypatch):
        # NSLOTS precedes SLURM_CPUS_PER_TASK in _GRANT_VARS
        monkeypatch.setenv("NSLOTS", "8")
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
        assert cores.resolve_cores() == 8


class TestGuards:

    def test_never_returns_below_one(self, clean_env, monkeypatch):
        monkeypatch.setattr(cores, "_affinity", lambda: 0)
        assert cores.resolve_cores() == 1

    def test_cap_is_applied(self, clean_env, monkeypatch):
        monkeypatch.setenv("NSLOTS", "64")
        assert cores.resolve_cores(cap=8) == 8

    def test_cap_does_not_raise_a_smaller_grant(self, clean_env, monkeypatch):
        monkeypatch.setenv("NSLOTS", "2")
        assert cores.resolve_cores(cap=8) == 2

    def test_explicit_zero_falls_through_rather_than_returning_zero(self, clean_env, monkeypatch):
        # 0 is not a usable worker count, so it must not be taken as an explicit request
        monkeypatch.setenv("NSLOTS", "4")
        assert cores.resolve_cores(0) == 4

    @pytest.mark.parametrize("junk", ["", "  ", "many", "-4", "2.5"])
    def test_non_numeric_scheduler_values_are_ignored(self, clean_env, monkeypatch, junk):
        monkeypatch.setenv("NSLOTS", junk)
        monkeypatch.setattr(cores, "_affinity", lambda: 3)
        assert cores.resolve_cores() == 3


class TestDescribe:

    def test_describe_reports_the_rejected_values_too(self, clean_env, monkeypatch):
        monkeypatch.setenv("NSLOTS", "8")
        monkeypatch.setattr(cores, "_affinity", lambda: 16)

        info = cores.describe()

        assert info["resolved"] == 8
        assert info["scheduler_grant"] == 8
        assert info["grant_var"] == "NSLOTS"
        # the rejected wider value stays visible, so oversubscription is diagnosable
        assert info["sched_getaffinity"] == 16
        assert info["cpu_count"] == os.cpu_count()

    def test_describe_without_a_scheduler(self, clean_env, monkeypatch):
        monkeypatch.setattr(cores, "_affinity", lambda: 4)

        info = cores.describe()

        assert info["scheduler_grant"] is None
        assert info["grant_var"] is None
        assert info["resolved"] == 4
