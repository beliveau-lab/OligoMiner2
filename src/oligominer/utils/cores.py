"""
# Core-count resolution

Resolve how many workers to use, correctly, when running under a batch scheduler.

Every parallel entry point in OligoMiner2 takes a ``cores`` or ``threads`` argument, and
those arguments used to default to ``1``. Passing a real value instead was measured to take
genome-wide mining from 4,649 s to 1,953 s (2.4x) with bit-identical results -- nothing had
to be rewritten, a default had to change. This module is that default.

## Why ``os.cpu_count()`` is the wrong answer

``cpu_count()`` reports the *machine's* cores. Under SGE, Slurm, or any cgroup-confined
shell that is not what the job was granted, and using it oversubscribes: N processes
fighting over the slots actually held, which is slower than running serially and inflates
every timing taken beside it. On one cluster login node ``cpu_count()`` reports 192 where
the real grant is 2 -- a 96x overshoot.

## Why affinity is only the fallback

``len(os.sched_getaffinity(0))`` reports the cores this process may run on, which is closer.
But it can still be wider than the slot grant: measured on an 8-slot SGE job, affinity
reported 16 because the cpuset was wider than the grant. Trusting affinity there runs 16
threads in 8 slots -- the same oversubscription one level down.

So the resolution order is: an explicit argument, then an environment override, then the
scheduler's own variable, and only then affinity.
"""

import os

# the scheduler's own statement of what this job was granted, in preference order
_GRANT_VARS = ("NSLOTS",                 # SGE / Altair Grid Engine
               "SLURM_CPUS_PER_TASK",    # Slurm
               "SLURM_CPUS_ON_NODE",
               "PBS_NP",                 # PBS / Torque
               "LSB_DJOB_NUMPROC")       # LSF

# environment variable a job script can set to fix the width without touching code
ENV_VAR = "OLIGOMINER_THREADS"


def _scheduler_grant():
    """
    Return the core count the batch scheduler granted this job.

    Returns:
        grant (int or None): the granted core count, or None if not running
            under a recognized scheduler.
    """
    grant = None
    for var in _GRANT_VARS:
        raw = os.environ.get(var)
        if raw and raw.strip().isdigit() and int(raw) >= 1:
            grant = int(raw)
            break

    # success
    return grant


def _affinity():
    """
    Return how many cores this process is allowed to run on.

    Returns:
        n (int): the affinity mask size, falling back to the machine core
            count on platforms without scheduler affinity.
    """
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:              # not Linux
        n = os.cpu_count() or 1

    # success
    return n


def resolve_cores(n=None, env_var=ENV_VAR, cap=None):
    """
    Return the number of workers to use.

    Resolution order: an explicit ``n`` always wins, then ``env_var``, then the
    batch scheduler's granted core count, then this process's CPU affinity.

    Args:
        n (int, optional): an explicit request, returned as-is when >= 1 so a
            caller that knows better always wins. None means "decide for me".
        env_var (str): environment variable consulted when n is None, so a job
            script can fix the width without touching code.
        cap (int, optional): upper bound, for stages whose speedup plateaus.

    Returns:
        cores (int): the number of workers to use, always at least 1.
    """
    if n is not None and int(n) >= 1:
        cores = int(n)
    else:
        raw = os.environ.get(env_var)
        if raw and raw.strip().isdigit() and int(raw) >= 1:
            cores = int(raw)
        else:
            cores = _scheduler_grant() or _affinity()

    if cap is not None:
        cores = min(cores, int(cap))

    # success
    return max(1, cores)


def describe():
    """
    Report every core count this process can see, for a provenance record.

    Useful in a run manifest: it records not just the value used but the values
    that were rejected, so an oversubscribed run can be diagnosed afterwards.

    Returns:
        info (dict): the scheduler grant and the variable it came from, the
            affinity mask size, the machine core count, any environment
            override, and the resolved value.
    """
    # success
    return {
        "scheduler_grant": _scheduler_grant(),
        "grant_var": next((v for v in _GRANT_VARS if os.environ.get(v)), None),
        "sched_getaffinity": _affinity(),
        "cpu_count": os.cpu_count(),
        "env_override": os.environ.get(ENV_VAR),
        "resolved": resolve_cores(),
    }
