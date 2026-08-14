"""
# Core-count resolution

Resolves how many workers a parallel stage should use when running under a batch
scheduler, and supplies the default for every ``cores`` and ``threads`` argument in the
package.

Resolution order:

1. an explicit argument passed by the caller
2. the ``OLIGOMINER_THREADS`` environment variable
3. the scheduler's granted core count (``NSLOTS``, ``SLURM_CPUS_PER_TASK``, and the
   equivalents for PBS and LSF)
4. this process's CPU affinity mask

The scheduler's own variable is preferred over the affinity mask because a cpuset can be
wider than the slot grant, and the affinity mask is preferred over ``os.cpu_count()``
because that reports the whole machine rather than the allocation. Using either of the
wider numbers oversubscribes the allocation, which runs slower than serial execution and
inflates any timing measured alongside it.
"""

import os

# environment variables by which each scheduler reports its core grant, in preference order
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

    Resolution order: an explicit ``n``, then ``env_var``, then the batch
    scheduler's granted core count, then this process's CPU affinity.

    Args:
        n (int, optional): an explicit request, returned as-is when >= 1. None
            resolves from the environment instead.
        env_var (str): environment variable consulted when n is None, letting a
            job script set the width without a code change.
        cap (int, optional): upper bound on the resolved value.

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

    Records the value used alongside the ones that were rejected, so an
    oversubscribed run can be diagnosed from its manifest.

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
