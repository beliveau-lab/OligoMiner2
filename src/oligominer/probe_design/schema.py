"""
# Probe set schema and provenance

Defines the on-disk format for a probe set: a version, the probe table, and the
record of how the probes were produced.

A probe CSV carries sequences and coordinates but nothing about the reference
they were designed against, the parameters used, the model that scored them, or
how many candidates were discarded on the way. Without that a probe set cannot be
reproduced, compared against another, or checked for whether it answers the
question that was asked.

## Versioning

`SCHEMA_VERSION` is the version this package writes. `read_manifest` accepts any
version this package knows how to read and upgrades older documents through
`MIGRATIONS`, so a probe set written by an earlier release stays loadable.

An unknown future version raises rather than being read optimistically, because a
document written by a newer release may carry columns whose meaning this release
would guess at.
"""

import json
from datetime import datetime, timezone

# the schema version this package writes
SCHEMA_VERSION = 1

# versions this package can read, oldest first
SUPPORTED_VERSIONS = (1,)


class SchemaError(Exception):
    """Raised when a probe set document cannot be read at this version."""


def utc_now():
    """
    Return the current UTC time as an ISO 8601 string.

    Returns:
        stamp (str): the timestamp, second resolution, with a Z suffix.
    """
    # success
    return datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def new_manifest(name=None):
    """
    Build an empty manifest at the current schema version.

    Args:
        name (str, optional): a human-readable name for the probe set.

    Returns:
        manifest (dict): the manifest.
    """
    # success
    return {
        'schema_version': SCHEMA_VERSION,
        'name': name,
        'created': utc_now(),
        'oligominer_version': _package_version(),
        'n_probes': 0,
        'columns': [],
        'target': {},
        'stages': [],
        'attrition': {},
    }


def _package_version():
    """
    Return the installed oligominer version.

    Returns:
        version (str): the version string, or 'unknown' if unavailable.
    """
    try:
        from oligominer._version import __version__
    except ImportError:
        return 'unknown'

    # success
    return __version__


def record_stage(manifest, name, params=None, n_in=None, n_out=None,
                 **details):
    """
    Append a pipeline stage to a manifest.

    Each stage records what ran, the parameters it ran under, and how many probes
    entered and left, which is what makes a short probe set explainable rather
    than merely small.

    A generative stage such as mining records n_out alone, since no probes
    entered it. A filtering stage records both.

    Args:
        manifest (dict): the manifest to append to.
        name (str): the stage name, e.g. 'mine', 'kmer_screen', 'align'.
        params (dict, optional): the parameters the stage ran under.
        n_in (int, optional): probes entering the stage. Omit for a stage that
            creates probes rather than filtering them.
        n_out (int, optional): probes leaving the stage.
        **details: any further fields to record with the stage.

    Returns:
        manifest (dict): the same manifest, for chaining.

    Raises:
        SchemaError: if a stage reports more probes leaving than entered.
    """
    if n_in is not None and n_out is not None and int(n_out) > int(n_in):
        raise SchemaError(
            f'stage {name!r} reports {n_out} probes leaving and {n_in} entering; '
            f'a filtering stage cannot emit more than it received. Omit n_in for '
            f'a stage that creates probes.')

    stage = {
        'stage': name,
        'at': utc_now(),
        'params': params or {},
    }
    if n_in is not None:
        stage['n_in'] = int(n_in)
    if n_out is not None:
        stage['n_out'] = int(n_out)
    if n_in is not None and n_out is not None:
        stage['n_dropped'] = int(n_in) - int(n_out)
    stage.update(details)

    manifest['stages'].append(stage)
    manifest['attrition'] = attrition_summary(manifest)

    # success
    return manifest


def attrition_summary(manifest):
    """
    Summarize where probes were lost across the recorded stages.

    Args:
        manifest (dict): a manifest carrying stages.

    Returns:
        summary (dict): n_in at the first stage, n_out at the last, and the
            number dropped per stage that reported both.
    """
    stages = manifest.get('stages', [])
    with_out = [s for s in stages if 'n_out' in s]
    if not with_out:
        return {}

    filtering = [s for s in stages if 'n_dropped' in s]

    # success
    return {
        'n_in': with_out[0]['n_out'],
        'n_out': with_out[-1]['n_out'],
        'by_stage': {s['stage']: s['n_dropped'] for s in filtering},
    }


def validate_manifest(manifest):
    """
    Check a manifest carries the fields this schema requires.

    Args:
        manifest (dict): the manifest to check.

    Returns:
        manifest (dict): the manifest, unchanged.

    Raises:
        SchemaError: if a required field is missing or the version is unreadable.
    """
    if not isinstance(manifest, dict):
        raise SchemaError(f'manifest must be a dict, got {type(manifest).__name__}')

    version = manifest.get('schema_version')
    if version is None:
        raise SchemaError('manifest carries no schema_version')
    if version not in SUPPORTED_VERSIONS:
        raise SchemaError(
            f'probe set schema version {version} cannot be read by this release, '
            f'which supports {list(SUPPORTED_VERSIONS)}. Upgrade oligominer.')

    for field in ('created', 'n_probes', 'columns', 'stages'):
        if field not in manifest:
            raise SchemaError(f'manifest is missing required field {field!r}')

    # success
    return manifest


# upgrades keyed by the version they read, each returning the next version up
MIGRATIONS = {}


def upgrade_manifest(manifest):
    """
    Upgrade a manifest to the current schema version.

    Args:
        manifest (dict): a manifest at any supported version.

    Returns:
        manifest (dict): the manifest at SCHEMA_VERSION.

    Raises:
        SchemaError: if no upgrade path reaches the current version.
    """
    version = manifest.get('schema_version')
    while version != SCHEMA_VERSION:
        migrate = MIGRATIONS.get(version)
        if migrate is None:
            raise SchemaError(
                f'no migration from probe set schema version {version} '
                f'to {SCHEMA_VERSION}')
        manifest = migrate(manifest)
        version = manifest['schema_version']

    # success
    return manifest


def dumps(manifest):
    """
    Serialize a manifest to JSON.

    Args:
        manifest (dict): the manifest.

    Returns:
        text (str): the JSON document.
    """
    # success
    return json.dumps(validate_manifest(manifest), indent=1, default=str)


def loads(text):
    """
    Read a manifest from JSON, upgrading it if it is older than this release.

    Args:
        text (str): the JSON document.

    Returns:
        manifest (dict): the manifest at the current schema version.

    Raises:
        SchemaError: if the document is unreadable at this version.
    """
    manifest = json.loads(text)
    validate_manifest(manifest)

    # success
    return upgrade_manifest(manifest)
