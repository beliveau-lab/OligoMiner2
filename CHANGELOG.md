# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!-- docs-include-start -->

## [Unreleased]

## [0.0.4] - 2026-09-17

First release under the packaging standard. `0.0.1.dev1` was a development release and could only
be installed with `--pre`; this one can be pinned normally.

`0.0.2` and `0.0.3` were tagged but never published: each release run failed before the approval
gate, and neither version number can be reused, because the build is not byte-reproducible and a
rebuilt artifact no longer matches the hash already recorded in the index.

### Added

- `pyarrow` is now a declared dependency. `models/retrain.py` and the l4t seqprops loader both read
  and write parquet, so those call sites raised `ImportError` for anyone who installed the package.

### Changed

- `requires-python` is `>=3.10`, replacing `>=3.8,<3.12`. The upper bound excluded three
  maintained Python versions, and nothing tested 3.8 or 3.9.
- Package metadata: the licence and trove classifiers are declared rather than commented out, so
  PyPI shows them; `biopython` is no longer listed twice; project URLs carry a repository, issues
  and changelog link.
- Development tooling moved from a published extra to a dependency group, so installing the
  package no longer offers it to users.

### Fixed

- Probe mining raised `OverflowError` on numpy 2. The dinucleotide sentinel `-1` was written into
  an unsigned array; earlier numpy silently wrapped it to 255.

[Unreleased]: https://github.com/beliveau-lab/OligoMiner2/compare/v0.0.4...HEAD
[0.0.4]: https://github.com/beliveau-lab/OligoMiner2/releases/tag/v0.0.4
