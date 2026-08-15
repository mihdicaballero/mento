# Changelog

All notable changes to mento are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). While the
version number is below 1.0, the public API may still change between minor releases.

Releases before 0.5.0 were not tracked in a changelog; the entries below were reconstructed
from the release history and are summaries rather than complete lists.

## [Unreleased]

### Added

- Public design results API: `beam.flexure_design` and `beam.shear_design` return plain,
  frozen data objects with the reinforcement a check or design produced, so results no
  longer have to be read from private attributes such as `_A_s_bot` or `_stirrup_s_l`.
  Reading either before running a check raises `DesignNotRunError`. Required areas and
  DCRs are the envelope over every load combination checked, so they describe the
  combination that governs each face. See
  [Design results](https://mento-docs.readthedocs.io/en/latest/user_guide/design_results.html).
- A [citing guide](https://mento-docs.readthedocs.io/en/latest/getting_started/citing.html)
  in the documentation, covering which version to cite and a BibTeX entry.
- A conda recipe under `conda-recipe/`, kept in step with `pyproject.toml`, ready to be
  submitted to conda-forge. The submission and update procedure, and the one-time Zenodo
  setup that gives releases a DOI, are documented in CONTRIBUTING.
- Two example notebooks in Spanish, design and check of a rectangular beam under
  CIRSOC 201-2025.

### Changed

- `mento.summary` was renamed to `mento.beam_summary`, matching `mento.shear_wall_summary`.
  The old module still works and emits a `DeprecationWarning`; `from mento import
  BeamSummary` is unaffected.

### Fixed

- **`pip install mento` on Google Colab.** The IPython requirement was `>=8.0`, while Colab
  pins ipython to 7.34.0, so installing mento there could not resolve without upgrading
  IPython out from under the running session. The floor is now `>=7.34`; mento only uses
  `IPython.display.Markdown` and `display`, which both bounds cover.

## [0.5.0] - 2026-08-01

Infrastructure release. No changes to the design calculations, and no changes to the public
API beyond one addition.

### Added

- `mento.__version__` exposes the installed package version.
- Contributing guide, code of conduct, security policy and citation metadata at the
  repository root, plus issue and pull request templates.
- Python 3.13 is tested and supported.
- Optional dependency groups: `mento[test]`, `mento[dev]` and `mento[docs]`.
- Continuous integration now runs mypy, `ruff format --check` and a documentation build in
  addition to the test suite.
- Releases are published to PyPI automatically from a GitHub Release, using trusted
  publishing.
- A disclaimer of professional responsibility in the README and the documentation.

### Changed

- Dependencies are declared in `pyproject.toml` instead of being read from
  `requirements.txt`, with lower bounds rather than exact pins.
- **numpy is no longer pinned to 1.26.4.** The pin blocked numpy 2.x and caused install
  conflicts in environments with other scientific packages. The test suite passes on both
  numpy 1.26 and numpy 2.5.
- The ruff pre-commit hook was updated from v0.6.2 to v0.15.22, so the hook, CI and a local
  `ruff format` all agree. Eleven files that had drifted were reformatted.
- mypy configuration moved from `mypy.ini` into `pyproject.toml` as a single source. Modules
  that do not pass strict checking yet are listed explicitly, so the modules that are clean
  cannot regress.
- Read the Docs builds on Python 3.12 and installs the `docs` extra.

### Fixed

- `RectangularSection._ax` was annotated as `plt.Axes`, which is not a valid type; it is now
  `matplotlib.axes.Axes`.
- The `LICENSE` file still contained the unfilled `<year>` and `<author>` placeholders.
- The test workflow checked for `requirements.txt` but installed `requirements_dev.txt`.

### Removed

- Dead `[tool.black]` and `[tool.pylint]` configuration, an empty `.github/workflow` file,
  and unused documentation dependencies (dask, xarray, sparse, mip and others that mento
  never imported).

## [0.4.1] - 2026-05-17

### Added

- Shear wall summary, with results across multiple walls and a Word export.

## [0.4.0] - 2026-05-15

### Fixed

- Minimum reinforcement calculation for ACI 318-19.

## [0.3.6] - 2026-05-15

### Added

- Shear wall check and design for ACI 318-19 and CIRSOC 201-25.

## [0.3.5] - 2026-01-04

## [0.3.4] - 2025-12-20

### Added

- One way slab check and design.

## [0.3.0] - 2025-11-09

### Added

- Flexure check and design for EN 1992-2004.

## [0.2.8] - 2025-10-21

## [0.2.7] - 2025-08-17

## [0.2.6] - 2025-07-27

## [0.2.5] - 2025-03-23

First public release on PyPI: rectangular concrete beam check and design for flexure and
shear under ACI 318-19 and CIRSOC 201-25, unit aware calculations, results as pandas
DataFrames, and Word calculation reports.

[Unreleased]: https://github.com/mihdicaballero/mento/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/mihdicaballero/mento/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com/mihdicaballero/mento/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/mihdicaballero/mento/compare/v0.3.6...v0.4.0
[0.3.6]: https://github.com/mihdicaballero/mento/compare/v0.3.5...v0.3.6
[0.3.5]: https://github.com/mihdicaballero/mento/compare/v0.3.4...v0.3.5
[0.3.4]: https://github.com/mihdicaballero/mento/compare/v0.3.0...v0.3.4
[0.3.0]: https://github.com/mihdicaballero/mento/compare/v0.2.8...v0.3.0
[0.2.8]: https://github.com/mihdicaballero/mento/compare/v0.2.7...v0.2.8
[0.2.7]: https://github.com/mihdicaballero/mento/compare/v0.2.6...v0.2.7
[0.2.6]: https://github.com/mihdicaballero/mento/compare/v0.2.5...v0.2.6
[0.2.5]: https://github.com/mihdicaballero/mento/releases/tag/v0.2.5
