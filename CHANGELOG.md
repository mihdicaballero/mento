# Changelog

All notable changes to mento are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). While the
version number is below 1.0, the public API may still change between minor releases.

Releases before 0.5.0 were not tracked in a changelog; the entries below were reconstructed
from the release history and are summaries rather than complete lists.

## [Unreleased]

### Removed

- **Support for Python 3.10 and 3.11.** `requires-python` is now `>=3.12`, so pip refuses
  to install mento on the older interpreters instead of installing it and failing later.
  The classifiers and the conda recipe's `python_min` follow. Stay on 0.5.2 if you are
  pinned to 3.10 or 3.11.
- Ubuntu from the test matrix. Tests run on windows-latest against Python 3.12 and 3.13,
  down from eight jobs to two. mento is pure Python and nothing in it is
  platform-specific, but Linux is no longer verified on every pull request. The lint and
  docs jobs still run on ubuntu-latest, and the PyPI distributions are still built there.

## [0.5.2] - 2026-08-25

### Added

- Spanish detailed reports. `mento.set_language("es")` switches
  `flexure_results_detailed()`, `shear_results_detailed()` and their `_doc()` counterparts
  to Spanish, for beams, one-way slabs and shear walls, in the console and in the generated
  Word documents. English remains the default, and `mento.get_language()` and
  `mento.available_languages()` report the current and the available choices. Variable
  names, units, the design code designation and the generated file names are not
  translated. A label with no translation is written in English rather than raising. See
  [Report language](https://mento-docs.readthedocs.io/en/latest/user_guide/language.html).
  Closes [#79](https://github.com/mihdicaballero/mento/issues/79) and
  [#126](https://github.com/mihdicaballero/mento/issues/126).
- A DOI. Releases are archived on Zenodo, and
  [10.5281/zenodo.21956634](https://doi.org/10.5281/zenodo.21956634) always resolves to the
  latest one. It is in `CITATION.cff`, in the README badge and in the citing guide.
- A [Theory section](https://mento-docs.readthedocs.io/en/latest/theory/index.html) in the
  documentation, with one page per element and design code, so the tool can be audited
  against the codes it implements. Each page names what is out of scope, the places where
  mento takes a position the code leaves open, and a table mapping every check to the test
  that pins it and the external source it was verified against.

### Changed

- Flexure output is labelled with the symbols of the active design code. An EN 1992-2004
  result was reported with ACI symbols — `M_u` and `\phi M_n` in the markdown line, and
  `Mu,top` / `Mu,bot` heading the design forces of the detailed table — where the Eurocode
  writes `M_Ed` and `M_Rd`. Only the limiting-case branch was affected, which is the one
  taken when no explicit force is passed. ACI 318-19 and CIRSOC 201-25 are unchanged.
- The conda recipe is pinned to 0.5.1 and its checksum.

### Fixed

- **EN 1992-2004 flexural design produced layouts that its own check then rejected.** Over
  a sweep of 288 combinations, 29 designs came back with DCR > 1. Several independent
  defects: the lever arm applied `lambda` twice, under-sizing the tension steel by about
  2%; the ductility limits mixed the neutral-axis and block-depth conventions, leaving
  `M_lim` about 20% too high; the redistribution limit ignored `k_3`/`k_4` above C50/60;
  `M_Rd` counted compression steel the section did not need, which made capacity
  non-monotonic; and the top face was sized with the bottom face's effective depth.
  Cross-checked against the Concise Eurocode 2 closed form and plain equilibrium. Flexural
  design is now driven by one shared strategy for both codes.
- **The ACI size effect factor was not capped at 1.0.** ACI 318-19 Eq. 22.5.5.1.3 writes
  `lambda_s = sqrt(2/(1 + d/10in)) <= 1.0`. Without the cap, sections with `d` below about
  250 mm got a factor above 1, inflating `V_c` instead of reducing it — the opposite of
  what the provision is for, and on the unsafe side. It only applies in the
  `A_v < A_v_min` branch, so in practice this moves slabs; beams under test are unaffected.
- Checking a section with no reinforcement for ACI shear raised `ZeroDivisionError`
  instead of reporting an insufficient section. With `rho_w` at zero, Table 22.5.5.1 puts
  `phi*V_n` at exactly zero and the DCR division blew up. It now reports `DCR = inf`,
  matching the guard already in the wall module. EN is unaffected, since `V_Rd,c` carries
  the `v_min` floor of 6.2.2(2).
- A one-way slab is no longer reported as a beam. `OneWaySlab` inherited the titles of
  `RectangularBeam`, so its detailed reports were headed `BEAM FLEXURE DETAILED RESULTS`
  and its Word file named `Beam S1 flexure check ACI 318-19.docx`. Slabs now use their own
  wording and file name; beams are unchanged.
- The forces table of the detailed flexure report is headed `Design forces`, matching the
  shear report and the Word output. It was `Design_forces` in the console output only.
- The Word reports spell the `Limit checks` heading the same way everywhere. Flexure used
  `Limit Checks` and shear used `Limit checks`, in both the element and the summary
  documents.

## [0.5.1] - 2026-08-15

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

[Unreleased]: https://github.com/mihdicaballero/mento/compare/v0.5.2...HEAD
[0.5.2]: https://github.com/mihdicaballero/mento/compare/v0.5.1...v0.5.2
[0.5.1]: https://github.com/mihdicaballero/mento/compare/v0.5.0...v0.5.1
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
