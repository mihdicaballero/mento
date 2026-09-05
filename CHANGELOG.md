# Changelog

All notable changes to mento are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). From 1.0.0
the public API is stable: breaking changes need a major release. Before 1.0.0 it could
change between minor releases, which is what ADR-0003 documented at the time.

Releases before 0.5.0 were not tracked in a changelog; the entries below were reconstructed
from the release history and are summaries rather than complete lists.

## [Unreleased]

## [1.1.0] - 2026-09-02

A footing is a section mento now designs as it is built, and the results a check returns
carry the resistance their ratio was formed from. Both came from mako, which had to
correct the engine's output for the first and divide it back out for the second.

### Added

- **`Footing`: a one-way slab bearing on the ground.** Everything about it is a
  `OneWaySlab`; what changes is the minimum longitudinal reinforcement, and that change is
  the design codes' rather than the class's. `Section.support` (`"free"` / `"soil"`) is
  a ClassVar the codes read, so the clause lives in `codes/`: ACI 318-19 §9.6.1.1(b)
  grants the exemption and §13.3.1.2 substitutes the shrinkage and temperature steel of
  §24.4.3.2 on the gross section (CIRSOC 201-25 shares the clause); EN 1992-1-1 takes the
  halved geometric minimum of a foundation on every face and, on a face that is bending,
  the larger of that and the crack-control minimum of §7.3.2(2). A footing is detailed
  between 100 and 300 mm, both faces set out at one spacing (or the top at twice the
  bottom), and a thin one warns. mento does no geotechnical calculation: bearing,
  settlement, sliding, overturning, plan size and thickness stay with the caller, and
  a face the moment does not put in tension is left unreinforced, because whether a
  footing carries steel there is the consumer's call. A theory page and a user guide
  page open on exactly that. (#150, #151)

- **A designed slab reads as a spacing, and its bars may only sit so far apart.** The
  rebar search answers in groups of bars, and `OneWaySlab` applied that answer the beam's
  way, so a designed strip came out as `2Ø12 + 4Ø12` with the spacing that actually
  details a slab left at zero. Each layer is now spread over the strip when the design is
  applied and the spacing is what the slab stores — `str(slab.reinforcement.bottom)` is
  `Ø12 mm/17 cm` — with `RebarLayer.s` carrying it and reading `None` on a beam. Both
  codes also cap the centre-to-centre spacing (ACI 318-19 §7.7.2.3, EN 1992-1-1 §9.3.1.1)
  through a registry hook, so a lightly loaded strip no longer covers its area with a bar
  every half metre. Areas and DCRs are unchanged. (#149)

- **The frozen results carry the capacity the DCR was formed from.** `FlexureFaceCheck`
  and `FlexureFaceDesign` gain `M_capacity`; `ShearCheck` and `ShearDesign` gain
  `V_capacity`. One neutral name for every code: `ØMn` and `ØVn` (capped by `ØVmax`)
  under ACI 318-19 and CIRSOC 201-25, `MRd` and `VRd` under EN 1992-1-1. A caller
  printing the resistance next to the demand no longer divides the demand by a rounded
  ratio, gets a real number for a face reinforced against no demand, and needs none of
  the `_phi_M_n_*`, `_M_Rd_*`, `_phi_V_n` or `_V_Rd` attributes of the compatibility
  layer. On a design the capacity is the governing combination's, so `demand / DCR`
  gives it back there too. Brief: `docs/architecture/capacity-in-design-results.md`.
  (#152)

- **`A_s_calc`: the steel the moment alone asks for, next to `A_s_req`.** `A_s_req` on
  `FlexureFaceCheck` and `FlexureFaceDesign` has the minimum folded in — `max(mechanical,
  minimum)`, or the 4/3 rule under ACI — which is right for choosing bars and wrong for
  anchoring them: a development length scaled by `A_s,nec / A_s,prov` needs the
  mechanical area, and a footing mat governed by its minimum carries little of the
  stress the minimum is sized for. Both design codes now return that area before the
  `max` — `reinforcement_for_moment` or `A_s1_lim + A_s2` under EN, `A_s_calc` under
  ACI — and it reaches the results as `A_s_calc`, enveloped over the combinations like
  `A_s_req`. Zero with no moment; equal to `A_s_req` once the moment governs, compression
  steel included. `A_s_req` does not change meaning.

### Fixed

- **`flexure_check_results()` returned `A_s_req`, `A_s_min` and `A_s_max` as bare
  floats** in the design code's canonical units (mm² or in²), where `check_flexure()`
  returned them as quantities in cm² or in². Both entry points now read the check state
  through the same path and return quantities. (#152)

- **The EN crack-control minimum governs the thin footings, not the thick ones.** Two
  docstrings had it inverted; the calculation and the theory page were right all along.
  Written on `A_ct = b·h/2`, its ratio goes with `k/2`, and `k` decays from 1.00 to 0.65
  between 300 and 800 mm, so the geometric minimum takes over around h = 640 mm for C25
  with B500S. The trend is now pinned by tests. (#151)

## [1.0.1] - 2026-08-31

### Changed

- **`BeamSummary.check()` returns the columns the report prints, and only those.** It
  carried the section, the reinforcement, the required areas, the demands, the capacities
  and the three DCRs, and the Word report then selected a shorter set on its way to the
  page — so a notebook and the document disagreed about what the summary is. Both now
  show `Beam`, `b`, `h`, `As,top`, `As,bot`, `Av`, the governing demands, the three DCRs
  and `Ok?`: the shape the shear-wall summary already had, ending on the DCRs and the
  verdict.

  Gone from `check(capacity_check=False)`: `cc`, `As,req,top`, `As,req,bot`, `Av,req`,
  `Av,real`, and the capacity columns (`ØMn,top`/`MRd,top`, `ØMn,bot`/`MRd,bot`,
  `ØVn`/`VRd`). None of it is lost — `flexure_results()` and `shear_results()` report the
  required areas and the capacities per combination, which is where they mean something,
  and `check(capacity_check=True)` still returns the capacities.

  This removes columns from a public DataFrame, which under the policy above would call
  for a major release. It ships as a patch deliberately: 1.0.0 was a day old and the
  change is narrow. Recorded here rather than passed over quietly.

## [1.0.0] - 2026-08-31

The API is stable from this release on (ADR-0003): breaking changes need a major version
from here. What made 1.0 worth declaring is the architecture work below — checks that
return values instead of writing them to the section, so a caller can run one over many
sections and trust the answer.

### Added

- **`section.reinforcement`** — what a section carries right now, as a frozen
  `SectionReinforcement(bottom, top, transverse)`. It never raises, so a layout that was
  just set with `set_longitudinal_rebar_bot` can be read back without a check having run.
  The design-result objects were previously the only way to ask, and they are gated behind
  `DesignNotRunError`.
- **`beam.shear_check_results(forces)` and `beam.flexure_check_results(forces)`** — one
  frozen result per load combination, and no report built at all. Three to five times
  faster than the reporting path, and the numbers are identical; a test asserts that for
  both design codes.
- **`beam.shear_checks` and `beam.flexure_checks`** — every combination checked so far, so
  a caller that wants to envelope them differently no longer has to reach into the element.

### Changed

- **A check no longer writes to the section.** `shear_check_results` and
  `flexure_check_results` leave the element untouched — measured across both design codes,
  both checks and with and without stirrups: zero attributes changed, zero created. The
  reporting entry points still copy their results back, because the report tables read
  them off the element; that compatibility layer is deprecated and goes in a later release.
- **The shear check keeps the assumed stirrup diameter, as the flexure check already did.**
  On a section with no stirrups configured, `check_shear` used to drop the diameter the
  settings assume and recompute the effective depths on the section, so a `check_flexure`
  run afterwards reported a different DCR purely because of call order — 0.593 against
  0.582 on the same beam. Both checks now read the same effective depth. If you check a
  section that has no transverse reinforcement, say so with
  `set_transverse_rebar(n_stirrups=0, d_b=0, s_l=0)` and the depths follow.
- **`BeamSummary.check()` names its verdict column `Ok?`**, not `Status`.
- Adding a design code no longer means editing an element. A code declares itself in
  `mento/codes/<code>/code.py` and is found by walking that directory; nothing under
  `beam.py`, `shear_wall.py`, `rebar.py` or the summaries names a code any more.
- The Word reports read better on the page: tables are sized to their content rather than
  stretched across the line, forces show one decimal and DCRs two, the pass/fail cells are
  shaded green or red, and the summaries drop the capacity ticks that repeat what the DCR
  column beside them already says.

### Performance

- **A check is roughly ten times faster.** Shear and flexure on one section went from
  2.35 ms to 0.23 ms, so 20,000 sections take 4.6 s instead of 47 s. The equations run on
  plain floats and a section publishes its geometry and materials converted once
  (ADR-0005); pint stays at the boundary, where the inputs and the results are.
- Flexural design is four times faster: 306 ms to 72 ms, from the reinforcement search no
  longer building a `Quantity` per candidate.

### Fixed

- The report tables no longer write to the section they describe.


### Removed

- **Support for Python 3.10 and 3.11.** `requires-python` is now `>=3.12`, so pip refuses
  to install mento on the older interpreters instead of installing it and failing later.
  The classifiers and the conda recipe's `python_min` follow. Stay on 0.5.2 if you are
  pinned to 3.10 or 3.11.
- Ubuntu from the test matrix. Tests run on windows-latest against Python 3.12 and 3.13,
  down from eight jobs to two. mento is pure Python and nothing in it is
  platform-specific, but Linux is no longer verified on every pull request. The lint and
  docs jobs still run on ubuntu-latest, and the PyPI distributions are still built there.

### Fixed

- **Shear design ignored the spacing limit across the width of the section.** The stirrups
  were sized from the required area alone, so a wide beam came back with a single
  two-legged stirrup whose legs sat far further apart than ACI 318-19 Table 9.7.6.2.2 or
  EN 1992-1-1 9.2.2(8) allow — 44 cm against a 28 cm limit on a 50 cm section. The check
  reported the violation, but the design would not avoid it, so `design_shear` handed back
  a layout its own detailed report then marked as not compliant. The design now starts from
  the fewest legs the width admits and adds stirrups rather than only tightening the
  longitudinal spacing. Over a sweep of 168 width and demand combinations across the three
  codes, 95 designs were in violation and none are now. Closes
  [#94](https://github.com/mihdicaballero/mento/issues/94).
- The spacing across the width was computed with whatever stirrup diameter the previous
  pass had left on the beam instead of the one being tried, so the value stored for each
  candidate was off by the difference between the two diameters.

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

[Unreleased]: https://github.com/mihdicaballero/mento/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/mihdicaballero/mento/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/mihdicaballero/mento/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/mihdicaballero/mento/compare/v0.5.2...v1.0.0
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
