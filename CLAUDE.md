# CLAUDE.md — mento

Guidance for Claude Code (claude.ai/code) when working in this repository.

Reinforced concrete design Python package. Covers beams, slabs, sections, materials, rebar, and design code implementations (ACI 318-19, EN 1992-2004, CIRSOC 201-25). Uses strict mypy typing, ruff formatting, and pytest with coverage.

---

## Python environment

mento requires Python >= 3.12. Which interpreter to use depends on where the session runs.

**Windows (local development):** use the `rame-env` conda environment. It is the only one
with mento's dependencies installed, so the base `anaconda3` interpreter either fails on
import or runs against stale packages.

```
C:\Users\mihdi\anaconda3\envs\rame-env\python.exe
```

`.vscode/settings.json` registers the project with the conda environment manager
(`python-envs.pythonProjects`), so the VS Code test runner and IntelliSense resolve the same
interpreter.

**Linux / Claude Code on the web:** there is no conda environment. The SessionStart hook
`.claude/hooks/session-start.sh` builds `.venv` on the first Python >= 3.12 it finds and
installs mento with `-e ".[dev]"`. Use that interpreter:

```bash
.venv/bin/python
```

If `import mento` fails, run the hook by hand: `bash .claude/hooks/session-start.sh`.

Every command below writes the interpreter as `$PY`. Substitute the one for your platform:
the `rame-env` path on Windows, `.venv/bin/python` on Linux.

---

## Running tests

```bash
# Full suite (uses pyproject.toml addopts: --cov, --cov-report=html, --cov-report=term-missing)
$PY -m pytest tests/

# Single file, fast iteration (strip addopts to avoid --cov conflicts)
$PY -m pytest tests/test_beam.py --override-ini="addopts=" -v

# Single file with coverage
$PY -m pytest tests/test_beam.py --override-ini="addopts=" --cov=mento --cov-report=term-missing -q

# Single test by name
$PY -m pytest tests/test_beam.py::test_my_function --override-ini="addopts=" -v

# Coverage for one module only (Linux)
$PY -m pytest tests/ --override-ini="addopts=" --cov=mento --cov-report=term-missing -q | grep -E "beam|slab|rebar"
```

On Windows, pipe to `Select-String` instead of `grep`, and prefix the interpreter path with
the call operator: `& "C:\...\python.exe" -m pytest ...`.

> `--override-ini="addopts="` strips the default `--cov` flags from pyproject.toml. Required when running single files or adding custom `--cov` arguments to avoid argument conflicts.

---

## Linting and type checking

```bash
$PY -m ruff check . --fix     # lint, auto-fix
$PY -m ruff format .          # format
$PY -m mypy mento/            # strict type check
```

`make all` runs mypy, ruff and pytest in one go when the tools are on `PATH`.

Ruff config: 120-char line limit. MyPy config: `strict = true`, `allow_any_generics = true`. Package checked: `mento/`.

---

## Project structure

```
mento/
├── __init__.py             Lazy-import public API via __getattr__
├── _version.py             Package version
├── units.py                Pint unit registry (m, cm, mm, kN, kNm, MPa, kip, psi, etc.)
├── material.py             Concrete_ACI_318_19, Concrete_EN_1992_2004, Concrete_CIRSOC_201_25, SteelBar, SteelStrand
├── rebar.py                Rebar — bar database and selection logic
├── section.py              Section base class
├── rectangular.py          RectangularSection — geometry and cover calculations
├── beam.py                 RectangularBeam — design, check, and result accessors
├── slab.py                 OneWaySlab and Footing — one-way slab design
├── shear_wall.py           ShearWall — in-plane shear check/design
├── column.py               Column — geometry (shape, position, edge distances) for punching shear
├── punching.py             PunchingSlab, PunchingNode, Opening, Capital — two-way punching shear
├── punching_results.py     PunchingCheck — frozen result of a punching check
├── design_results.py       Public read-only result dataclasses (FlexureDesign, ShearDesign, SectionReinforcement, ...) — ADR-0001
├── precompute.py           SectionFloats — section geometry and materials as plain floats, converted once (ADR-0005)
├── forces.py               Forces(label, N_x, V_z, M_y, M_x) with pint units
├── node.py                 Node(section, forces) — drives check/design
├── settings.py             BeamSettings — metric/imperial defaults for design rules
├── results.py              Formatter, TablePrinter, DocumentBuilder — display helpers
├── beam_summary.py         BeamSummary — aggregate results for multiple beams
├── shear_wall_summary.py   ShearWallSummary — the same for walls
├── summary.py              Deprecated shim re-exporting BeamSummary (emits DeprecationWarning)
├── i18n.py                 set_language() — language of the detailed report output
├── plots/                  Matplotlib drawings: sections.py, walls.py, punching.py
├── reports/                Presentation layer: tables, views, documents, summaries, table_style, headings, punching
└── codes/
    ├── registry.py               DesignCode dataclass and the registered codes — every element dispatches through it
    ├── check_state.py            Check states held off the element; the report path copies them back
    ├── flexure_design.py         Flexure design engine shared by the codes (private)
    ├── aci_318_19/               code.py (registry entry) + equations/{flexure,shear,wall,punching}.py — floats only
    ├── en_1992_2004/             code.py + equations/
    ├── ACI_318_19_beam.py        Beam shear/flexure checks and design (functions typed `self: RectangularBeam`)
    ├── ACI_318_19_wall.py        Wall shear per ACI 318-19
    ├── ACI_318_19_punching.py    Punching checker — preconditions enforced, equations pending (Phase 2)
    ├── EN_1992_2004_beam.py      Beam shear/flexure per EN 1992-2004
    └── EN_1992_2004_punching.py  Punching checker — equations pending (Phase 5)
```

---

## Test structure

```
tests/
├── conftest.py           Shared fixtures + Agg matplotlib backend for the whole suite
├── test_architecture_boundaries.py   Enforces the layer rules (no units in equations, no report tables in codes/, ...)
├── test_aci_318_19_flexure_equations.py
├── test_aci_318_19_shear_equations.py
├── test_aci_318_19_wall_equations.py
├── test_en_1992_2004_equations.py
├── test_beam.py
├── test_beam_summary.py
├── test_design_results.py
├── test_slab.py
├── test_footing.py
├── test_material.py
├── test_rebar.py
├── test_section.py
├── test_rectangular.py
├── test_punching.py
├── test_shear_wall.py
├── test_shear_wall_summary.py
├── test_forces.py
├── test_node.py
├── test_units.py
├── test_settings.py
├── test_results.py
├── test_headings.py
├── test_table_style.py
├── test_i18n.py
├── test_documents.py               Single-document report, report targets (path / buffer), safe file names
├── test_lazy_imports.py            Subprocess: import + design must not import matplotlib/seaborn/IPython/docx
├── test_optional_presentation.py   Fallbacks when a presentation library is missing
└── test_init.py

scripts/
└── modules_testing.py    Manual/exploratory script; not part of the test suite
```

`pyproject.toml` is the only pytest config. A `tests/pytest.ini` would change the
rootdir and silently disable the coverage `addopts`, so there is none.

**Fixtures:** `conftest.py` holds the fixtures that several modules define
identically (`concrete_c25`, `steel_b500s`, `steel`) and sets the `Agg` backend
once, so no test module needs the `matplotlib.use("Agg")` incantation. Fixtures
tied to a specific validated example stay in the module that asserts against
them. Note that `beam_example_imperial` exists in both `test_beam.py` and
`test_rebar.py` with **different** geometry and settings — they are per-module on
purpose; do not hoist them.

---

## Architecture & key patterns

**Class hierarchy:**
```
Section → RectangularSection → RectangularBeam → OneWaySlab → Footing
                                               → ShearWall
PunchingSlab (standalone dataclass); PunchingNode(slab, column, forces) pairs it with a Column
```

**Unit-system detection:** `Concrete` auto-detects metric vs. imperial from `f_c` units (MPa → metric, psi → imperial). This propagates through `BeamSettings` and all `Forces` objects — never hard-code unit assumptions.

**Design code delegation:** elements never compare `concrete.design_code` against a string; they look the code up in `codes/registry.py` (`DesignCode`) and call its hooks (`check_shear`, `design_flexure`, `check_punching`, ...). Each code's entry lives in `codes/aci_318_19/code.py` / `codes/en_1992_2004/code.py`; the hooks are module-level functions typed as `self: RectangularBeam` in `codes/ACI_318_19_beam.py` and friends, which convert on entry and call the float-only clause functions in `equations/` (ADR-0002, ADR-0005). `tests/test_architecture_boundaries.py` fails the build if these rules are broken. `Concrete_CIRSOC_201_25` subclasses `Concrete_ACI_318_19` (same formulas, metric only, `design_code = "CIRSOC 201-25"`).

**`BeamSettings` sentinel pattern:** Unset fields use `_NOT_SET` so `__post_init__` can apply metric or imperial defaults conditionally based on the detected unit system.

**Lazy presentation imports:** matplotlib, seaborn, IPython and python-docx are imported *inside* the functions that use them (`plot()`, `reports/_notebook.py`, `DocumentBuilder`, the `*_doc` functions), never at module top level outside `mento/plots/`. `tests/test_lazy_imports.py` fails if `import mento` + a design pulls any of them in — mento-web (Pyodide) depends on it.

**`__init__.py` lazy loading:** Public API uses `__getattr__` so submodules are only imported on first attribute access. `TYPE_CHECKING` guards prevent circular imports.

---

## Key module notes

### beam

- `RectangularBeam` is a `@dataclass` that extends `RectangularSection`.
- Design code logic is factored into `codes/ACI_318_19_beam.py` and `codes/EN_1992_2004_beam.py`; beam delegates to these via direct function calls.
- `BeamSettings` controls spacing, bar diameter limits, and unit system — constructed first and passed to the beam.

### material

- Three concrete classes for different codes: `Concrete_ACI_318_19`, `Concrete_EN_1992_2004`, `Concrete_CIRSOC_201_25`.
- `SteelBar` defines yield strength and modulus.
- Unit system (metric vs. imperial) is determined by the concrete instance passed to `BeamSettings`.

### slab

- `OneWaySlab` covers one-way slab flexure and shear design.
- Shares the same material and unit infrastructure as beam.

### punching

- `PunchingSlab` takes `(concrete, steel_bar, h, c_c, outer_direction="x")`. Plan and status: [docs/architecture/punching-roadmap.md](docs/architecture/punching-roadmap.md).
- The top reinforcement is declared as bars, and `d` and ρ are **derived** from it — position 1 is the base mat, position 3 the extra bars over the column:
  ```python
  slab.set_rebar_x(d_b1=12*mm, s_b1=15*cm, d_b3=16*mm, s_b3=15*cm)
  slab.set_rebar_y(d_b1=12*mm, s_b1=15*cm)
  slab.d_x, slab.d_y, slab.d_avg, slab.rho_x, slab.rho_y, slab.rho_l, slab.A_s_x
  ```
- `rho_x` / `rho_y` / `d_avg` are read-only. `slab.d_avg = ...` raises; use `slab.set_effective_depth(d_x=..., d_y=...)`, which re-derives ρ against the depths it sets.
- With no bars declared, the depths fall back to a two-mat Ø16 / #5 layout (so `d_avg = h - c_c - 16 mm`) and `rho_x`/`rho_y` are `None`, not zero.
- `slab.data` / `node.data` render the inputs as Markdown (delegate to `mento/reports/punching.py`), like `beam.data`.
- The slab is checked through `PunchingNode(slab, column, forces, openings=None, capital=None)`. `Column` describes `shape` (`"rectangular"` / `"circular"`), `position` (`"interior"` / `"edge"` / `"corner"`), and edge distances when applicable.
- Forces at a node: **`V_z`** is the punching load (`Vu` / `VEd`), plus `M_x` and `M_y`. `N_x` is not used — see `docs/source/user_guide/local_axes.rst`.
- `node.check()` is wired: it dispatches through the code registry and returns a frozen `PunchingCheck` (`mento/punching_results.py`), not a DataFrame. The **equations are stubs that raise** — they are being recreated from a validated Calcpad sheet (Phase 2 ACI, Phase 5 EN). The checkers already refuse a force with no `V_z`, a capital or opening (Phase 3), and — for EN only — a slab with no declared ρ.
- `node.design()` raises via the registry: `design_punching` is Phase 4.

### units

- Import units directly: `from mento import m, cm, mm, kN, MPa` etc.
- `ureg` is the shared `UnitRegistry` — do not create additional registries.

### results / beam_summary

- `Formatter`: formats pint quantities for display.
- `TablePrinter`: renders pandas DataFrames as styled tables (Markdown/IPython).
- `DocumentBuilder`: builds Word (python-docx) report documents.
- `BeamSummary` (in `mento/beam_summary.py`; `mento/summary.py` is a deprecated shim): aggregates design results across multiple `RectangularBeam` instances.
  - `.check(capacity_check=False)` — DCR summary table for all beams; set `capacity_check=True` to zero forces and report capacities (MRd,top/bot or ØMn,top/bot) instead.
  - `.design()` — runs flexure + shear design for every beam and fills rebar columns.
  - `.flexure_results(capacity_check=False)` / `.shear_results(capacity_check=False)` — per-beam detailed check tables; `capacity_check=True` adds code-specific capacity columns.
  - `.results_detailed_doc(index=1)` — exports a Word document (`Beam_Summary_{design_code}.docx`) with full flexure/shear detail for the selected beam (1-based index) followed by summary tables for all beams. Saves to the current working directory.
  - `.export_design(path)` / `.import_design(path)` — round-trip the designed rebar to/from Excel.

---

## Running mento interactively (design / check from the CLI)

The full API walkthrough — materials, section, forces, node, reading the designed rebar,
printing results outside Jupyter, and the gotchas — lives in a rule file that Claude Code
imports:

@.claude/rules/interactive-usage.md

---

## Claude Code setup in this repo

```
.claude/
├── settings.json                  Shared permissions and the SessionStart hook (committed)
├── settings.local.json            Personal overrides, machine-specific paths (gitignored)
├── hooks/session-start.sh         Builds .venv and installs mento on Linux/web sessions
├── rules/interactive-usage.md     Imported by this file via @ — how to drive mento from a script
└── skills/
    ├── design-section/            Design or check a real section with mento
    └── prepare-release/           Bump version + CHANGELOG for a release branch
```

Keep this file under ~2000 words. Anything longer, or only needed for one kind of task,
belongs in `.claude/rules/` (imported with `@`) or in a skill.

---

## CI / GitHub Actions

Workflow: `.github/workflows/tests.yml`
- Runs on push/PR to `main`
- `tests` job: Python 3.12 and 3.13 on windows-latest, installing `-e ".[test]"`.
  Runs `pytest --cov=mento --cov-config=.coveragerc --cov-report=xml`; the 3.12 job
  uploads coverage to Codecov.
- `lint` job: ubuntu-latest, `ruff check .`, `ruff format --check .` and `mypy mento/`
  (no auto-fix in CI).
- `docs` job: ubuntu-latest, installs pandoc via apt, then
  `sphinx-build -b html -W --keep-going` so warnings fail the build.

The lint and docs jobs stay on ubuntu-latest because they are build infrastructure, not
platform coverage — and the docs job installs pandoc with apt. Only the test matrix
describes the platforms mento is verified on.

---

## Common patterns

**Bare object construction (unit-testing methods without __init__):**
```python
obj = object.__new__(MyClass)
obj.some_attr = value
obj.some_method()
```

**Checking coverage for one file:** see *Running tests* above — filter the
`--cov-report=term-missing` output with `grep` (Linux) or `Select-String` (Windows).
