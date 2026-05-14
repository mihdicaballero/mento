# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# CLAUDE.md — mento

Reinforced concrete design Python package. Covers beams, slabs, sections, materials, rebar, and design code implementations (ACI 318-19, EN 1992-2004, CIRSOC 201-25). Uses strict mypy typing, ruff formatting, and pytest with coverage.

---

## Python environment

**Always use the `py312` conda environment.** The base `anaconda3` env should never be used for mento.

```
C:\Users\mihdi\anaconda3\envs\py312\python.exe
```

This is also set in `.vscode/settings.json` as `python.defaultInterpreterPath` so the VS Code test runner and IntelliSense use it automatically.

---

## Running tests

```powershell
# Full suite (uses pyproject.toml addopts: --cov, --cov-report=html, --cov-report=term-missing)
& "C:\Users\mihdi\anaconda3\envs\py312\python.exe" -m pytest tests/

# Single file, fast iteration (strip addopts to avoid --cov conflicts)
& "C:\Users\mihdi\anaconda3\envs\py312\python.exe" -m pytest tests/test_beam.py --override-ini="addopts=" -v

# Single file with coverage
& "C:\Users\mihdi\anaconda3\envs\py312\python.exe" -m pytest tests/test_beam.py --override-ini="addopts=" --cov=mento --cov-report=term-missing -q

# Check coverage for a specific module only
& "C:\Users\mihdi\anaconda3\envs\py312\python.exe" -m pytest tests/ --override-ini="addopts=" --cov=mento --cov-report=term-missing -q 2>&1 | Select-String "beam|slab|rebar"
```

> `--override-ini="addopts="` strips the default `--cov` flags from pyproject.toml. Required when running single files or adding custom `--cov` arguments to avoid argument conflicts.

---

## Linting and type checking

```powershell
# Ruff lint (auto-fix)
& "C:\Users\mihdi\anaconda3\envs\py312\python.exe" -m ruff check . --fix

# Ruff format
& "C:\Users\mihdi\anaconda3\envs\py312\python.exe" -m ruff format .

# MyPy (strict)
& "C:\Users\mihdi\anaconda3\envs\py312\python.exe" -m mypy mento/
```

Ruff config: 120-char line limit. MyPy config: `strict = true`, `allow_any_generics = true`. Package checked: `mento/`.

---

## Project structure

```
mento/
├── __init__.py           Lazy-import public API via __getattr__
├── _version.py           Package version
├── units.py              Pint unit registry (m, cm, mm, kN, kNm, MPa, kip, psi, etc.)
├── material.py           Concrete_ACI_318_19, Concrete_EN_1992_2004, Concrete_CIRSOC_201_25, SteelBar, SteelStrand
├── rebar.py              Rebar — bar database and selection logic
├── section.py            Section base class
├── rectangular.py        RectangularSection — geometry and cover calculations
├── beam.py               RectangularBeam — design, check, and visualization
├── slab.py               OneWaySlab — one-way slab design
├── forces.py             Forces(Vu, Mu, Nu, ...) with pint units
├── node.py               Node(x, y, z) with associated forces
├── settings.py           BeamSettings — metric/imperial defaults for design rules
├── results.py            Formatter, TablePrinter, DocumentBuilder — output and plotting
├── summary.py            BeamSummary — aggregate results for multiple beams
├── column.py             Column — geometry (shape, position, edge distances) for punching shear
├── punching.py           PunchingSlab — two-way punching shear check per ACI/EN
└── codes/
    ├── ACI_318_19_beam.py    Shear and flexure checks/design per ACI 318-19
    └── EN_1992_2004_beam.py  Shear and flexure checks/design per EN 1992-2004
```

---

## Test structure

```
tests/
├── test_beam.py
├── test_slab.py
├── test_material.py
├── test_rebar.py
├── test_section.py
├── test_rectangular.py
├── test_forces.py
├── test_node.py
├── test_units.py
├── test_settings.py
├── test_results.py
├── test_summary.py
├── test_init.py
└── modules_testing.py    Manual/exploratory testing helpers (not collected by pytest)
```

---

## Architecture & key patterns

**Class hierarchy:**
```
Section → RectangularSection → RectangularBeam
                             → OneWaySlab
PunchingSlab (standalone dataclass, uses Column)
```

**Unit-system detection:** `Concrete` auto-detects metric vs. imperial from `f_c` units (MPa → metric, psi → imperial). This propagates through `BeamSettings` and all `Forces` objects — never hard-code unit assumptions.

**Design code delegation:** `codes/ACI_318_19_beam.py` and `codes/EN_1992_2004_beam.py` contain module-level functions typed as `self: RectangularBeam`; `RectangularBeam` imports and calls them directly. `Concrete_CIRSOC_201_25` subclasses `Concrete_ACI_318_19` (same formulas, metric only, `design_code = "CIRSOC 201-25"`).

**`BeamSettings` sentinel pattern:** Unset fields use `_NOT_SET` so `__post_init__` can apply metric or imperial defaults conditionally based on the detected unit system.

**`__init__.py` lazy loading:** Public API uses `__getattr__` so submodules are only imported on first attribute access. `TYPE_CHECKING` guards prevent circular imports.

---

## Key module notes

### beam

- `RectangularBeam` is a `@dataclass` that extends `RectangularSection`.
- Design code logic is factored into `codes/ACI_318_19_beam.py` and `codes/EN_1992_2004_beam.py`; beam delegates to these via direct function calls.
- `BeamSettings` controls spacing, bar diameter limits, and unit system — constructed first and passed to the beam.
- Tests use `MPLBACKEND=Agg` implicitly (matplotlib doesn't need a display). If plots fail in CI, set `MPLBACKEND=Agg`.

### material

- Three concrete classes for different codes: `Concrete_ACI_318_19`, `Concrete_EN_1992_2004`, `Concrete_CIRSOC_201_25`.
- `SteelBar` defines yield strength and modulus.
- Unit system (metric vs. imperial) is determined by the concrete instance passed to `BeamSettings`.

### slab

- `OneWaySlab` covers one-way slab flexure and shear design.
- Shares the same material and unit infrastructure as beam.

### punching

- `PunchingSlab` takes `(concrete, steel_bar, h, c_c, rho_x, rho_y)`.
- `d_avg` is computed as `h - c_c - 16 mm` (metric) or `h - c_c - 5/8 in` (imperial); override after construction if needed: `slab.d_avg = custom_value`.
- Requires a `Column` instance describing `shape` (`"rectangular"` / `"circular"`), `position` (`"interior"` / `"edge"` / `"corner"`), and edge distances when applicable.

### units

- Import units directly: `from mento import m, cm, mm, kN, MPa` etc.
- `ureg` is the shared `UnitRegistry` — do not create additional registries.

### results / summary

- `Formatter`: formats pint quantities for display.
- `TablePrinter`: renders pandas DataFrames as styled tables (Markdown/IPython).
- `DocumentBuilder`: builds Word (python-docx) report documents.
- `BeamSummary`: aggregates design results across multiple `RectangularBeam` instances.
  - `.check(capacity_check=False)` — DCR summary table for all beams; set `capacity_check=True` to zero forces and report capacities (MRd,top/bot or ØMn,top/bot) instead.
  - `.design()` — runs flexure + shear design for every beam and fills rebar columns.
  - `.flexure_results(capacity_check=False)` / `.shear_results(capacity_check=False)` — per-beam detailed check tables; `capacity_check=True` adds code-specific capacity columns.
  - `.results_detailed_doc(index=1)` — exports a Word document (`Beam_Summary_{design_code}.docx`) with full flexure/shear detail for the selected beam (1-based index) followed by summary tables for all beams. Saves to the current working directory.
  - `.export_design(path)` / `.import_design(path)` — round-trip the designed rebar to/from Excel.

---

## CI / GitHub Actions

Workflow: `.github/workflows/tests.yml`
- Runs on push/PR to `main`
- Matrix: Python 3.10, 3.11, 3.12 on ubuntu-latest and windows-latest
- Steps: `pip install pytest pytest-cov`, `pip install -r requirements_dev.txt`
- Lint: `ruff check .` (no auto-fix in CI)
- Tests: `pytest --cov=mento --cov-config=.coveragerc --cov-report=xml`
- Coverage uploaded to Codecov

---

## Common patterns

**Bare object construction (unit-testing methods without __init__):**
```python
obj = object.__new__(MyClass)
obj.some_attr = value
obj.some_method()
```

**Checking coverage for one file:**
```powershell
& "C:\Users\mihdi\anaconda3\envs\py312\python.exe" -m pytest tests/ --override-ini="addopts=" --cov=mento --cov-report=term-missing -q 2>&1 | Select-String "filename_or_module"
```

**Running a single test by name:**
```powershell
& "C:\Users\mihdi\anaconda3\envs\py312\python.exe" -m pytest tests/test_beam.py::test_my_function --override-ini="addopts=" -v
```
