# mento architecture roadmap

Status: accepted (2026-08-26). Companion decision records live in [`adr/`](adr/).

This document describes where mento's internal architecture is today, what the two
closest open-source reference projects do, the target architecture we are converging
to, and the migration phases to get there without ever breaking the test suite.

It exists because mento needs to scale in two directions at once:

- **More design codes and more modules** (punching, columns, more of CIRSOC, future
  codes) contributed by structural engineers who should not need to understand the
  whole package to add one verification.
- **mako**: running mento as a computation engine over ~1000 stations of a continuous
  beam, per load combination, accumulating per-section results and then unifying them
  into a continuous rebar layout that is itself checked at the end. That loop needs
  stateless checks returning immutable results, with no plotting or Word machinery
  loaded along the way.

mento is deliberately production-oriented: pint units end to end, Word reports, and
design plots are differentiators for office use and they stay. What changes is *where*
they live, not *whether* they exist.

---

## 1. Where mento is today

The package works and is well validated — `tests/` is mostly a numeric regression
suite against spreadsheet/ETABS/code-worked examples (hundreds of `pytest.approx`
asserts with relative tolerances). That suite is the safety net that makes everything
below feasible. The structural problems are three, and they are specific:

### 1.1 `codes/` is not a layer — it is `RectangularBeam` methods exiled to other files

Every function in `codes/ACI_318_19_beam.py` and `codes/EN_1992_2004_beam.py` takes
`self: "RectangularBeam"` and writes private attributes back onto the beam
(81 `self._x = ...` assignments in the ACI module, 77 in the EN module). Across
`beam.py` and the two code modules, **~164 distinct private attribute names are
written on the same instance from three different files**. The import cycle is broken
with `TYPE_CHECKING`, but the coupling is total: the underscore does not mean
"private" here, it is the data bus of the system. Adding a new code today means
adding another file that mutates the same object.

### 1.2 Results are mutable state, not values

Each check overwrites the attributes left by the previous combination. Two artifacts
in the codebase are the internal diagnosis of this, written by ourselves:

- The *envelope* accumulators (`_update_flexure_envelope` / `_update_shear_envelope`
  in `beam.py`), whose docstring says it plainly: *"a reader of `_A_s_req_bot` after
  checking several combinations gets whichever one happened to run last"*.
- `design_results.py`, a clean frozen-dataclass facade — but implemented as a
  *read-only adapter* over the beam's private state, complete with a
  `DesignNotRunError` to detect "the beam has not been mutated yet".

Both are patches on top of the cause rather than its removal.

### 1.3 Presentation is smeared across every layer

- `beam.py` (2135 lines) is ~30 % matplotlib plotting and ~360 lines of Word
  document assembly.
- The `_initialize_dicts_*` functions inside the code modules build report tables
  (UI labels, rounding, ✅ emojis): ~29 % of the ACI module and ~38 % of the EN
  module is table formatting, not code equations. The EN/ACI modules even call a
  string formatter that lives on the beam.
- The one clean boundary today is `results.py` (`Formatter`, `TablePrinter`,
  `DocumentBuilder`), which does not know about beam.

### 1.4 Hygiene findings (cheap to fix, tracked in Phase 0)

- `build/lib/mento/` (a stale copy of the package) and `docs/build/` are committed.
- Two contradictory pytest configs: `pyproject.toml` (`testpaths = ["tests"]`) and
  `tests/pytest.ini` (`testpaths = test`, a nonexistent directory).
- No `conftest.py`; beam fixtures are re-declared per test file.
- `tests/modules_testing.py` (829 lines) contains no tests — it is a manual
  exploration script living inside `tests/`.
- 17 generated `.docx`/`.xlsx` artifacts are versioned under `docs/source/examples/`.

---

## 2. Reference architecture: what the two closest projects converge on

We studied [concrete-properties](https://github.com/robbievanleeuwen/concrete-properties)
(Robbie van Leeuwen) and [structuralcodes](https://github.com/fib-international/structuralcodes)
(fib International). They diverge in mechanisms but **agree on the structural
decisions** — the strongest available signal that these are the right ones:

| Decision | concrete-properties | structuralcodes | mento today |
| --- | --- | --- | --- |
| Code-agnostic engine | `ConcreteSection` knows no code; `DesignCode` wraps it | Section/fiber engine knows no code | `if design_code == ...` repeated in 4 beam methods |
| Code equations as pure functions | (equations live in per-code classes) | `codes/mc2010/`, `codes/ec2_2004/`: `float → float` functions, clause in docstring | functions mutate the beam |
| Results = immutable dataclasses, one module | `results.py`, 8 dataclasses (`UltimateBendingResults`, ...) | `core/_section_results.py` | `design_results.py` exists, but as read-only facade |
| Presentation separated | plots on result objects, centralized via `post.py` `plotting_context` | zero `plot_*` in core; matplotlib not even a dependency | smeared across beam.py and codes/ |
| Tests = code-table values, parametrized | validation + physical-invariance tests | `parametrize` with tabulated expected values, `rel_tol=0.001` | **already there** (biggest existing asset) |
| Code registry | explicit object per code | `_DESIGN_CODES` dict + module protocol `__title__`/`__year__`/`__materials__` | `elif` chains |

What mento takes from each:

- **From structuralcodes**: the layered architecture (code equations → materials →
  geometry/sections → calculators; each layer depends only on lower ones), pure
  per-clause equation functions, and the minimal module protocol + registry for
  adding codes.
- **From concrete-properties**: units as a *presentation* concern (`UnitDisplay`
  injected into result objects, engine works in consistent base units), the
  `plotting_context` pattern that centralizes matplotlib boilerplate, and the
  per-code wrapper object for safety factors and check orchestration.

What mento does **not** copy: structuralcodes' academic minimalism (no reports, no
plots, floats-only results). mento's users need assembled, report-ready results.

---

## 3. Target architecture

Five layers, dependency arrows only pointing downward:

```
presentation      plots/, reports/, results.py      — consumes result dataclasses only
elements          beam.py, slab.py, shear_wall.py,  — geometry + thin orchestration;
                  punching.py                          no code logic, no presentation
codes             codes/<code>/equations/*.py       — pure functions, clause in docstring
                  codes/<code>/checker.py           — per-code class: φ/γ factors,
                                                       orchestration, result assembly
results           design_results.py (extended)      — ALL results as frozen dataclasses,
                                                       RETURNED by checks
foundation        units.py, material.py, forces.py, — leaves; already clean today
                  settings.py, section.py, i18n.py
```

### 3.1 The code layer (hybrid pattern — see ADR-0002)

Per design code, a subpackage with two kinds of modules:

```python
# codes/aci_318_19/equations/shear.py — ONLY formulas; cite the clause; tested
# against the code's own tables.
def V_c(f_c, lambda_, rho_w, b_w, d, sigma_Nu) -> Quantity:
    """Concrete shear strength — ACI 318-19, Table 22.5.5.1(c)."""

def s_max(d: Quantity, V_s: Quantity, threshold: Quantity) -> Quantity:
    """Maximum stirrup spacing — ACI 318-19 §9.7.6.2.2."""
```

```python
# codes/aci_318_19/checker.py — ONLY factors + orchestration + result assembly.
class ACI318_19:
    phi_v = 0.75

    def check_shear(self, section, force, settings) -> ShearCheckResult:
        Vc = eq.V_c(f_c=section.concrete.f_c, b_w=section.width, ...)
        Vs = eq.V_s(...)
        return ShearCheckResult(V_c=Vc, phi_V_n=self.phi_v * (Vc + Vs), ...)
```

The discipline that keeps this honest: **the checker does not compute — it calls
equations and assembles the result.** If a checker method starts containing algebra,
that algebra belongs in `equations/`.

Why hybrid and not pure-functions-only (structuralcodes) or class-only
(concrete-properties): the equations get the testability and the low contribution
threshold of pure functions — an engineer who knows ACI but little Python can read,
verify, and contribute one — while the rich result assembly that Word reports and
plots need gets an explicit home instead of being spread across loose functions.

### 3.2 Results as values (see ADR-0001)

Checks and designs **return** frozen dataclasses (`ShearCheckResult`,
`FlexureCheckResult`, `FlexureDesign`, `ShearDesign`, ...). Nothing mutates the
element. Consequences:

- Envelopes across combinations become plain code over a list of results
  (`max(r.DCR for r in results)`) instead of stateful accumulators.
- `DesignNotRunError` disappears by construction: if you have a result object, the
  design ran.
- Result objects carry everything reports need (demands, capacities, ratios,
  governing equations, intermediate values), so presentation needs no private-state
  access.

### 3.3 The mako scenario as design validation

The target architecture must make this loop trivial, stateless, and
parallelizable — with no matplotlib/docx import anywhere on the path:

```python
code = ACI318_19()
results = [
    code.check_shear(section, force, settings)
    for section in stations          # ~1000 stations of a continuous beam
    for force in combos              # per load combination
]
# mako then: envelope per station → propose continuous rebar layout →
# re-check the unified layout with the same check functions.
```

If a change makes this loop need a mutable element, a display import, or per-call
state cleanup, the change is wrong.

### 3.4 Presentation as its own layer (see ADR-0004)

Word reports, plots, styled tables, and pint display formatting are kept as
first-class features — implemented in modules that consume result dataclasses:

- `_initialize_dicts_*` table builders move out of `codes/` into the report layer.
- Section/rebar plotting moves out of `beam.py` into a plotting module (adopting a
  `plotting_context`-style helper to centralize matplotlib boilerplate).
- Display units follow concrete-properties: calculations run in pint quantities;
  *rounding and unit choice for display* belong to the presentation layer.

---

## 4. Migration phases

Incremental; the package imports and the full test suite passes after every phase.
Each phase is independently mergeable.

### Phase 0 — Hygiene (no risk, do first)

- Remove `build/` and `docs/build/` from version control; extend `.gitignore`.
- Delete `tests/pytest.ini`; `pyproject.toml` is the single pytest config.
- Add `tests/conftest.py`; move the shared beam/material fixtures there.
- Move `tests/modules_testing.py` out of `tests/` (e.g. `scripts/` or delete).
- Stop versioning generated `.docx`/`.xlsx` under `docs/source/examples/`.

### Phase 1 — Extract pure equations

Inside `codes/`, split code formulas out of the mutating orchestrators into
`equations/` modules (pure in / pure out, clause in docstring). The existing
orchestrators keep mutating the beam **but now call the pure functions**. Add
parametrized tests against code tables for each extracted equation. No public
behavior changes; this phase only *creates the layer that does not exist yet*.

### Phase 2 — Checks return result dataclasses

Extend `design_results.py` with check-result types; make check/design functions
build and return them. Private beam attributes are kept synchronized as a
deprecated compatibility layer (`DeprecationWarning` on documented ones). The
envelope accumulators are deleted — enveloping becomes a pure operation over the
returned list. This is the phase that unblocks mako.

### Phase 3 — Extract presentation

Move `_initialize_dicts_*` out of the code modules into the report layer; move
matplotlib and docx code out of `beam.py`. After this phase, `codes/` contains only
equations + checkers, and elements contain only geometry + thin orchestration.

### Phase 4 — Code registry

Formalize the per-code subpackage protocol (structuralcodes-style
`__title__` / `__year__` / `__materials__` metadata + a registry dict) so CIRSOC
and future codes plug in as subpackages without touching elements. Only worth doing
once phases 1–3 exist; a registry over today's structure would formalize the
coupling instead of removing it.

---

## 5. Versioning policy (see ADR-0003)

Pre-1.0, following the policy structuralcodes documents publicly:

- **Patch** releases: bug fixes and backwards-compatible features.
- **Minor** releases: **may contain breaking changes**, always listed in the
  changelog with migration notes and, where cheap, a deprecation shim for one
  release cycle.
- **1.0.0** marks the stable API; SemVer from then on.

This replaces the previous informal rule of never breaking the public API. It is
what makes phases 2–3 affordable: the documented-private attributes
(`beam._A_s_bot` and friends) can be deprecated and later removed in a minor
release instead of being supported forever.

---

## 6. What we deliberately do NOT do

- No hexagonal architecture, no dependency injection framework, no repository
  pattern, no speculative ABCs. If a pattern forces three classes to implement one
  shear check, the pattern is wrong for this package.
- No registry machinery before the layers exist (Phase 4 is last on purpose).
- No rewrite. Every phase lands on green tests; the numeric validation suite is the
  contract.

**Success criterion** for the whole roadmap: a contributor adds a new verification
(say, torsion per ACI 318-19) by copying the pattern of an existing one — an
`equations/` module with clause-cited functions, a checker method returning a
result dataclass, parametrized tests against the code's tables — without asking the
maintainer how the rest of the system works.
