# mento architecture roadmap

Status: accepted (2026-08-26); amended 2026-08-29 with the measured performance
baseline (§1.5), the equation-signature decision (ADR-0005), the Phase 2 split,
per-phase exit criteria, and the performance track. Companion decision records
live in [`adr/`](adr/).

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

### 1.4 Hygiene findings (resolved — Phase 0, 2026-08-29)

- ~~`build/lib/mento/` and `docs/build/` are committed.~~ Already untracked when
  Phase 0 ran; `.gitignore`'s `build/` pattern matches at any depth. The stale
  copies remain on disk as local build artifacts only.
- Two contradictory pytest configs. **Worse than diagnosed:** `tests/pytest.ini`
  won — running `pytest tests/` set rootdir to `tests/` and used it, so
  `pyproject.toml`'s coverage `addopts` were silently ignored. Only CI, which
  passes `--cov` explicitly, was measuring coverage. `pytest.ini` is deleted;
  `pyproject.toml` is now the configfile.
- No `conftest.py`; fixtures re-declared per test file. Three fixtures were
  declared identically in two modules each (`concrete_c25`, `steel_b500s`,
  `steel`) and now live in `tests/conftest.py`, which also sets the `Agg`
  backend that four modules were each setting by hand.
- `tests/modules_testing.py` (829 lines) contains no tests. Moved to
  `scripts/modules_testing.py`.
- ~~17 generated `.docx`/`.xlsx` artifacts are versioned.~~ Already untracked; the
  `*.docx`/`*.xlsx`/`*.pdf` ignore rules cover them.

One trap found and deliberately left alone: `beam_example_imperial` is defined in
both `test_beam.py` and `test_rebar.py` with **different** geometry and settings.
Same name, different fixture — hoisting either into `conftest.py` would silently
change the other module's expected values.

### 1.5 Performance baseline (measured 2026-08-29)

Profiled on rame-env / Python 3.12, commit 6299b78 (ACI 318-19, 30×50 section,
M = 120 kNm); consistent with an independent measurement of ~331 ms on another
machine:

| Call | Cost today |
| --- | --- |
| `design_flexure` | 375 ms |
| `check_shear` | 8.3 ms |
| `check_flexure` (designed beam) | 7.3 ms |

Two facts anchor the plan:

- **pint is the multiplier.** ~65 % of profiler self-time is pint-internal
  (~967k pint calls vs ~2.8k mento calls per design). A representative flexure
  chain runs 1071 µs/call with pint + string units, 313 µs with prebuilt unit
  objects, **0.9 µs with plain floats** (~1200×). This is what decides the
  equation signature — see ADR-0005.
- **The design cost is rebar selection, not equations.** 97 % of
  `design_flexure` is the combinatorial search in `rebar.py` (hundreds of
  spacing/penalty evaluations, each doing pint arithmetic). The flexure
  equations cost ~6 ms. So the mako loop — which *checks* per station and
  designs rebar once — costs ~8 ms/station today, not 375: ~2.7 min serial for
  1000 stations × 2 faces × 10 combos. Slow, but an optimization target, not a
  blocker. *(Addressed 2026-08-29 — see the performance track: 306 → 72 ms.)*

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
# codes/aci_318_19/equations/shear.py — ONLY formulas; cite the clause; floats in
# N/mm/MPa per ADR-0005; tested against the code's own tables.
def V_c(f_c: float, lambda_: float, rho_w: float, b_w: float, d: float, sigma_Nu: float) -> float:
    """Concrete shear strength [N] — ACI 318-19, Table 22.5.5.1(c).

    f_c [MPa], b_w [mm], d [mm], sigma_Nu [MPa].
    """

def s_max(d: float, V_s: float, threshold: float) -> float:
    """Maximum stirrup spacing [mm] — ACI 318-19 §9.7.6.2.2."""
```

```python
# codes/aci_318_19/checker.py — ONLY factors + unit boundary + orchestration +
# result assembly. Converts pint → float once on entry (prebuilt unit objects,
# never strings), wraps float → pint once in the result.
class ACI318_19:
    phi_v = 0.75

    def check_shear(self, section, force, settings) -> ShearCheckResult:
        Vc = eq.V_c(f_c=section.concrete.f_c.to(MPa).magnitude, b_w=section.width.to(mm).magnitude, ...)
        Vs = eq.V_s(...)
        return ShearCheckResult(V_c=Vc * N, phi_V_n=self.phi_v * (Vc + Vs) * N, ...)
```

The discipline that keeps this honest: **the checker does not compute — it calls
equations and assembles the result — and the equations do not convert.** If a
checker method starts containing algebra, that algebra belongs in `equations/`.
This boundary is enforced mechanically, not by review: a test walks the AST of
every `checker.py` and fails on arithmetic between quantities, and another fails
on any pint import inside `equations/`.

Equation signatures are floats in a fixed N/mm/MPa/N·mm system (ADR-0005): the
measured pint overhead (~1000× on the hot chain, see §1.5) rules out Quantities
inside the equations layer, and pint stays in everything users touch — inputs,
results, reports.

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

**Performance target:** 20,000 checks (1000 stations × 2 faces × 10 combinations)
in **under 5 seconds, serial**, measured at the close of Phase 2b against the
§1.5 baseline. With float equations (ADR-0005) the per-check cost becomes boundary
conversion plus dataclass assembly, so the target is comfortable — it exists so we
know when to stop optimizing, and so a regression is a failed exit criterion
rather than an opinion.

### 3.4 Presentation as its own layer (see ADR-0004)

Word reports, plots, styled tables, and pint display formatting are kept as
first-class features — implemented in modules that consume result dataclasses:

- `_initialize_dicts_*` table builders move out of `codes/` into the report layer.
- Section/rebar plotting moves out of `beam.py` into a plotting module (adopting a
  `plotting_context`-style helper to centralize matplotlib boilerplate).
- Units follow concrete-properties more closely than first stated: the equations
  layer is float-only in a consistent base system, and pint lives at the
  boundaries — inputs, result dataclasses, and display (ADR-0005). *Rounding and
  unit choice for display* belong to the presentation layer.

---

## 4. Migration phases

Incremental; the package imports and the full test suite passes after every phase.
Each phase is independently mergeable, and each has an **exit criterion** beyond
green tests — green only says nothing broke, not that the phase achieved its goal.

### Phase 0 — Hygiene (no risk, do first) — **done 2026-08-29**

- ~~Remove `build/` and `docs/build/` from version control; extend
  `.gitignore`.~~ Already clean; no change needed (see §1.4).
- Delete `tests/pytest.ini`; `pyproject.toml` is the single pytest config. ✔
- Add `tests/conftest.py`; move the shared fixtures there. ✔
- Move `tests/modules_testing.py` out of `tests/` → `scripts/`. ✔
- ~~Stop versioning generated `.docx`/`.xlsx`.~~ Already clean.

**Exit met:** one pytest config (`configfile: pyproject.toml`, rootdir at the
repo root); shared fixtures and the `Agg` backend live in `conftest.py`; no
generated artifacts tracked. Suite unchanged at **733 passed, 1 skipped**
before and after; `ruff check`, `ruff format --check` and `mypy mento/` clean.

### Performance track — cheap wins, parallel to any phase

Independent of the architecture work, driven by the §1.5 measurements.

**Rebar selection precompute — done 2026-08-29.** `rebar.py` was 97 % of design
cost. The candidate search now strips units once and runs entirely on floats,
re-applying pint only to the combinations it keeps; the scoring pass dropped
`df.apply(axis=1)`, which was building a Series per candidate row.

| | before | after |
| --- | --- | --- |
| `design_flexure`, wall clock (min of 25) | 306 ms | **72 ms** (4.3×) |
| Python function calls, one design | 2,879,215 | **709,745** (4.1×) |
| of those, inside pint | 1,001,507 | **304,603** (3.3×) |

`check_shear` and `check_flexure` are unchanged at ~12 ms, as expected — they do
not run the search. Suite unchanged at 773 passed, coverage still 100 %.

**Mechanical `.to()` fix — deprioritized, and this is the interesting part.**
The measurement that motivated it counted 2,906 `.to()` calls per
`design_flexure`, which put the fix at roughly 47 % of runtime. After the
precompute there are **192** `.to()` calls left, 111 of them string-based: the
conversions the fix would have optimized are no longer executed at all. At
~54 µs saved per call the whole fix is now worth about **7 ms of 72**. It stays
worth doing opportunistically when touching a file, but it is not a milestone.

The general lesson for the remaining phases: measure the cheap fix *after* the
structural one, not before — removing work beats speeding it up.

**Exit:** partially met. `design_flexure` is 72 ms against a 50 ms target. What
remains is not a hot spot but ~7,600 Quantity constructions spread thin across
`flexure_design.py` and the ACI module — no profile peak is left in mento
itself. Closing the last 22 ms is Phase 1's float equations (ADR-0005), not more
micro-optimization here.

### Phase 1 — Extract pure equations *(in progress)*

Inside `codes/`, split code formulas out of the mutating orchestrators into
`equations/` modules — **plain floats per ADR-0005**, clause in docstring. The
existing orchestrators keep mutating the beam **but now call the pure
functions**, converting at the call site with prebuilt unit objects. Add
parametrized tests against code tables for each extracted equation. No public
behavior changes; this phase only *creates the layer that does not exist yet*.
The boundary gets teeth from day one: the AST test of §3.1 (no algebra in
checkers, no pint in `equations/`) lands with the first extracted module.

**Done (2026-08-29): ACI 318-19 one-way shear.**
`codes/aci_318_19/equations/shear.py` holds nine clause-cited functions
(Eq. 22.5.5.1.3, Table 22.5.5.1, §22.5.5.1.1, §22.5.1.2, §9.6.3.1, §20.2.2.4,
Table 9.6.3.4, §22.5.8.5.3). The shear orchestrator now calls them and
`ACI_318_19_beam.py` no longer imports `math` at all — every inline formula in
the module is gone. `tests/test_architecture_boundaries.py` enforces the
boundary by AST and was verified to fail on a deliberate violation, not just to
pass. Packaging moved to automatic discovery so a new code subpackage cannot be
left out of the wheel; verified by building one and listing its contents.

**Remaining:** ACI flexure, EN 1992-2004 (shear and flexure), the wall module.

**Exit:** zero `self.` and zero pint imports inside `equations/`; every function
cites its clause and has a parametrized test against the code's table; the
boundary test is in CI. *(Met for the shear module; the exit applies per code
module as each is extracted.)*

### Phase 2a — Migrate tests to the public API

Before touching production code, move the tests that read private attributes and
DataFrame cells (~278 references, plus ~40 distinct private attributes read from
tests) to `flexure_design` / `shear_design` and the public result API. The suite
keeps passing against **unchanged** production code — if any expected value has
to change, the migration was not faithful and it shows immediately. This
de-risks the transcription-error failure mode of Phase 2 at zero cost.

**Exit:** no test reads a private beam attribute; suite green against unmodified
production code.

### Phase 2b — Checks return result dataclasses

Only now does production code change. Extend `design_results.py` with
check-result types; make check/design functions build and return them. Because
the tests already assert through the public API (Phase 2a), this phase does not
touch them: if they pass, behavior was preserved. Private beam attributes are
kept synchronized as a deprecated compatibility layer (`DeprecationWarning` on
documented ones). The envelope accumulators are deleted — enveloping becomes a
pure operation over the returned list. Together with Phase 1's float equations,
this is what unblocks mako: immutability for correctness and parallelism, float
equations for speed.

**Exit:** the §3.3 mako loop runs without mutating any element; the §1.5
benchmark is rerun and meets the < 5 s / 20,000-checks target.

### Phase 3 — Extract presentation

Move `_initialize_dicts_*` out of the code modules into the report layer; move
matplotlib and docx code out of `beam.py`. After this phase, `codes/` contains only
equations + checkers, and elements contain only geometry + thin orchestration.

**Exit:** `import mento.codes` imports neither matplotlib nor docx — verified by
a test.

### Phase 4 — Code registry

Formalize the per-code subpackage protocol (structuralcodes-style
`__title__` / `__year__` / `__materials__` metadata + a registry dict) so CIRSOC
and future codes plug in as subpackages without touching elements. Only worth doing
once phases 1–3 exist; a registry over today's structure would formalize the
coupling instead of removing it.

**Exit:** adding a code requires editing no file under `elements`.

### Features in flight during the migration

Open `help wanted` issues (#123, #124, #125, #127) are not frozen, but they
follow the layer that exists when they land: once Phase 1 creates `equations/`
for a code, new verifications for that code write their formulas there from the
start; until then they follow the current style and are listed in the issue as
migration debt. Nothing is written in the old style *silently*.

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
