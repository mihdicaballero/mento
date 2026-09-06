# Punching shear — module roadmap

Status: **accepted** (2026-09-06); Phase 1 done the same day. Supersedes the eight-phase
plan drafted in August 2026, which was written against the pre-migration architecture and
is no longer buildable as written.

Companion documents: [`ROADMAP.md`](ROADMAP.md) (Phases 0–4, all done) and
[`adr/`](adr/) — in particular [ADR-0001](adr/0001-immutable-result-dataclasses.md),
[ADR-0002](adr/0002-hybrid-code-layer.md),
[ADR-0004](adr/0004-presentation-as-separate-layer.md) and
[ADR-0005](adr/0005-float-equations-pint-at-the-boundary.md).

---

## 0. Why the August plan had to be rewritten

The original plan put the ACI punching code in `mento/codes/ACI_318_19_punching.py` as
functions typed `self: PunchingNode` that write private attributes back onto the node,
returning a `DataFrame` from `check()`. That was the correct shape *in August*. It is
now the shape the architecture migration removed:

| August plan said | The repository now enforces |
| --- | --- |
| `codes/ACI_318_19_punching.py`, pint end to end | `codes/aci_318_19/equations/punching.py`, **floats only** — `test_equations_do_not_import_units` fails otherwise (ADR-0005) |
| `check()` returns a `DataFrame` | Checks return **frozen result dataclasses**; DataFrames are built in `reports/` — `test_code_modules_do_not_build_report_tables` (ADR-0001, ADR-0004) |
| Each code file mutates the node | A checker converts on entry, calls float equations, returns a state (ADR-0002) |
| Plot method on the element | `mento/plots/` — `test_elements_do_not_import_presentation_libraries` |
| `if concrete.design_code == "ACI 318-19"` | The **registry** — `test_no_module_outside_codes_names_a_design_code` |
| Phase 0: add `M_x` to `Forces` | **Done.** `Forces` already carries `M_x` |

So the phase *content* survives; the file layout and the return types do not.

---

## 1. What exists today

Geometry, reinforcement and containers — no design-code calculation anywhere yet.

- `mento/column.py` — `Column(shape, position, b, h, edge_distance_x, edge_distance_y)`,
  with validation that edge/corner columns declare their edge distances.
- `mento/punching.py` — `PunchingSlab`, `Opening`, `Capital`, `PunchingNode`.
  `PunchingNode.check()` and `.design()` raise `NotImplementedError`.
- `mento/plots/punching.py` — `plot_punching_node`, a plan view of column, capital,
  openings and free edges.
- `mento/reports/punching.py` — the Markdown views `slab_data` / `node_data`
  (§3), added by Phase 1.
- All five names exported lazily from `mento/__init__.py`.
- `Forces` carries `M_x` alongside `M_y`, `V_z`, `N_x`.
- **`shapely` is not a dependency yet.**

Since Phase 1 (§5), `PunchingSlab` also carries the top reinforcement and derives `d` and
ρ from it — §2 is the record of why.

---

## 2. The open question: reinforcement, ρ and *d*

This is the part of the module that is currently wrong, and it has to be settled before
any calculation is written, because every capacity equation reads from it.

### 2.1 What the API did before Phase 1

```python
slab = PunchingSlab(concrete=conc, steel_bar=steel,
                    h=25*cm, c_c=25*mm, rho_x=0.012, rho_y=0.012)
# d_avg = h - c_c - 16 mm  (metric)  |  h - c_c - 5/8 in (imperial)
slab.d_avg = (d_x + d_y) / 2          # "override after construction if needed"
```

Three things were wrong with it:

1. **ρ and *d* are supplied independently, but they are not independent.**
   ρ = A_s /(b·d). A user who overrides `d_avg` has silently changed the denominator of
   a ρ they passed in by hand against a different *d*. Nothing detects it.
2. **The 16 mm is a guess about a bar the user knows.** It is not a modelling
   simplification; it is a placeholder for an input the API refuses to accept.
3. **It cannot express what EN 1992 actually asks for.** EN needs ρ in *each*
   direction, measured over a defined width around the column, and it needs *d* in each
   direction because the two mats sit at different depths. One `rho_x`, one `rho_y` and
   one `d_avg` cannot carry "base mat Ø12@15 plus Ø16@15 extra over the column, x
   outside y".

### 2.2 What each code actually needs

**EN 1992-1-1:2004 §6.4.4** — ρ enters the concrete resistance directly:

```
v_Rd,c = C_Rd,c · k · (100 · ρ_l · f_ck)^(1/3) + k_1·σ_cp   ≥   v_min + k_1·σ_cp
C_Rd,c = 0.18/γ_c        k = 1 + √(200/d) ≤ 2.0   [d in mm]
ρ_l    = √(ρ_ly · ρ_lz) ≤ 0.02
v_min  = 0.035 · k^(3/2) · f_ck^(1/2)
```

with, per §6.4.4(1), ρ_ly and ρ_lz the **bonded tension reinforcement ratios in the two
orthogonal directions, averaged over a slab width equal to the column width plus 3·d on
each side**, and per eq. (6.32) `d = (d_y + d_z)/2`. So EN wants: two ratios, two
effective depths, and a width over which to average — which is exactly where the extra
top bars over the column ("refuerzo") change the answer, since they are concentrated in
that band and the base mat alone underestimates ρ_l.

**ACI 318-19 §22.6.5** — ρ does **not** appear in `v_c`:

```
v_c = least of   0.33·λ_s·λ·√f'c
                 0.17·(1 + 2/β)·λ_s·λ·√f'c
                 0.083·(2 + α_s·d/b_o)·λ_s·λ·√f'c      [MPa]
λ_s = √(2/(1+0.004·d)) ≤ 1.0        [d in mm, §22.5.5.1.3]
```

ACI needs only `d` (and `b_o`, `β`, `α_s`). ρ becomes relevant to ACI only if we later
implement the flexural half of the unbalanced-moment transfer, §8.4.2.3, which
concentrates γ_f·M_sc into the strip `c_2 + 3h` — the same "extra bars over the column"
the EN width is about. Worth stating plainly in the docs, because a user who fills in
`rho_x`/`rho_y` for an ACI check gets no effect from them whatsoever and has no way to
know that.

### 2.3 Proposal (implemented in Phase 1)

**Reinforcement is declared as bars, not as ratios. ρ and *d* are derived and read-only.**
Naming mirrors `OneWaySlab.set_slab_longitudinal_rebar_top`, where position 1 is the base
mat and position 3 is a second, interleaved set — i.e. armadura base y refuerzo:

```python
slab = PunchingSlab(concrete=conc, steel_bar=steel, h=25*cm, c_c=25*mm)

slab.set_rebar_x(d_b1=12*mm, s_b1=15*cm,      # base mat, x
                 d_b3=16*mm, s_b3=15*cm)      # refuerzo over the column, x
slab.set_rebar_y(d_b1=12*mm, s_b1=15*cm)      # base mat only, y

slab.d_x, slab.d_y, slab.d_avg                # derived
slab.rho_x, slab.rho_y, slab.rho_l            # derived
```

Derivation, with `outer_direction: Literal["x","y"] = "x"` on the slab (which mat is
closest to the top face — a real detailing choice that changes *d* by one bar diameter):

```
d_x = h - c_c - d_b,x/2                     (x outer)
d_y = h - c_c - d_b,x - d_b,y/2
d_avg = (d_x + d_y)/2                        EN eq. (6.32); also ACI's d
ρ_x = A_s,x / (b_ρ · d_x),  ρ_y likewise,   capped at 0.02 for EN
ρ_l = √(ρ_x · ρ_y)                           EN §6.4.4(1)
```

where `d_b,x` is the governing (largest) diameter present in x, `A_s,x` is the base mat
plus the refuerzo per unit width, and `b_ρ` is the EN averaging width
`c + 3d` each side — so the refuerzo counts only if it actually extends across that band.

Three consequences worth being explicit about:

- **`d_avg` stops being a settable field.** The escape hatch becomes
  `slab.set_effective_depth(d_x=..., d_y=...)`, which sets both depths *and* re-derives
  ρ from the same bars. Assigning to `slab.d_avg` is what produced the confusion in the
  first place; it should raise, with a message naming the replacement.
- **The 16 mm default survives, demoted.** With no `set_rebar_*` call and no explicit
  depth, `d_avg = h - c_c - 16 mm` (metric) / `- 5/8 in` (imperial), exactly as today —
  it is a reasonable first pass for an ACI check, which needs no ρ. An EN check with no
  ρ declared should **raise**, not silently use ρ = 0, because `(100·ρ·f_ck)^(1/3) = 0`
  drops `v_Rd,c` to `v_min` and quietly under-reports capacity.
- **`rho_x` / `rho_y` are removed as constructor arguments** (decision D2 below). A ρ
  supplied next to a *d* derived from something else is the inconsistency this whole
  section is about; the module is pre-release, so the second way of saying it goes rather
  than being carried as a legacy path. `set_effective_depth` covers the user who has a
  depth from elsewhere, and re-derives ρ against it.

### 2.4 Decisions (settled 2026-09-06)

| # | Question | Decision |
| --- | --- | --- |
| D1 | Is `outer_direction="x"` the right default? | Yes, x outermost — overridable per slab. |
| D2 | Keep `rho_*` constructor args as an override, or remove them? | **Removed.** ρ only ever comes from declared bars. |
| D3 | ρ averaging width: EN's `c + 3d` each side, or the full tributary width? | EN's `c + 3d`. Under D4 it cancels out of the ratio, so it costs nothing to state. |
| D4 | Does the refuerzo need an extent (`l_x`, `l_y` from the column face)? | Assumed to span the ρ band in Phase 1; add the extent in Phase 3 with the rest of the perimeter geometry. |
| D5 | ACI: implement §8.4.2.3 flexural moment transfer? | Shear only in Phase 2; the flexural transfer is its own phase. |

---

## 2bis. Which force components a punching node reads

Settled 2026-09-06, while writing the worked notebook — the geometry was right and the
forces were not. Revisited the same day, and the second reading is the one that stands.

`docs/source/user_guide/local_axes.rst` defines the convention for a *member*: local x is
the element's longitudinal axis, `N_x` is the axial force on it, `V_z` a shear across its
cross-section, `M_y` a bending moment about the section's y-axis.

A punching node is not a member. It is the **column-to-slab connection**:

| Component | At a punching node |
| --- | --- |
| `V_z` | **The design punching load** — the vertical force transferred at the connection. Positive magnitude. |
| `M_x`, `M_y` | **The unbalanced moments** transferred to the slab, about its two in-plane axes. Both read; biaxial transfer is normal at a corner column. |
| `N_x` | **Not used.** |

### Why `V_z` and not `N_x`

The first attempt used `N_x`, reasoning that the punching load *is* the column's axial
force and should be named as the normal force it is. That was wrong for two reasons, and
the second is the one that settles it:

1. **The codes name this quantity `V_u` / `V_Ed`.** ACI 318-19 §22.6 and EN 1992-1-1 §6.4
   both call the punching demand a shear, and an engineer reading the clause alongside the
   API should find the same symbol on both. `V_z` maps to it; `N_x` does not.
2. **Equating it to the column axial load is a modelling simplification, not a
   definition.** Both codes let the load acting *inside* the control perimeter be deducted
   — EN writes it `V_Ed,red`, §6.4.3(3) — and at an edge or corner column the two are not
   the same number anyway. Baking the simplification into the input name would make mento
   assert something the codes do not.

So mento takes `V_z` as *the design punching load already worked out*, however it was
worked out. Where it came from is the engineer's modelling step, not the API's business.

One naming collision to state rather than discover later: in the member convention, `M_x`
about a longitudinal axis is **torsion**. At a punching node there is no longitudinal
axis, and `M_x` is an in-plane moment on the slab. The `Forces` container is shared; the
meaning is per element.

A `Forces` carrying only moments is **refused with a message naming `V_z`**, not silently
checked at zero load — the failure mode otherwise is a DCR of 0.00 that looks like a pass.
That guard is in place; see §4.

---

## 3. `slab.data` and the presentation surface

The request is the same affordance `RectangularBeam` and `ShearWall` already have:
a Markdown view of the object's inputs. It goes in `mento/reports/punching.py`
(new module, matching `reports/walls.py`), and the elements delegate — never import
IPython or matplotlib themselves (ADR-0004, enforced by
`test_elements_do_not_import_presentation_libraries`).

```python
slab.data          # h, c_c, d_x, d_y, d_avg, ρ_x, ρ_y, ρ_l, rebar strings, materials
node.data          # slab.data + column + capital + openings + the forces table
node.results       # data + the punching check summary, once check() has run
```

`node.data` is the one that is actually useful at the REPL: a `PunchingSlab` on its own
does not know its column, and the whole point of reading the data back is checking the
geometry you just typed.

---

## 4. Where the code goes

```
mento/
├── punching.py                       PunchingSlab, Opening, Capital, PunchingNode
│                                       — geometry + orchestration only
├── column.py                         Column (unchanged)
├── punching_results.py               frozen PunchingCheck + envelope_punching
├── plots/punching.py                 plan view (exists; grows the critical perimeter)
├── reports/punching.py               data views (exist); results views with Phase 2
└── codes/
    ├── registry.py                   check_punching / design_punching hooks
    ├── ACI_318_19_punching.py        the ACI checker
    ├── EN_1992_2004_punching.py      the EN checker
    ├── aci_318_19/
    │   ├── code.py                   registers check_punching
    │   └── equations/punching.py     floats: v_c, λ_s, α_s, γ_v, J_c …
    └── en_1992_2004/
        ├── code.py                   registers check_punching
        └── equations/punching.py     floats: v_Rd,c, k, ρ_l, β, v_Rd,max …
```

**The skeleton of all of that is in place as of 2026-09-06** — every file above
exists, `node.check()` dispatches through the registry, and the preconditions each
code needs are written and tested. What is deliberately absent is the arithmetic:
every function in the two `equations/punching.py` modules is a signature with its
clause, its argument units and a body that raises, because the formulas are being
recreated from a validated Calcpad sheet and each one lands with the worked example
that checks it. A test asserts no stub has quietly grown a body that returns a
number; it shrinks as the module fills in.

Preconditions already enforced by the checkers, so Phase 2 is only equations:

- **A force with no `V_z` is refused, naming `V_z`** (§2bis). Checking at zero load
  would report a DCR of 0.00, which reads as a pass.
- **A capital or an opening is refused** until Phase 3. Ignoring either is not
  conservative — it is wrong in the unsafe direction for an opening and the safe one
  for a capital, with nothing on the result to say which.
- **EN refuses a slab with no declared ρ**, naming `set_rebar_x()` / `set_rebar_y()`.
  ACI does not, and must not: its `v_c` never reads ρ.

`PunchingCheck` carries `label`, `b_0`, `d`, `v_u`, `v_c` and `DCR` — the intersection
both codes report — with `v_c` being the resistance its own `DCR` was formed from, the
same contract `ShearCheck` keeps. Its envelope *is* one of the combinations rather than
a mixture of several: punching is a single stress against a single resistance on one
perimeter, unlike flexure, where each face envelopes its quantities independently.

Two boundary consequences to plan for:

- **`shapely` enters at Phase 3, not Phase 2.** The perimeter geometry for an
  unobstructed column is closed-form; shapely earns its place when openings and capitals
  start clipping the perimeter. It must be imported from `punching.py` (geometry), not
  from `codes/*/equations/` (floats), so the perimeter arrives at the equations as a
  length and a section modulus, already computed.
- **The registry gains optional punching hooks**, the way it already carries optional
  wall hooks (`check_shear_wall = None`), so a code with no punching implementation says
  so by name instead of failing obscurely.

---

## 5. Phases

Each phase is a commit on `feat/punching-general`, and each has an exit criterion that is
a passing test, not a judgement.

### Phase 1 — Reinforcement, ρ and *d* (the §2.3 proposal) — **done 2026-09-06**

No design code involved. `set_rebar_x` / `set_rebar_y`, derived `d_x`/`d_y`/`d_avg`,
derived `ρ_x`/`ρ_y`/`ρ_l`/`A_s_x`/`A_s_y`, `outer_direction`, `set_effective_depth`, the
raise on `slab.d_avg = ...`, `rho_x`/`rho_y` gone from the constructor, and
`slab.data` / `node.data` in a new `reports/punching.py`.

*Exit met.* A slab built from bars reproduces hand-computed `d_x`, `d_y`, ρ_x, ρ_y in
metric and imperial; a slab with no bars still gives the historical
`d_avg = h - c_c - 16 mm` and returns `None` rather than zero for ρ; `node.data` renders
with `punching.py` importing no IPython, so `test_architecture_boundaries.py` stays green.
`mento/punching.py` and `mento/reports/punching.py` are at 100 % line coverage.

Two things Phase 2 inherits from it:

- `PunchingSlab.has_rebar` is the flag the EN checker reads to refuse a check with no ρ.
  ACI does not need it: `v_c` never uses ρ.
- The largest declared bar in a direction governs that mat's depth. Two interleaved sets
  of different diameters do not physically sit at one depth; the larger is the reading a
  drawing dimensions to, and the conservative one.

### Phase 2 — ACI 318-19 check, no capital, no openings

**Skeleton landed 2026-09-06** (see §4): `PunchingCheck`, both checkers, the registry
hooks, `node.check()`, the preconditions, and both `equations/punching.py` modules as
clause-cited signatures that raise.

**What is left is the arithmetic**, recreated from the validated Calcpad sheet:

1. The critical section — `b_0` and the section property `J_c` — for the three column
   positions, rectangular and circular. This is geometry, so it lives in `punching.py`
   and reaches the equations as numbers; only the `d/2` offset itself is a clause.
2. `v_c` from §22.6.5.2, least of the three expressions, with `λ_s` (§22.5.5.1.3) and
   `α_s` (§22.6.5.3).
3. `γ_v` from §8.4.4.2.2 and the combined stress at the critical point from §8.4.4.2.3,
   biaxial.
4. The `PunchingCheck` assembled from them, and the report tables in `reports/punching.py`.

*Exit:* a CRSI / ACI worked example within the suite's usual tolerance; a biaxial case;
`test_architecture_boundaries.py` green — which means the equations import no units and
the code layer builds no tables.

### Phase 3 — Capital and openings (shapely enters)

Capital first — two perimeters, at d/2 from the capital edge and at d/2 from the column
face with the reduced *d*, governing one controls — then openings: the proximity rule
(ACI ≤ 10h, EN ≤ 6d), tangents from the column centroid, arc subtraction, and the shifted
centroid a clipped perimeter produces. That centroid shift is the piece of this module
most likely to be wrong silently, so it gets its own tests, separate from any check.

*Exit:* a distant opening changes nothing; a large near opening reduces `b_o` by a
hand-verified amount; the modified-centroid function is tested on its own.

### Phase 4 — Shear reinforcement design

Closed stirrups and headed studs; ACI §22.6.6 / §22.6.7 upper limits, `A_v` per
perimeter, iteration outward until an unreinforced perimeter passes, `s_r ≤ 0.5d`,
`s_t ≤ 2d`. Returns a `PunchingDesign` frozen result.

### Phase 5 — EN 1992-1-1 check and design

Skeleton landed with Phase 2's: the checker, its ρ precondition, and the equation
signatures. Left: `v_Rd,c` (eq. 6.47) with the ρ of Phase 1, `u_1` at 2d, β per position
(eqs. 6.38–6.46), `v_Rd,max` at `u_0`, reinforcement per eq. (6.52) with `s_r ≤ 0.75d`,
and `u_out`.

*Exit:* a published EN worked example; and `test_a_new_design_code_needs_no_element_edited`
still passes — i.e. EN was added without touching `punching.py`.

### Phase 6 — Plot, Word report, summary

Critical perimeter and reinforcement perimeters on the existing plan view;
`punching_results_detailed_doc()` via `DocumentBuilder`; `PunchingSummary` mirroring
`BeamSummary` for a batch of nodes.

---

## 6. What this document does not decide

- Whether punching belongs on `OneWaySlab` / a future `TwoWaySlab` rather than on its own
  `PunchingSlab`. Standalone is the current bet; revisit if a two-way slab element lands.
- Post-tensioned slabs (EN §6.4.5 decompression, ACI §22.6.5.5) — out of scope.
- Seismic detailing for slab-column connections (ACI §18.14.5) — out of scope.
