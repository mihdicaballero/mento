# Punching shear — module roadmap

Status: **draft, for review** (2026-09-06). Supersedes the eight-phase plan drafted in
August 2026, which was written against the pre-migration architecture and is no longer
buildable as written. Nothing in this document has been implemented yet beyond §1.

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

## 1. What exists today (Phase 1 of the old plan, landed)

Geometry and containers only — no calculation anywhere.

- `mento/column.py` — `Column(shape, position, b, h, edge_distance_x, edge_distance_y)`,
  with validation that edge/corner columns declare their edge distances.
- `mento/punching.py` — `PunchingSlab`, `Opening`, `Capital`, `PunchingNode`.
  `PunchingNode.check()` and `.design()` raise `NotImplementedError`.
- `mento/plots/punching.py` — `plot_punching_node`, a plan view of column, capital,
  openings and free edges.
- `tests/test_punching.py` — 380 lines, geometry and validation only.
- All five names exported lazily from `mento/__init__.py`.
- `Forces` carries `M_x` alongside `M_y`, `V_z`, `N_x`.
- **`shapely` is not a dependency yet.**

---

## 2. The open question: reinforcement, ρ and *d*

This is the part of the module that is currently wrong, and it has to be settled before
any calculation is written, because every capacity equation reads from it.

### 2.1 What the API does today

```python
slab = PunchingSlab(concrete=conc, steel_bar=steel,
                    h=25*cm, c_c=25*mm, rho_x=0.012, rho_y=0.012)
# d_avg = h - c_c - 16 mm  (metric)  |  h - c_c - 5/8 in (imperial)
slab.d_avg = (d_x + d_y) / 2          # "override after construction if needed"
```

Three things are wrong with it:

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
`rho_x`/`rho_y` for an ACI check today gets no effect from them whatsoever and has no
way to know that.

### 2.3 Proposal

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
- **`rho_x` / `rho_y` as constructor arguments stay accepted** as a direct override for
  the user who has ρ from a FE model and no bar schedule, but they are then the *source*
  of ρ and `set_rebar_*` refuses to co-exist with them (raise, don't silently pick one).

### 2.4 Open decisions for you

| # | Question | Default if you don't say |
| --- | --- | --- |
| D1 | Is `outer_direction="x"` the right default? | yes, x outermost |
| D2 | Should `rho_*` constructor args be kept at all, or removed now while the module is pre-release? | kept, as an override |
| D3 | ρ averaging width: EN's `c + 3d` each side, or the full tributary width? | EN's `c + 3d` |
| D4 | Does the refuerzo need an extent (`l_x`, `l_y` from the column face), so we can tell whether it reaches across the ρ band, or do we assume it always does? | assume it always does in Phase 1; add extent in Phase 3 |
| D5 | ACI: implement §8.4.2.3 flexural moment transfer (which is what makes ρ matter to ACI), or leave punching as a shear-stress check only? | shear only in Phase 2; flexural transfer is its own phase |

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
├── punching.py                     PunchingSlab, Opening, Capital, PunchingNode
│                                     — geometry + orchestration only
├── column.py                       Column (unchanged)
├── punching_results.py             frozen dataclasses: PunchingCheck, PunchingDesign
├── plots/punching.py               plan view (exists; grows the critical perimeter)
├── reports/punching.py             data / results / detailed views + Word doc
└── codes/
    ├── registry.py                 + check_punching / design_punching hooks
    ├── aci_318_19/
    │   ├── code.py                 registers the hooks
    │   └── equations/punching.py   floats: v_c, λ_s, γ_v, J_c, b_o …
    └── en_1992_2004/
        ├── code.py
        └── equations/punching.py   floats: v_Rd,c, k, β, u_1, u_out …
```

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

### Phase 1 — Reinforcement, ρ and *d* (the §2.3 proposal)

No design code involved. `set_rebar_x` / `set_rebar_y`, derived `d_x`/`d_y`/`d_avg`,
derived `ρ_x`/`ρ_y`/`ρ_l`, `outer_direction`, `set_effective_depth`, the raise on
`slab.d_avg = ...`, and `slab.data` / `node.data` in `reports/punching.py`.

*Exit:* a slab built from bars reproduces hand-computed `d_x`, `d_y`, ρ_x, ρ_y for a
metric and an imperial case; the legacy `rho_x`/`rho_y` path still gives today's numbers;
`node.data` renders without `punching.py` importing IPython.

### Phase 2 — ACI 318-19 check, no capital, no openings

`equations/punching.py` (floats), the checker in `aci_318_19/code.py`, a `PunchingCheck`
frozen result, registry hooks, `node.check()`. Covers the three column positions,
rectangular and circular, uniaxial and biaxial via γ_v and J_c.

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

`v_Rd,c` with the ρ of Phase 1, `u_1` at 2d, β per position (eqs. 6.38–6.46),
`v_Rd,max` at `u_0`, reinforcement per eq. (6.52) with `s_r ≤ 0.75d`, and `u_out`.

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
