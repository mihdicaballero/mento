# ShearWall Transverse-Rebar Design (ACI 318-19)

## Context

`_design_shear_ACI_318_19_wall` ([mento/codes/ACI_318_19_wall.py:374](mento/codes/ACI_318_19_wall.py)) is currently a one-line stub that just delegates to the check function and leaves the user to read `_rho_t_req` and call `set_horizontal_rebar` / `set_vertical_rebar` manually.

The user wants a self-contained design routine that mirrors the beam's `_design_shear_ACI_318_19` ([mento/codes/ACI_318_19_beam.py](mento/codes/ACI_318_19_beam.py)) — i.e., the design call should:

1. Compute required reinforcement (already done by the check function).
2. **Select an actual bar diameter + spacing pair** for both directions (horizontal and vertical mesh) that satisfies the ratio requirement and the §11.7.3 spacing caps.
3. Apply the chosen mesh to the wall instance so a subsequent `check_shear` succeeds without further user input.

This matches the beam pattern where `Rebar(self).transverse_rebar(...)` ([mento/beam.py:634-637](mento/beam.py)) picks the stirrup bar/spacing from `A_v_req` and `V_s_req`.

---

## Scope

Wall-side equivalent for the **distributed mesh** in walls (not stirrups). Phase 0 still: shear-only design, no flexure. No interaction with `Rebar` class — wall mesh is simpler (one bar size + one spacing per direction, no legs).

---

## What's already in place (reused, not rewritten)

- `_check_shear_ACI_318_19_wall` ([ACI_318_19_wall.py:307](mento/codes/ACI_318_19_wall.py:307)) — already computes:
  - `self._Acv`, `self._f_yt_wall`, `self._alpha_c`, `self._hw_lw`
  - `self._V_c_wall`, `self._V_n_wall`, `self._V_n_max`, `self._phi_V_n_wall`, `self._phi_V_n_max_wall`
  - `self._s_h_max`, `self._s_v_max` (via `_calculate_spacing_limits_wall`)
  - `self._rho_t_min` (= 0.0025), `self._rho_l_min` (per §11.6.2, depends on `_rho_t_req`)
  - `self._rho_t_req` (required ρt for strength, floored at ρt,min)
- `ShearWall.set_horizontal_rebar(d_b, s)` / `set_vertical_rebar(d_b, s)` ([shear_wall.py:160-178](mento/shear_wall.py:160)) — apply the chosen mesh and recompute `_rho_t` / `_rho_l`.
- `BeamSettings.minimum_longitudinal_diameter` ([settings.py](mento/settings.py)) — already exposes the smallest standard bar for the active unit system (8 mm metric / 3/8 in imperial).
- `Concrete.unit_system` for metric/imperial branching.

---

## Implementation

### Reference: CalcPad design routine

For reference, the user's CalcPad calc follows ACI §11.5.4/§11.6/§11.7.3 with these steps:

```
V_s,req     = V_u/φ_v - V_c
ρ_t,req     = max( V_s,req/(f_yt·A_cv) ,  ρ_t,min )
A_vh,req    = ρ_t,req · t_w               # area per unit height

s_h,max     = min(l_w/5, 3·t_w, 450 mm)

# CalcPad bar pick (we DIVERGE here — see below):
s_h,design  = s_h,max                     # start from max spacing
A_b,req     = ρ_t,req · t_w · s_h,design / n_c
d_b,req     = sqrt(4·A_b,req/π)
d_b,sel     = first bar in [8,10,12,16,20,25] mm with d_b ≥ d_b,req
s_h,sel     = floor( n_c·A_b,sel/(ρ_t,req·t_w) / 10mm ) · 10mm

# Vertical (§11.6.2)
r_hw        = clamp(hw/lw, 0.5, 2.5)
ρ_l,eq      = 0.0025 + 0.5·(2.5 − r_hw)·(ρ_t,req − 0.0025)
ρ_l,req     = max(0.0025, ρ_l,eq)
# Same bar-pick procedure for the vertical mesh
```

The CalcPad approach starts from `s = s_max` so it picks the **largest acceptable bar**. The user wants the **opposite**: prefer the **smallest** standard bar that can still meet the ratio with a practical spacing. We keep the rest of the structure (formulas, 10 mm rounding, `n_c`) and invert the search.

### File 1: `mento/codes/ACI_318_19_wall.py` — new bar-selection helper

```python
from mento.units import cm, inch
import math

def _select_wall_mesh(
    self: "ShearWall",
    rho_req: float,
    s_max: Quantity,
    n_curtains: int = 1,
) -> tuple[Quantity, Quantity]:
    """
    Pick the SMALLEST standard bar diameter d_b and the largest practical
    spacing s such that  ρ_provided = n_c · A_b / (t · s)  ≥  rho_req  and
    s ≤ s_max.

    Mirrors the CalcPad reference formulas but inverts the bar search so
    Ø8 mm is preferred over Ø10 / Ø12 / … whenever it can still meet the
    ratio within the spacing limits.

    Spacing is rounded DOWN to 10 mm increments (metric) / 1 in (imperial),
    consistent with the CalcPad routine. A practical floor of 5 cm
    (metric) / 2 in (imperial) is enforced so the algorithm never returns
    absurdly tight spacing — when even the smallest bar at the floor can't
    meet rho_req, we move up to the next standard size.

    Raises ValueError if no candidate in the bar list can meet rho_req
    within s_max even at the largest size.
    """
    t = self.thickness
    metric = self.concrete.unit_system == "metric"
    if metric:
        bar_candidates = [8*mm, 10*mm, 12*mm, 16*mm, 20*mm, 25*mm]
        step = 10 * mm
        s_floor = 5 * cm
    else:
        bar_candidates = [3/8*inch, 1/2*inch, 5/8*inch, 3/4*inch, 7/8*inch, 1.0*inch]
        step = 1 * inch
        s_floor = 2 * inch

    s_cap = min(s_max, ...)  # see below; algorithm caps at s_max
    for d_b in bar_candidates:
        A_b = math.pi / 4 * d_b**2
        # Largest spacing that meets the ratio for this bar size
        s_for_rho = (n_curtains * A_b / (t * rho_req)).to(step.units)
        # Round DOWN to the practical step (matches CalcPad's floor(.../10mm))
        n_steps = math.floor((s_for_rho / step).to("").magnitude)
        s = n_steps * step
        # Cap at code §11.7.3 limit
        if s > s_max:
            s = math.floor((s_max / step).to("").magnitude) * step
        # Honour the practical floor — if this bar would require tighter
        # spacing than the floor, try the next (larger) candidate
        if s < s_floor:
            continue
        return d_b, s
    raise ValueError(
        f"No standard bar in {bar_candidates!r} can deliver ρ ≥ {rho_req:.5f} "
        f"within spacing ≤ {s_max:~P}"
    )
```

Notes:
- **n_curtains default = 1**, preserving existing test/setter semantics (`set_horizontal_rebar` already computes `ρ = A_b/(t·s)`). For double-curtain walls the user can pass `n_curtains=2` to the design routine; full plumbing into the setters is a Phase 1 follow-up.
- **10 mm rounding** matches CalcPad and is finer than the previous 2.5 cm draft.
- **5 cm floor** stops the algorithm from returning Ø8 @ 1 cm for very low ρ_req cases (which would be mathematically valid but practically absurd).

### File 1 (cont'd): `mento/codes/ACI_318_19_wall.py` — rewrite `_design_shear_ACI_318_19_wall`

Currently a stub. Replace with a real design step that:

```python
def _design_shear_ACI_318_19_wall(self: "ShearWall", force: Forces) -> None:
    """
    Wall transverse + vertical rebar design per ACI 318-19 Section 11.

    Mirrors the beam's `_design_shear_ACI_318_19`:
      1. Run the check so all wall-state attributes (_rho_t_req, _rho_l_min,
         _s_h_max, _s_v_max, capacities) are populated for this force.
      2. Select a standard bar + practical spacing for the horizontal mesh
         such that ρt ≥ _rho_t_req and s_h ≤ _s_h_max.
      3. Select the same for the vertical mesh against _rho_l_min.
      4. Apply via set_horizontal_rebar / set_vertical_rebar (which also
         updates _rho_t / _rho_l).
    """
    # 1. Populate required ratios and limits
    _check_shear_ACI_318_19_wall(self, force)

    # 2. Horizontal mesh (resists in-plane shear)
    rho_t_req = float(self._rho_t_req.to("").magnitude)
    d_b_h, s_h = _select_wall_mesh(self, rho_t_req, self._s_h_max)

    # 3. Vertical mesh (minimum vertical reinforcement, §11.6.2)
    # NOTE: _rho_l_min was computed by the check using the current _rho_t_req,
    # which is correct for the design step too.
    rho_l_req = float(self._rho_l_min.to("").magnitude)
    d_b_v, s_v = _select_wall_mesh(self, rho_l_req, self._s_v_max)

    # 4. Apply mesh (updates _d_b_h/_s_h/_rho_t and _d_b_v/_s_v/_rho_l)
    self.set_horizontal_rebar(d_b_h, s_h)
    self.set_vertical_rebar(d_b_v, s_v)
```

### File 2: `mento/shear_wall.py` — tighten `design_shear`

Currently iterates forces, tracks the worst-case `_rho_t_req`, and re-runs `check_shear`. After the above change, the per-force `_design_shear_ACI_318_19_wall` call already picks and applies a mesh; we want the FINAL mesh to satisfy the worst-case combination across all forces. Approach:

1. First pass: run `_check_shear_ACI_318_19_wall` for every force, track the worst-case `_rho_t_req` and the corresponding `_rho_l_min`.
2. Single design step at worst case: call `_design_shear_ACI_318_19_wall` with a synthetic worst-case force OR (cleaner) just call `_select_wall_mesh` directly with the tracked worst-case ratios.
3. Apply the mesh once.
4. Re-run `check_shear(forces)` so the returned DataFrame and detail dicts reflect the chosen mesh for all forces.

Replacement body:

```python
def design_shear(self, forces: list[Forces]) -> DataFrame:
    """Design the horizontal and vertical mesh for the worst-case force."""
    if not forces:
        raise ValueError("design_shear requires at least one Forces object.")

    code = self.concrete.design_code
    if code not in ("ACI 318-19", "CIRSOC 201-25"):
        raise NotImplementedError(
            f"Shear wall design not implemented for design code: {code}"
        )

    # Pass 1: find the worst-case required ratios across all forces
    max_rho_t_req: float = 0.0
    max_rho_l_min: float = 0.0
    for force in forces:
        _check_shear_ACI_318_19_wall(self, force)
        rho_t_req = float(self._rho_t_req.to("").magnitude)
        rho_l_min = float(self._rho_l_min.to("").magnitude)
        if rho_t_req > max_rho_t_req:
            max_rho_t_req = rho_t_req
        if rho_l_min > max_rho_l_min:
            max_rho_l_min = rho_l_min

    # Pass 2: pick mesh for the worst case
    d_b_h, s_h = _select_wall_mesh(self, max_rho_t_req, self._s_h_max)
    d_b_v, s_v = _select_wall_mesh(self, max_rho_l_min, self._s_v_max)
    self.set_horizontal_rebar(d_b_h, s_h)
    self.set_vertical_rebar(d_b_v, s_v)

    # Pass 3: re-check for all forces with the designed mesh
    return self.check_shear(forces)
```

Add `_select_wall_mesh` to the existing import from `mento.codes.ACI_318_19_wall`.

---

## Critical files

- `mento/codes/ACI_318_19_wall.py` — add `_select_wall_mesh`, replace `_design_shear_ACI_318_19_wall` body.
- `mento/shear_wall.py` — rewrite `design_shear` to use the new helper and update imports.

## Files unchanged but worth re-checking

- `tests/test_shear_wall.py` — the existing `TestDesignShear` cases (`test_design_returns_dataframe`, `test_rho_t_req_at_least_rho_t_min`) will keep passing. Add one new test that verifies the wall has assigned mesh after `design_shear` (i.e. `_s_h > 0` and `_rho_t >= _rho_t_req`).
- `docs/source/examples/shear_wall_check_ACI_318-19.ipynb` — no required changes; the design cell will now actually pick bars instead of leaving `_rho_t_req` for the user. Optionally extend the notebook to show `node_1.design_shear()` followed by `node_1.check_shear()`.

---

## Verification

1. **Targeted test** — new pytest case in `tests/test_shear_wall.py::TestDesignShear`:
   ```python
   def test_design_assigns_mesh_satisfying_requirements(self, wall_metric):
       wall_metric.design_shear([Forces(V_z=1200 * kN)])
       assert wall_metric._s_h.magnitude > 0
       assert wall_metric._s_v.magnitude > 0
       assert wall_metric._rho_t >= wall_metric._rho_t_req
       assert wall_metric._rho_l >= wall_metric._rho_l_min
       assert wall_metric._s_h <= wall_metric._s_h_max
       assert wall_metric._s_v <= wall_metric._s_v_max
   ```

2. **Full suite**:
   ```
   python -m pytest tests/ --override-ini="addopts=" -q
   ```
   Expect: 434 + 1 = 435 passing.

3. **Manual smoke check** (paste in REPL):
   ```python
   from mento import ShearWall, Concrete_ACI_318_19, SteelBar, Forces, Node
   from mento import MPa, cm, mm, m, kN
   concrete = Concrete_ACI_318_19(name="C25", f_c=25*MPa)
   steel    = SteelBar(name="ADN420", f_y=420*MPa)
   wall = ShearWall(label="W1", concrete=concrete, steel_bar=steel,
                    thickness=25*cm, length=4*m, height=3.5*m, c_c=2*cm)
   node = Node(section=wall, forces=[Forces(V_z=1200*kN)])
   node.design_shear()
   print(wall._d_b_h, wall._s_h, "→ ρt =", wall._rho_t)
   print(wall._d_b_v, wall._s_v, "→ ρl =", wall._rho_l)
   node.check_shear()   # DCR should be < 1.0
   ```

   **Expected**: for the 25 cm × 4 m wall with Vu = 1 200 kN, hw/lw = 0.875 → α_c = 0.25 → φVc = 937.5 kN. Vu > φVc by 262.5 kN, so ρt,req ≈ (1 200/0.75 − 1 250)/(420 · 1.0 m²·MPa) ≈ 0.000833 → clamped to ρt,min = 0.0025. With the small-bar-first algorithm: Ø8 mm needs s ≤ Ab/(t·ρ) = 50.27/(250·0.0025) ≈ 80 mm → rounded down to 80 mm → returned. Result: **Ø8 / 8 cm** for both directions, not Ø12 or Ø20 that the CalcPad approach would pick.

4. **Boundary case** — call `design_shear` with a `Vu` so high that no candidate bar fits within `s_h_max`. Confirm `_select_wall_mesh` raises a clear `ValueError`.


revise:

lets go up to 2.5 cm steps but always rounded to below integer (like 12.5 to 12 cm)
if the case of 8 mm gives spacing of 1 cm, dont show a solution with 8 mm, as an exmaple.
the lowest rebar should be picked after sorting for minimum rebar ratio. i dont want Ø8/5 cm when Ø10/15 is cheaper because or very low rebar ratio. you can consider some wiggle as a functional that considers rebar ratio and spacing, with 80% points for ratio and 20% for spacing, if that makes sense.
