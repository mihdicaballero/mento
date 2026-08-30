from __future__ import annotations

import math
from typing import TYPE_CHECKING

from pint import Quantity

from mento.codes.aci_318_19.equations import shear as shear_eq
from mento.codes.aci_318_19.equations import wall as wall_eq
from mento.codes.check_state import WallShearCheckState, apply_wall_shear_state, new_wall_shear_state
from mento.material import Concrete_ACI_318_19
from mento.units import MPa, cm, mm, psi, inch, dimensionless
from mento.forces import Forces

if TYPE_CHECKING:
    from mento.shear_wall import ShearWall


##########################################################
# WALL MESH BAR CATALOGUE
##########################################################

# Standard reinforcing bar diameters for wall distributed mesh design.
# ACI 318-19 uses the same catalogue for both mesh directions.
_ACI_WALL_BARS_METRIC = [10 * mm, 12 * mm, 16 * mm, 20 * mm, 25 * mm]
_ACI_WALL_BARS_IMPERIAL = [
    0.375 * inch,  # #3
    0.5 * inch,  # #4
    0.625 * inch,  # #5
    0.75 * inch,  # #6
    0.875 * inch,  # #7
    1.0 * inch,  # #8
]
# CIRSOC 201-25 (metric only) reuses the ACI shear provisions but allows
# smaller bars for the design: Ø6 mm minimum for the transverse mesh and
# Ø10 mm minimum for the vertical mesh.
_CIRSOC_WALL_BARS_TRANSVERSE = [6 * mm, 8 * mm, 10 * mm, 12 * mm, 16 * mm, 20 * mm, 25 * mm]
_CIRSOC_WALL_BARS_VERTICAL = [10 * mm, 12 * mm, 16 * mm, 20 * mm, 25 * mm]

# Crack-control cap — the selector scores bars at/below this first; larger
# bars are used only as a fallback when nothing capped can meet the demand.
_WALL_BAR_CAP_METRIC = 12 * mm
_WALL_BAR_CAP_IMPERIAL = 0.5 * inch  # #4


##########################################################
# HELPERS
##########################################################


def _wall_units(self: "ShearWall") -> tuple[Quantity, Quantity]:
    """(stress, length) units for the wall's unit system."""
    return (psi, inch) if self.concrete.is_imperial else (MPa, mm)


def _calculate_f_yt_wall(self: "ShearWall") -> Quantity:
    """Cap fyt at 420 MPa (metric) / 60 ksi (imperial) per ACI 318-19 §20.2.2.4."""
    stress_unit, _ = _wall_units(self)
    return (
        shear_eq.max_yield_strength_for_shear(
            self.steel_bar.f_y.to(stress_unit).magnitude, is_imperial=self.concrete.is_imperial
        )
        * stress_unit
    )


def _calculate_alpha_c(self: "ShearWall", st: WallShearCheckState) -> float:
    """α_c per ACI 318-19 §11.5.4.6. Stores hw/lw as a side effect on st.hw_lw."""
    hw_lw = (self.height / self.length).to("").magnitude
    st.hw_lw = hw_lw
    return wall_eq.alpha_c(hw_lw, is_imperial=self.concrete.is_imperial)


def _calculate_wall_Acv(self: "ShearWall", st: WallShearCheckState) -> None:
    """Acv = lw × t  (ACI 318-19 §11.5.4.6)."""
    st.Acv = self.length * self.thickness


def _calculate_wall_shear_strength(
    self: "ShearWall",
    st: WallShearCheckState,
    concrete: Concrete_ACI_318_19,
) -> None:
    """
    ACI 318-19 §11.5.4.6:
        Vc      = α_c × λ × √f'c × Acv
        Vs      = ρt × fyt × Acv
        Vn      = Vc + Vs
        Vn,max  = 0.66 × λ × √f'c × Acv  (metric)
                  8   × λ × √f'c × Acv  (imperial)
    Acv
    """
    lam = concrete.lambda_factor
    phi_v = concrete.phi_v
    is_imperial = concrete.is_imperial
    stress_unit, _ = _wall_units(self)
    f_c_mag = concrete.f_c.to(stress_unit).magnitude

    Vc = wall_eq.concrete_shear_stress(f_c_mag, st.alpha_c, lam) * stress_unit * st.Acv
    Vn_max = wall_eq.max_shear_stress(f_c_mag, lam, is_imperial=is_imperial) * stress_unit * st.Acv
    Vs = (
        wall_eq.reinforcement_shear_stress(float(self._rho_t), st.f_yt_wall.to(stress_unit).magnitude)
        * stress_unit
        * st.Acv
    )
    Vn = Vc + Vs

    st.V_c_wall = Vc.to("kN")  # type:ignore
    st.V_s_wall = Vs.to("kN")  # type:ignore
    st.V_n_wall = Vn.to("kN")  # type:ignore
    st.V_n_max = Vn_max.to("kN")  # type:ignore
    st.phi_V_n_wall = (phi_v * min(Vn, Vn_max)).to("kN")  # type:ignore
    st.phi_V_n_max_wall = (phi_v * Vn_max).to("kN")  # type:ignore


def _calculate_rho_min_wall(self: "ShearWall", st: WallShearCheckState) -> None:
    """
    ACI 318-19 §11.6.1 / §11.6.2:
        ρt_min = 0.0025 (horizontal, always)
        ρl_min = max(0.0025, ρl_eq)  (vertical)

    Eq. (11.6.2):
        ρl ≥ 0.0025 + 0.5·(2.5 − hw/lw)·(ρt − 0.0025)

    The hw/lw ratio is clamped to [0.5, 2.5]:
      - hw/lw ≥ 2.5 → ρl_eq = 0.0025 (only minimum vertical)
      - hw/lw ≤ 0.5 → ρl_eq = ρt (vertical equals horizontal)

    ρl_req need not exceed ρt required for strength (§11.5.4.3).
    """
    st.rho_t_min = wall_eq.MIN_REINFORCEMENT_RATIO * dimensionless
    st.rho_l_min = wall_eq.min_vertical_reinforcement_ratio(float(st.hw_lw), float(st.rho_t_req)) * dimensionless


def _calculate_spacing_limits_wall(self: "ShearWall", st: WallShearCheckState) -> None:
    """
    ACI 318-19 §11.7.3 spacing limits.
    Horizontal: s_h,max = min(lw/5, 3t, 450 mm / 18 in)
    Vertical:   s_v,max = min(lw/3, 3t, 450 mm / 18 in)
    """
    is_imperial = self.concrete.is_imperial
    _, length_unit = _wall_units(self)
    lw = self.length.to(length_unit).magnitude
    t = self.thickness.to(length_unit).magnitude

    st.s_h_max = wall_eq.max_horizontal_spacing(lw, t, is_imperial=is_imperial) * length_unit
    st.s_v_max = wall_eq.max_vertical_spacing(lw, t, is_imperial=is_imperial) * length_unit


##########################################################
# MAIN CHECK FUNCTION
##########################################################


def _check_shear_ACI_318_19_wall(self: "ShearWall", force: Forces) -> WallShearCheckState:
    """
    ACI 318-19 Section 11 shear check for a structural wall.

    Calculation only: the result is returned as a value, and only the reporting
    path copies it back onto the wall. See the beam's shear check for the same
    split.
    """
    if not isinstance(self.concrete, Concrete_ACI_318_19):
        raise TypeError("ACI 318-19 wall shear check requires Concrete_ACI_318_19.")

    concrete = self.concrete
    st = new_wall_shear_state(self)

    # 1. Demand
    st.V_u = force._V_z.to("kN")
    st.N_u = force._N_x.to("kN")

    # 2. Geometry: Acv = lw × t
    _calculate_wall_Acv(self, st)

    # 3. Material: fyt cap
    st.f_yt_wall = _calculate_f_yt_wall(self)

    # 4. α_c based on hw/lw (also sets st.hw_lw)
    st.alpha_c = _calculate_alpha_c(self, st)

    # 5. Shear strength components
    _calculate_wall_shear_strength(self, st, concrete)

    # 6. Spacing limits
    _calculate_spacing_limits_wall(self, st)

    # 7. Required ρt for design
    phi_v = concrete.phi_v
    f_c = concrete.f_c
    lam = concrete.lambda_factor

    stress_unit, _ = _wall_units(self)
    Vc_intensity = wall_eq.concrete_shear_stress(f_c.to(stress_unit).magnitude, st.alpha_c, lam) * stress_unit

    st.rho_t_min = wall_eq.MIN_REINFORCEMENT_RATIO * dimensionless
    rho_t_req_raw = ((st.V_u / phi_v) / st.Acv - Vc_intensity) / st.f_yt_wall
    rho_t_req_raw = rho_t_req_raw.to("")
    st.rho_t_req = max(rho_t_req_raw, st.rho_t_min)

    # 8. Minimum reinforcement ratios (ρl,min depends on ρt,req per §11.6.2)
    _calculate_rho_min_wall(self, st)

    # 9. DCR
    phi_Vn_eff = min(st.phi_V_n_wall, st.phi_V_n_max_wall)
    if phi_Vn_eff.magnitude == 0:  # pragma: no cover - defensive: Vc > 0 for valid concrete
        st.DCR = float("inf")
    else:
        st.DCR = float((st.V_u / phi_Vn_eff).to("").magnitude)

    return st


##########################################################
# DESIGN — BAR SELECTION
##########################################################


def _select_wall_mesh(
    self: "ShearWall",
    rho_req: float,
    s_max: Quantity,
    bar_list: list,
) -> tuple:
    """
    Pick (d_b, s) for one wall mesh direction.

    Two-tier search:
      1. Apply the 80/20 scoring functional to bars up to the crack-control
         cap (Ø12 mm / #4). Return the best-scoring capped candidate.
      2. Only if no capped bar yields a valid candidate, retry on the full
         ``bar_list`` (up to Ø25).

    Scoring functional (per candidate (d_b, s)):
        rho_provided   = n_curtains · A_b / (t · s)   (mesh on both faces, E.F.)
        ratio_score    = rho_req / rho_provided       ∈ (0, 1]  (minimise steel)
        diameter_score = d_min / d_b                  ∈ (0, 1]  (prefer small bar)
        score          = 0.80 * ratio_score + 0.20 * diameter_score
    `d_min` is the smallest bar in the tier being scored. The functional thus
    prefers the lowest reinforcement ratio (least steel) and, for near-equal
    ratios, the smaller diameter. Ties are broken toward the smaller diameter.

    Spacing grid: 2.5 cm multiples floored to whole cm (metric) / integer
    inches (imperial), with a practical floor of 5 cm / 2 in.
    """
    t = self.thickness
    n_c = self._n_curtains  # mesh on both faces (E.F.)
    metric = self.concrete.unit_system == "metric"
    step = 2.5 * cm if metric else 1.0 * inch
    s_floor = 5.0 * cm if metric else 2.0 * inch
    cap = _WALL_BAR_CAP_METRIC if metric else _WALL_BAR_CAP_IMPERIAL
    unit = 1 * cm if metric else 1 * inch

    # Build the spacing grid up to s_max
    grid: list = []
    k = 2
    while True:
        s = math.floor((k * step).to(unit.units).magnitude) * unit
        if s > s_max:
            break
        if s >= s_floor:
            grid.append(s)
        k += 1

    def _best(candidate_bars: list[Quantity]) -> tuple[tuple[float, float], Quantity, Quantity] | None:
        best = None  # ((score, -d_b_mm), d_b, s)
        d_min = min(candidate_bars)  # smallest bar in this tier
        for d_b in candidate_bars:
            A_b = math.pi / 4 * d_b**2
            feasible = [s for s in grid if (n_c * A_b / (t * s)).to("").magnitude >= rho_req]
            if not feasible:
                continue
            s = min(max(feasible), s_max)
            rho_prov = (n_c * A_b / (t * s)).to("").magnitude
            ratio_score = rho_req / rho_prov
            diameter_score = (d_min / d_b).to("").magnitude
            score = 0.80 * ratio_score + 0.20 * diameter_score
            key = (score, -d_b.to("mm").magnitude)  # tie-break: smaller bar
            if best is None or key > best[0]:
                best = (key, d_b, s)
        return best

    capped = [b for b in bar_list if b <= cap]
    result = _best(capped) or _best(bar_list)  # tier 1, then fallback
    if result is None:
        raise ValueError(
            f"No standard bar can satisfy ρ ≥ {rho_req:.5f} with spacing in [{s_floor:.0f~P}, {s_max:.0f~P}]."
        )
    return result[1], result[2]


##########################################################
# DESIGN FUNCTION
##########################################################


def _design_shear_wall_core(
    self: "ShearWall",
    forces: list,
    transverse_bars: list,
    vertical_bars: list,
) -> None:
    """
    Worst-case wall mesh design across all force combinations.

      1. Run the shear check for every force; track the worst-case ρt,req.
      2. Derive ρl,min for that worst case per ACI 318-19 §11.6.2.
      3. Select the horizontal mesh (ρt,req / s_h,max) and the vertical mesh
         (ρl,min / s_v,max). The vertical mesh is ALWAYS the code minimum —
         no flexure design in Phase 0.
      4. Apply both meshes via set_horizontal_rebar / set_vertical_rebar.
    """
    if not forces:
        raise ValueError("Wall shear design requires at least one Forces object.")

    max_rho_t_req = 0.0
    state = None
    for force in forces:
        state = _check_shear_ACI_318_19_wall(self, force)
        max_rho_t_req = max(max_rho_t_req, state.rho_t_req.to("").magnitude)
    # Designing is meant to change the wall, so the last state is applied.
    assert state is not None  # the empty-forces case raised above
    apply_wall_shear_state(self, state)

    # ρl,min per §11.6.2 using the worst-case ρt,req (geometry-only hw/lw).
    r_hw = max(0.5, min(state.hw_lw, 2.5))
    rho_l_eq = 0.0025 + 0.5 * (2.5 - r_hw) * (max_rho_t_req - 0.0025)
    max_rho_l_min = max(0.0025, rho_l_eq)

    d_b_h, s_h = _select_wall_mesh(self, max_rho_t_req, self._s_h_max, transverse_bars)
    d_b_v, s_v = _select_wall_mesh(self, max_rho_l_min, self._s_v_max, vertical_bars)

    self.set_horizontal_rebar(d_b_h, s_h)
    self.set_vertical_rebar(d_b_v, s_v)


def _design_shear_ACI_318_19_wall(self: "ShearWall", forces: list) -> None:
    """Wall mesh design per ACI 318-19 Section 11.

    Also serves CIRSOC 201-25, which reuses the ACI shear provisions verbatim
    and differs only in the design bar catalogue (Ø6 mm minimum transverse,
    Ø10 mm minimum vertical). The catalogue is selected here by design code —
    mirroring how ``rebar.py`` keeps the CIRSOC beam logic in one place.
    """
    if self.concrete.design_code == "CIRSOC 201-25":
        transverse_bars = _CIRSOC_WALL_BARS_TRANSVERSE
        vertical_bars = _CIRSOC_WALL_BARS_VERTICAL
    elif self.concrete.unit_system == "metric":
        transverse_bars = vertical_bars = _ACI_WALL_BARS_METRIC
    else:
        transverse_bars = vertical_bars = _ACI_WALL_BARS_IMPERIAL
    _design_shear_wall_core(self, forces, transverse_bars=transverse_bars, vertical_bars=vertical_bars)
