from __future__ import annotations

import math
from typing import TYPE_CHECKING

import pandas as pd
from pint import Quantity

from mento.material import Concrete_ACI_318_19
from mento.units import MPa, cm, mm, ksi, inch, dimensionless
from mento.forces import Forces

if TYPE_CHECKING:
    from mento.shear_wall import ShearWall


##########################################################
# WALL MESH BAR CATALOGUE
##########################################################

# Standard reinforcing bar diameters for wall distributed mesh design.
# ACI 318-19 uses the same catalogue for both mesh directions.
_ACI_WALL_BARS_METRIC = [8 * mm, 10 * mm, 12 * mm, 16 * mm, 20 * mm, 25 * mm]
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


def _calculate_f_yt_wall(self: "ShearWall") -> Quantity:
    """Cap fyt at 420 MPa (metric) / 60 ksi (imperial) per ACI 318-19 §20.2.2.4."""
    if self.concrete.unit_system == "metric":
        return min(self.steel_bar.f_y, 420 * MPa)
    else:
        return min(self.steel_bar.f_y, 60 * ksi)


def _calculate_alpha_c(self: "ShearWall") -> float:
    """
    ACI 318-19 §11.5.4.6 — α_c as a function of hw/lw.
    Metric:   0.25 (hw/lw ≤ 1.5) to 0.17 (hw/lw ≥ 2.0), linear interpolation.
    Imperial: 3.0  (hw/lw ≤ 1.5) to 2.0  (hw/lw ≥ 2.0), linear interpolation.
    Stores hw_lw ratio as side-effect on self._hw_lw.
    """
    hw_lw = (self.height / self.length).to("").magnitude
    self._hw_lw = hw_lw

    if self.concrete.unit_system == "metric":
        alpha_hi, alpha_lo = 0.25, 0.17
    else:
        alpha_hi, alpha_lo = 3.0, 2.0

    if hw_lw <= 1.5:
        return alpha_hi
    elif hw_lw >= 2.0:
        return alpha_lo
    else:
        t = (hw_lw - 1.5) / 0.5
        return alpha_hi + t * (alpha_lo - alpha_hi)


def _calculate_wall_Acv(self: "ShearWall") -> None:
    """Acv = lw × t  (ACI 318-19 §11.5.4.6)."""
    self._Acv = self.length * self.thickness


def _calculate_wall_shear_strength(
    self: "ShearWall",
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
    f_c = concrete.f_c
    lam = concrete.lambda_factor
    phi_v = concrete.phi_v

    from mento.units import psi

    if concrete.unit_system == "metric":
        # take magnitude for sqrt to avoid applying math.sqrt to a Quantity/Unit
        sqrt_fc_term = math.sqrt((f_c / MPa).magnitude) * MPa
        Vc = self._alpha_c * lam * sqrt_fc_term * self._Acv
        Vn_max = 0.66 * lam * sqrt_fc_term * self._Acv
    else:
        # imperial branch: use psi magnitude for sqrt
        sqrt_fc_term = math.sqrt((f_c / psi).magnitude) * psi
        Vc = self._alpha_c * lam * sqrt_fc_term * self._Acv
        Vn_max = 8.0 * lam * sqrt_fc_term * self._Acv

    Vs = self._rho_t * self._f_yt_wall * self._Acv
    Vn = Vc + Vs

    self._V_c_wall = Vc.to("kN")  # type:ignore
    self._V_s_wall = Vs.to("kN")  # type:ignore
    self._V_n_wall = Vn.to("kN")  # type:ignore
    self._V_n_max = Vn_max.to("kN")  # type:ignore
    self._phi_V_n_wall = (phi_v * min(Vn, Vn_max)).to("kN")  # type:ignore
    self._phi_V_n_max_wall = (phi_v * Vn_max).to("kN")  # type:ignore


def _calculate_rho_min_wall(self: "ShearWall") -> None:
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
    self._rho_t_min = 0.0025 * dimensionless

    r_hw = max(0.5, min(float(self._hw_lw), 2.5))
    rho_l_eq = 0.0025 + 0.5 * (2.5 - r_hw) * (float(self._rho_t_req) - 0.0025)
    self._rho_l_min = max(0.0025, rho_l_eq) * dimensionless


def _calculate_spacing_limits_wall(self: "ShearWall") -> None:
    """
    ACI 318-19 §11.7.3 spacing limits.
    Horizontal: s_h,max = min(lw/5, 3t, 450 mm / 18 in)
    Vertical:   s_v,max = min(lw/3, 3t, 450 mm / 18 in)
    """
    lw = self.length
    t = self.thickness
    if self.concrete.unit_system == "metric":
        abs_max = 450 * mm
    else:
        abs_max = 18 * inch

    self._s_h_max = min(lw / 5, 3 * t, abs_max)
    self._s_v_max = min(lw / 3, 3 * t, abs_max)


def _compile_wall_shear_dicts(self: "ShearWall", force: Forces) -> None:
    """Populate result dicts used by detailed output methods."""
    phi_v = self.concrete.phi_v
    unit = "kN" if self.concrete.unit_system == "metric" else "kip"

    # Materials
    self._materials_shear_wall = {
        "Materials": [
            "Section Label",
            "Concrete strength",
            "Steel reinforcement yield strength",
            "Normalweight concrete",
            "Safety factor for shear",
        ],
        "Variable": ["", "fc", "fy", "λ", "Øv"],
        "Value": [
            self.label,
            round(self.concrete.f_c.to("MPa").magnitude, 2)
            if self.concrete.unit_system == "metric"
            else round(self.concrete.f_c.to("psi").magnitude, 0),
            round(self.steel_bar.f_y.to("MPa").magnitude, 2)
            if self.concrete.unit_system == "metric"
            else round(self.steel_bar.f_y.to("ksi").magnitude, 2),
            self.concrete.lambda_factor,
            phi_v,
        ],
        "Unit": [
            "",
            "MPa" if self.concrete.unit_system == "metric" else "psi",
            "MPa" if self.concrete.unit_system == "metric" else "ksi",
            "",
            "",
        ],
    }
    # Geometry
    self._geometry_shear_wall = {
        "Geometry": [
            "Wall thickness",
            "Wall length",
            "Wall height",
            "Aspect ratio",
            "Gross shear area",
        ],
        "Variable": ["t", "lw", "hw", "hw/lw", "Acv"],
        "Value": [
            round(self.thickness.to("cm").magnitude, 1),
            round(self.length.to("cm").magnitude, 1),
            round(self.height.to("cm").magnitude, 1),
            round(self._hw_lw, 3),
            round(self._Acv.to("cm**2").magnitude, 1),
        ],
        "Unit": ["cm", "cm", "cm", "", "cm²"],
    }

    rho_t_ok = bool(self._rho_t >= self._rho_t_min)
    rho_l_ok = bool(self._rho_l >= self._rho_l_min)
    s_h_ok = bool(self._s_h <= self._s_h_max) if self._s_h > 0 * mm else True
    s_v_ok = bool(self._s_v <= self._s_v_max) if self._s_v > 0 * mm else True
    Vn_max_ok = bool(self._V_u <= self._phi_V_n_max_wall)
    Vn_ok = bool(self._V_u <= self._phi_V_n_wall)

    self._all_wall_shear_checks_passed = all([rho_t_ok, rho_l_ok, Vn_max_ok, Vn_ok])

    def _v(q: Quantity) -> float:
        return (
            round(q.to("kN").magnitude, 2) if self.concrete.unit_system == "metric" else round(q.to("kip").magnitude, 2)
        )

    self._forces_shear_wall = {
        "Design forces": ["Shear"],
        "Variable": ["Vu"],
        "Value": [_v(self._V_u)],
        "Unit": [unit],
    }
    self._shear_capacity_wall = {
        "Shear strength": [
            "Concrete shear strength",
            "Steel shear strength",
            "Total shear strength",
            "Maximum shear strength",
            "Demand Capacity Ratio",
        ],
        "Variable": ["ØVc", "ØVs", "ØVn", "ØVn,max", "DCR"],
        "Value": [
            round((phi_v * self._V_c_wall).to("kN").magnitude, 2),
            round((phi_v * self._V_s_wall).to("kN").magnitude, 2),
            round(self._phi_V_n_wall.to("kN").magnitude, 2),
            round(self._phi_V_n_max_wall.to("kN").magnitude, 2),
            round(self._DCRv_wall, 3),
        ],
        "Unit": [unit, unit, unit, unit, ""],
    }
    self._data_min_max_wall = {
        "Check": [
            "Horizontal reinforcement ratio",
            "Minimum vertical reinf. ratio",
            "Horizontal bar spacing (E.F.)",
            "Vertical bar spacing (E.F.)",
            "Maximum shear capacity",
            "Total shear capacity",
        ],
        "Unit": ["", "", "mm", "mm", unit, unit],
        "Value": [
            round(self._rho_t.to("").magnitude, 5),
            round(self._rho_l.to("").magnitude, 5),
            round(self._s_h.to("mm").magnitude, 1) if self._s_h > 0 * mm else 0.0,
            round(self._s_v.to("mm").magnitude, 1) if self._s_v > 0 * mm else 0.0,
            _v(self._V_u),
            _v(self._V_u),
        ],
        "Min.": [
            round(self._rho_t_min.to("").magnitude, 5),
            round(self._rho_l_min.to("").magnitude, 5),
            "",
            "",
            "",
            "",
        ],
        "Max.": [
            "",
            "",
            round(self._s_h_max.to("mm").magnitude, 1),
            round(self._s_v_max.to("mm").magnitude, 1),
            _v(self._phi_V_n_max_wall),
            _v(self._phi_V_n_wall),
        ],
        "Ok?": [
            "✅" if rho_t_ok else "❌",
            "✅" if rho_l_ok else "❌",
            "✅" if s_h_ok else "❌",
            "✅" if s_v_ok else "❌",
            "✅" if Vn_max_ok else "❌",
            "✅" if Vn_ok else "❌",
        ],
    }


def _compile_results_wall_shear(self: "ShearWall", force: Forces) -> pd.DataFrame:
    """Build a one-row DataFrame for this force combination."""
    phi_v = self.concrete.phi_v

    def _v(q: Quantity) -> float:
        return (
            round(q.to("kN").magnitude, 2) if self.concrete.unit_system == "metric" else round(q.to("kip").magnitude, 2)
        )

    row = {
        "Label": self.label,
        "Comb.": force.label,
        "ρt,min": round(self._rho_t_min.to("").magnitude, 5),
        "ρt,req": round(self._rho_t_req.to("").magnitude, 5),
        "ρt": round(self._rho_t.to("").magnitude, 5),
        "ρl,min": round(self._rho_l_min.to("").magnitude, 5),
        "ρl": round(self._rho_l.to("").magnitude, 5),
        "Vu": _v(self._V_u),
        "ØVc": _v(phi_v * self._V_c_wall),
        "ØVs": _v(phi_v * self._V_s_wall),
        "ØVn": _v(self._phi_V_n_wall),
        "ØVn,max": _v(self._phi_V_n_max_wall),
        "Vu≤ØVn,max": bool(self._V_u <= self._phi_V_n_max_wall),
        "Vu≤ØVn": bool(self._V_u <= self._phi_V_n_wall),
        "DCR": round(self._DCRv_wall, 3),
    }
    return pd.DataFrame([row])


##########################################################
# MAIN CHECK FUNCTION
##########################################################


def _check_shear_ACI_318_19_wall(self: "ShearWall", force: Forces) -> pd.DataFrame:
    """
    ACI 318-19 Section 11 shear check for a structural wall.
    Returns a one-row DataFrame with shear check results.
    """
    if not isinstance(self.concrete, Concrete_ACI_318_19):
        raise TypeError("ACI 318-19 wall shear check requires Concrete_ACI_318_19.")

    concrete = self.concrete

    # 1. Demand
    self._V_u = force._V_z.to("kN")
    self._N_u = force._N_x.to("kN")

    # 2. Geometry: Acv = lw × t
    _calculate_wall_Acv(self)

    # 3. Material: fyt cap
    self._f_yt_wall = _calculate_f_yt_wall(self)

    # 4. α_c based on hw/lw (also sets self._hw_lw)
    self._alpha_c = _calculate_alpha_c(self)

    # 5. Shear strength components
    _calculate_wall_shear_strength(self, concrete)

    # 6. Spacing limits
    _calculate_spacing_limits_wall(self)

    # 7. Required ρt for design
    phi_v = concrete.phi_v
    f_c = concrete.f_c
    lam = concrete.lambda_factor

    if self.concrete.unit_system == "metric":
        Vc_intensity = self._alpha_c * lam * math.sqrt(f_c / MPa) * MPa
    else:
        from mento.units import psi

        Vc_intensity = self._alpha_c * lam * math.sqrt(f_c.to("psi").magnitude) * psi

    rho_t_req_raw = ((self._V_u / phi_v) / self._Acv - Vc_intensity) / self._f_yt_wall
    rho_t_req_raw = rho_t_req_raw.to("")
    self._rho_t_req = max(rho_t_req_raw, self._rho_t_min)

    # 8. Minimum reinforcement ratios (ρl,min depends on ρt,req per §11.6.2)
    _calculate_rho_min_wall(self)

    # 9. DCR
    phi_Vn_eff = min(self._phi_V_n_wall, self._phi_V_n_max_wall)
    if phi_Vn_eff.magnitude == 0:  # pragma: no cover - defensive: Vc > 0 for valid concrete
        self._DCRv_wall = float("inf")
    else:
        self._DCRv_wall = float((self._V_u / phi_Vn_eff).to("").magnitude)

    # 10. Compile detail dicts and result row
    _compile_wall_shear_dicts(self, force)
    return _compile_results_wall_shear(self, force)


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
            f"No standard bar can satisfy ρ ≥ {rho_req:.5f} " f"with spacing in [{s_floor:.0f~P}, {s_max:.0f~P}]."
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
    for force in forces:
        _check_shear_ACI_318_19_wall(self, force)
        max_rho_t_req = max(max_rho_t_req, self._rho_t_req.to("").magnitude)

    # ρl,min per §11.6.2 using the worst-case ρt,req (geometry-only hw/lw).
    r_hw = max(0.5, min(self._hw_lw, 2.5))
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
