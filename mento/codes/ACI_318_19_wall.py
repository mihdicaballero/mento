from __future__ import annotations

import math
from typing import TYPE_CHECKING

import pandas as pd
from pint import Quantity

from mento.material import Concrete_ACI_318_19
from mento.units import MPa, mm, ksi, inch, dimensionless
from mento.forces import Forces

if TYPE_CHECKING:
    from mento.shear_wall import ShearWall


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


def _calculate_wall_shear_strength(self: "ShearWall") -> None:
    """
    ACI 318-19 §11.5.4.6:
        Vc      = α_c × λ × √f'c × Acv
        Vs      = ρt × fyt × Acv
        Vn      = Vc + Vs
        Vn,max  = 0.66 × λ × √f'c × Acv  (metric)
                  8   × λ × √f'c × Acv  (imperial)
    """
    f_c = self.concrete.f_c
    lam = self.concrete.lambda_factor
    phi_v = self.concrete.phi_v

    if self.concrete.unit_system == "metric":
        sqrt_fc_term = math.sqrt(f_c / MPa) * MPa
        Vc = self._alpha_c * lam * sqrt_fc_term * self._Acv
        Vn_max = 0.66 * lam * sqrt_fc_term * self._Acv
    else:
        sqrt_fc_term = math.sqrt(f_c.to("psi").magnitude) * (1 * MPa / MPa * 1 * ksi / ksi)
        # Imperial: coefficients in psi units
        fc_psi = f_c.to("psi").magnitude
        sqrt_fc_psi = math.sqrt(fc_psi)
        from mento.units import psi

        Vc = self._alpha_c * lam * sqrt_fc_psi * psi * self._Acv
        Vn_max = 8.0 * lam * sqrt_fc_psi * psi * self._Acv

    Vs = self._rho_t * self._f_yt_wall * self._Acv
    Vn = Vc + Vs

    self._V_c_wall = Vc.to("kN")
    self._V_s_wall = Vs.to("kN")
    self._V_n_wall = Vn.to("kN")
    self._V_n_max = Vn_max.to("kN")
    self._phi_V_n_wall = (phi_v * min(Vn, Vn_max)).to("kN")
    self._phi_V_n_max_wall = (phi_v * Vn_max).to("kN")


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
            "Vertical reinforcement ratio",
            "Horizontal bar spacing",
            "Vertical bar spacing",
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
    _calculate_wall_shear_strength(self)

    # 6. Spacing limits and min reinforcement
    _calculate_spacing_limits_wall(self)
    _calculate_rho_min_wall(self)

    # 7. Required ρt for design
    phi_v = self.concrete.phi_v
    f_c = self.concrete.f_c
    lam = self.concrete.lambda_factor

    if self.concrete.unit_system == "metric":
        Vc_intensity = self._alpha_c * lam * math.sqrt(f_c / MPa) * MPa
    else:
        from mento.units import psi

        Vc_intensity = self._alpha_c * lam * math.sqrt(f_c.to("psi").magnitude) * psi

    rho_t_req_raw = ((self._V_u / phi_v) / self._Acv - Vc_intensity) / self._f_yt_wall
    rho_t_req_raw = rho_t_req_raw.to("")
    self._rho_t_req = max(rho_t_req_raw, self._rho_t_min)

    # 8. DCR
    phi_Vn_eff = min(self._phi_V_n_wall, self._phi_V_n_max_wall)
    if phi_Vn_eff.magnitude == 0:
        self._DCRv_wall = float("inf")
    else:
        self._DCRv_wall = float((self._V_u / phi_Vn_eff).to("").magnitude)

    # 9. Compile detail dicts and result row
    _compile_wall_shear_dicts(self, force)
    return _compile_results_wall_shear(self, force)


##########################################################
# DESIGN FUNCTION
##########################################################


def _design_shear_ACI_318_19_wall(self: "ShearWall", force: Forces) -> None:
    """
    Compute _rho_t_req without applying rebar selection.
    The check function is called as a side-effect; _rho_t_req is set on self.
    Rebar selection (set_horizontal_rebar) is left to the caller for Phase 0.
    """
    _check_shear_ACI_318_19_wall(self, force)
