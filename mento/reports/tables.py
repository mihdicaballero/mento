"""Report tables for beam checks — the presentation half of a check.

These builders turn the state a check leaves on a beam into the labelled,
rounded, unit-annotated tables the Word reports and notebook views display.
They lived inside the design-code modules, which made those modules roughly a
third UI text; Phase 3 of the architecture roadmap moves them here so
``codes/`` holds equations and orchestration only.

They still read and write the beam's private attributes. That is the coupling
Phase 2b could not remove while presentation was the thing storing results —
confining it to this layer is the step that makes removing it possible.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any, Dict, cast

import pandas as pd
from pint import Quantity

from mento.material import Concrete_ACI_318_19, Concrete_EN_1992_2004
from mento.units import inch, kN, kNm, mm

if TYPE_CHECKING:
    from mento.beam import RectangularBeam
    from mento.forces import Forces
    from mento.settings import BeamSettings


def _settings(beam: RectangularBeam) -> BeamSettings:
    """The beam's settings, which ``__post_init__`` always fills in.

    They are typed ``Optional`` on the section because the constructor accepts
    ``None``; by the time a check has run and there is a report to build, they
    exist.
    """
    return cast("BeamSettings", beam.settings)


def _compile_results_EN_1992_2004_shear(self: "RectangularBeam", force: Forces) -> Dict[str, Any]:
    """One row of the EN shear results table."""
    return {
        "Label": self.label,  # Beam label
        "Comb.": force.label,
        "Av,min": round(self._A_v_min.to("cm ** 2 / m").magnitude, 2),  # Minimum shear reinforcement area
        "Av,req": round(self._A_v_req.to("cm ** 2 / m").magnitude, 2),  # Required shear reinforcing area
        "Av": round(self._A_v.to("cm ** 2 / m").magnitude, 2),  # Provided stirrup reinforcement per unit length
        "NEd": self._N_Ed.to("kN").magnitude,
        "VEd,1": self._V_Ed_1.to("kN").magnitude,  # Max Vu for the design at the support
        "VEd,2": self._V_Ed_2.to("kN").magnitude,  # Max Vu for the design at d from the support
        "VRd,c": round(self._V_Rd_c.to("kN").magnitude, 2),  # Concrete contribution to shear capacity
        "VRd,s": round(self._V_Rd_s.to("kN").magnitude, 2),  # Reinforcement contribution to shear capacity
        "VRd": round(self._V_Rd.to("kN").magnitude, 2),  # Total shear capacity
        "VRd,max": round(self._V_Rd_max.to("kN").magnitude, 2),  # Maximum shear capacity
        "VEd,1≤VRd,max": self._max_shear_ok,  # Check if applied shear is within max shear capacity
        "VEd,2≤VRd": self._V_Ed_2.to("kN").magnitude
        <= self._V_Rd.to("kN").magnitude,  # Check if applied shear is within total capacity
        "DCR": round(self._DCRv, 3),
    }


def build_shear_report(self: "RectangularBeam", force: Forces) -> pd.DataFrame:
    """Result row and detail tables for the shear combination just checked."""
    if self.concrete.design_code in ("ACI 318-19", "CIRSOC 201-25"):
        row = _compile_results_ACI_shear(self, force)
        _initialize_dicts_ACI_318_19_shear(self)
        return row
    row_en = _compile_results_EN_1992_2004_shear(self, force)
    _initialize_dicts_EN_1992_2004_shear(self)
    return pd.DataFrame([row_en], index=[0])


def build_flexure_report(self: "RectangularBeam", force: Forces) -> pd.DataFrame:
    """Result row and detail tables for the flexure combination just checked."""
    if self.concrete.design_code in ("ACI 318-19", "CIRSOC 201-25"):
        results = _compile_results_ACI_flexure_metric(self, force)
        _initialize_dicts_ACI_318_19_flexure(self)
    else:
        results = _compile_results_EN_1992_2004_flexure_metric(self, force)
        _initialize_dicts_EN_1992_2004_flexure(self)
    return pd.DataFrame([results], index=[0])


def _compile_results_ACI_shear(self: "RectangularBeam", force: Forces) -> pd.DataFrame:
    results = {
        "Label": self.label,
        "Comb.": force.label,
        "Av,min": round(self._A_v_min.to("cm²/m").magnitude, 2),
        "Av,req": round(self._A_v_req.to("cm²/m").magnitude, 2),
        "Av": round(self._A_v.to("cm²/m").magnitude, 2),
        "Vu": self._V_u.to("kN").magnitude,
        "Nu": self._N_u.to("kN").magnitude,
        "ØVc": round(self._phi_V_c.to("kN").magnitude, 2),
        "ØVs": round(self._phi_V_s.to("kN").magnitude, 2),
        "ØVn": round(self._phi_V_n.to("kN").magnitude, 2),
        "ØVmax": round(self._phi_V_max.to("kN").magnitude, 2),
        "Vu≤ØVmax": self._max_shear_ok,
        "Vu≤ØVn": self._V_u <= self._phi_V_n,
        "DCR": round(self._DCRv, 3),
    }
    return pd.DataFrame([results], index=[0])


def _compile_results_ACI_flexure_metric(self: "RectangularBeam", force: Forces) -> Dict[str, Any]:
    # Create dictionaries for bottom and top rows
    if self._M_u >= 0:
        result = {
            "Label": self.label,
            "Comb.": force.label,
            "Position": "Bottom",
            "As,min": round(self._A_s_min_bot.to("cm ** 2").magnitude, 2),
            "As,req top": round(self._A_s_req_top.to("cm ** 2").magnitude, 2),
            "As,req bot": round(self._A_s_req_bot.to("cm ** 2").magnitude, 2),
            "As": round(self._A_s_bot.to("cm ** 2").magnitude, 2),
            # 'c/d': self._c_d_bot,
            "Mu": round(self._M_u_bot.to("kN*m").magnitude, 2),
            "ØMn": round(self._phi_M_n_bot.to("kN*m").magnitude, 2),
            "Mu≤ØMn": self._M_u_bot <= self._phi_M_n_bot,
            "DCR": round(self._DCRb_bot, 3),
        }
    else:
        result = {
            "Label": self.label,
            "Comb.": force.label,
            "Position": "Top",
            "As,min": round(self._A_s_min_top.to("cm ** 2").magnitude, 2),
            "As,req top": round(self._A_s_req_top.to("cm ** 2").magnitude, 2),
            "As,req bot": round(self._A_s_req_bot.to("cm ** 2").magnitude, 2),
            "As": round(self._A_s_top.to("cm ** 2").magnitude, 2),
            # 'c/d': self._c_d_top,
            "Mu": round(self._M_u_top.to("kN*m").magnitude, 2),
            "ØMn": round(self._phi_M_n_top.to("kN*m").magnitude, 2),
            "Mu≤ØMn": -self._M_u_top <= self._phi_M_n_top,
            "DCR": round(self._DCRb_top, 3),
        }
    return result


def _initialize_dicts_ACI_318_19_shear(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        """Initialize the dictionaries used in check and design methods."""
        self._materials_shear = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
                "Concrete density",
                "Normalweight concrete",
                "Safety factor for shear",
            ],
            "Variable": ["", "fc", "fy", "wc", "λ", "Øv"],
            "Value": [
                self.label,
                round(self.concrete.f_c.to("MPa").magnitude, 2),
                round(self.steel_bar.f_y.to("MPa").magnitude, 2),
                round(self.concrete.density.to("kg/m**3").magnitude, 1),
                self.concrete.lambda_factor,
                self.concrete.phi_v,
            ],
            "Unit": ["", "MPa", "MPa", "kg/m³", "", ""],
        }
        self._geometry_shear = {
            "Geometry": [
                "Section height",
                "Section width",
                "Clear cover",
                "Longitudinal tension rebar",
            ],
            "Variable": ["h", "b", "cc", "As"],
            "Value": [
                self.height.to("cm").magnitude,
                self.width.to("cm").magnitude,
                self.c_c.to("cm").magnitude,
                round(self._A_s_tension.to("cm**2").magnitude, 2),
            ],
            "Unit": ["cm", "cm", "cm", "cm²"],
        }
        self._forces_shear = {
            "Design forces": [
                "Axial, positive for compression",
                "Shear",
            ],
            "Variable": ["Nu", "Vu"],
            "Value": [
                round(self._N_u.to("kN").magnitude, 2),
                round(self._V_u.to("kN").magnitude, 2),
            ],
            "Unit": ["kN", "kN"],
        }
        # Min max lists
        # With no stirrup contribution there is no stirrup to report. This used to
        # be done by zeroing the section's own diameter, which made building a
        # report change the section it describes; it is a local now.
        zero_d_b = 0 * mm if self.concrete.unit_system == "metric" else 0 * inch
        d_b_shown = self._stirrup_d_b
        if self._phi_V_s == 0 * kN:
            db_min = zero_d_b
            d_b_shown = zero_d_b
        else:
            if self.concrete.design_code == "ACI 318-19":
                db_min = 10 * mm if self.concrete.unit_system == "metric" else 3 / 8 * inch
            else:
                db_min = 6 * mm
        min_values = [
            None,
            None,
            self._A_v_min,
            db_min,
        ]  # Use None for items without a minimum constraint
        max_values = [
            self._stirrup_s_max_l,
            self._stirrup_s_max_w,
            None,
            None,
        ]  # Use None for items without a maximum constraint
        current_values = [
            self._stirrup_s_l,
            self._stirrup_s_w,
            self._A_v,
            d_b_shown,
        ]  # Current values to check
        # Generate check marks based on the range conditions
        checks = [
            "✅" if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else "❌"
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_shear_checks_passed = all(check == "✅" for check in checks)
        self._data_min_max_shear = {
            "Check": [
                "Stirrup spacing along length",
                "Stirrup spacing along width",
                "Minimum shear reinforcement",
                "Minimum rebar diameter",
            ],
            "Unit": ["cm", "cm", "cm²/m", "mm"],
            "Value": [
                round(self._stirrup_s_l.to("cm").magnitude, 2),
                round(self._stirrup_s_w.to("cm").magnitude, 2),
                round(self._A_v.to("cm**2/m").magnitude, 2),
                round(d_b_shown.to("mm").magnitude, 0),
            ],
            "Min.": [
                "",
                "",
                round(self._A_v_min.to("cm**2/m").magnitude, 2),
                round(db_min.to("mm").magnitude, 0),
            ],
            "Max.": [
                round(self._stirrup_s_max_l.to("cm").magnitude, 2),
                round(self._stirrup_s_max_w.to("cm").magnitude, 2),
                "",
                "",
            ],
            "Ok?": checks,
        }
        self._shear_reinforcement = {
            "Shear reinforcement strength": [
                "Number of stirrups",
                "Stirrup diameter",
                "Stirrup spacing",
                "Effective height",
                "Minimum shear reinforcing",
                "Required shear reinforcing",
                "Defined shear reinforcing",
                "Shear rebar strength",
            ],
            "Variable": ["ns", "db", "s", "d", "Av,min", "Av,req", "Av", "ØVs"],
            "Value": [
                self._stirrup_n,
                round(d_b_shown.to("mm").magnitude, 3),
                round(self._stirrup_s_l.to("cm").magnitude, 3),
                round(self._d_shear.to("cm").magnitude, 2),
                round(self._A_v_min.to("cm**2/m").magnitude, 2),
                round(self._A_v_req.to("cm**2/m").magnitude, 2),
                round(self._A_v.to("cm**2/m").magnitude, 2),
                round(self._phi_V_s.to("kN").magnitude, 2),
            ],
            "Unit": ["", "mm", "cm", "cm", "cm²/m", "cm²/m", "cm²/m", "kN"],
        }
        check_max = "✅" if self._max_shear_ok else "❌"
        check_FU = "✅" if self._DCRv < 1 else "❌"
        self._shear_concrete = {
            "Shear strength": [
                "Effective shear area",
                "Longitudinal reinforcement ratio",
                "Size modification factor",
                "Axial stress",
                "Concrete effective shear stress",
                "Concrete strength",
                "Maximum shear strength",
                "Total shear strength",
                "Max shear check",
                "Demand Capacity Ratio",
            ],
            "Variable": [
                "Acv",
                "ρw",
                "λs",
                "σNu",
                "kc",
                "ØVc",
                "ØVmax",
                "ØVn",
                "",
                "DCR",
            ],
            "Value": [
                round(self._A_cv.to("cm**2").magnitude, 2),
                round(self._rho_w.magnitude, 5),
                round(self._lambda_s, 3),
                round(self._sigma_Nu.to("MPa").magnitude, 2),
                round(self._k_c_min.to("MPa").magnitude, 2),
                round(self._phi_V_c.to("kN").magnitude, 2),
                round(self._phi_V_max.to("kN").magnitude, 2),
                round(self._phi_V_n.to("kN").magnitude, 2),
                check_max,
                round(self._DCRv, 2),
            ],
            "Unit": ["cm²", "", "", "MPa", "MPa", "kN", "kN", "kN", "", check_FU],
        }
        self._shear_all_checks = self._all_shear_checks_passed and (check_max == "✅") and (check_FU == "✅")


def _initialize_dicts_ACI_318_19_flexure(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        # Update longitudinal rebar attributes
        self._update_longitudinal_rebar_attributes()
        """Initialize the dictionaries used in check and design methods."""
        self._materials_flexure = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
            ],
            "Variable": ["", "fc", "fy"],
            "Value": [
                self.label,
                round(self.concrete.f_c.to("MPa").magnitude, 2),
                round(self.steel_bar.f_y.to("MPa").magnitude, 2),
            ],
            "Unit": ["", "MPa", "MPa"],
        }
        self._geometry_flexure = {
            "Geometry": [
                "Section height",
                "Section width",
                "Clear cover",
                "Mechanical top cover",
                "Mechanical bottom cover",
            ],
            "Variable": ["h", "b", "cc", "cm,top", "cm,bot"],
            "Value": [
                self.height.to("cm").magnitude,
                self.width.to("cm").magnitude,
                self.c_c.to("cm").magnitude,
                round(self._c_mec_top.to("cm").magnitude, 2),
                round(self._c_mec_bot.to("cm").magnitude, 2),
            ],
            "Unit": ["cm", "cm", "cm", "cm", "cm"],
        }
        self._forces_flexure = {
            "Design forces": [
                "Top max moment",
                "Bottom max moment",
            ],
            "Variable": ["Mu,top", "Mu,bot"],
            "Value": [
                round(self._M_u_top.to("kN*m").magnitude, 2),
                round(self._M_u_bot.to("kN*m").magnitude, 2),
            ],
            "Unit": ["kNm", "kNm"],
        }
        # Min max lists
        settings = _settings(self)
        min_spacing_top: Quantity = max(settings.clear_spacing, settings.vibrator_size, self._d_b_max_top)
        min_spacing_bot: Quantity = max(settings.clear_spacing, self._d_b_max_bot)
        min_values = [
            self._A_s_min_top,
            min_spacing_top,
            self._A_s_min_bot,
            min_spacing_bot,
        ]  # Use None for items without a minimum constraint
        max_values = [
            self._A_s_max_top,
            None,
            self._A_s_max_bot,
            None,
        ]  # Use None for items without a maximum constraint
        current_values = [
            self._A_s_top,
            self._available_s_top,
            self._A_s_bot,
            self._available_s_bot,
        ]  # Current values to check

        ARTICLE_STR = "9.6.1.3"

        checks = []
        for i, (curr, min_val, max_val) in enumerate(zip(current_values, min_values, max_values)):
            # --- EXCEPTION FOR DOUBLY REINFORCED SECTIONS ---
            # If doubly reinforced, ignore maximum limits for top (i=0) and bottom (i=2)
            if self._doubly_reinforced and i in (0, 2):
                # If it passes min, we give the special tag
                if min_val is None or curr >= min_val:
                    checks.append("✅ D.R.")
                    continue
                # If it fails min, let the normal logic handle it (fall through)
            # -------------------------------------------------

            passed = (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val)
            if passed:
                checks.append("✅")
                continue

            # Detect if fails by min or max
            failed_min = (min_val is not None) and (curr < min_val)

            if i == 0 and failed_min and getattr(self, "_A_s_bool_top", True):
                # Position 0 -> _A_s_top vs _A_s_min_top
                checks.append(ARTICLE_STR)
            elif i == 2 and failed_min and getattr(self, "_A_s_bool_bot", True):
                # Position 2 -> _A_s_bot vs _A_s_min_bot
                checks.append(ARTICLE_STR)
            else:
                # Any other failure (includes max or no flags)
                checks.append("❌")

        self._all_flexure_checks_passed = not any(check in ("❌") for check in checks)
        self._data_min_max_flexure = {
            "Check": [
                "Min/Max As rebar top",
                "Minimum spacing top",
                "Min/Max As rebar bottom",
                "Minimum spacing bottom",
            ],
            "Unit": ["cm²", "mm", "cm²", "mm"],
            "Value": [
                round(self._A_s_top.to("cm**2").magnitude, 2),
                round(self._available_s_top.to("mm").magnitude, 2),
                round(self._A_s_bot.to("cm**2").magnitude, 2),
                round(self._available_s_bot.to("mm").magnitude, 2),
            ],
            "Min.": [
                round(self._A_s_min_top.to("cm**2").magnitude, 2),
                round(min_spacing_top.to("mm").magnitude, 2),
                round(self._A_s_min_bot.to("cm**2").magnitude, 2),
                round(min_spacing_bot.to("mm").magnitude, 2),
            ],
            "Max.": [
                round(self._A_s_max_top.to("cm**2").magnitude, 2),
                "",
                round(self._A_s_max_bot.to("cm**2").magnitude, 2),
                "",
            ],
            "Ok?": checks,
        }
        check_DCR_top = "✅" if self._DCRb_top < 1 else "❌"
        check_DCR_bot = "✅" if self._DCRb_bot < 1 else "❌"
        self._flexure_capacity_top = {
            "Top reinforcement check": [
                "First layer bars",
                "Second layer bars",
                "Effective height",
                "Depth of equivalent strength block ratio",
                "Minimum rebar reinforcing",
                "Required rebar reinforcing top",
                "Required rebar reinforcing bottom",
                "Defined rebar reinforcing top",
                "Longitudinal reinforcement ratio",
                "Total flexural strength",
                "Demand Capacity Ratio",
            ],
            "Variable": [
                "n1+n2",
                "n3+n4",
                "d",
                "c/d",
                "As,min",
                "As,req top",
                "As,req bot",
                "As",
                "ρl",
                "ØMn",
                "DCR",
            ],
            "Value": [
                self._format_longitudinal_rebar_string(self._n1_t, self._d_b1_t, self._n2_t, self._d_b2_t),
                self._format_longitudinal_rebar_string(self._n3_t, self._d_b3_t, self._n4_t, self._d_b4_t),
                round(self._d_top.to("cm").magnitude, 2),
                self._c_d_top,
                round(self._A_s_min_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_bot.to("cm**2").magnitude, 2),
                round(self._A_s_top.to("cm**2").magnitude, 2),
                round(self._rho_l_top.magnitude, 5),
                round(self._phi_M_n_top.to("kN*m").magnitude, 2),
                round(self._DCRb_top, 2),
            ],
            "Unit": [
                "",
                "",
                "cm",
                "",
                "cm²",
                "cm²",
                "cm²",
                "cm²",
                "",
                "kNm",
                check_DCR_top,
            ],
        }
        self._flexure_capacity_bot = {
            "Bottom reinforcement check": [
                "First layer bars",
                "Second layer bars",
                "Effective height",
                "Depth of equivalent strength block ratio",
                "Minimum rebar reinforcing",
                "Required rebar reinforcing top",
                "Required rebar reinforcing bottom",
                "Defined rebar reinforcing bottom",
                "Longitudinal reinforcement ratio",
                "Total flexural strength",
                "Demand Capacity Ratio",
            ],
            "Variable": [
                "n1+n2",
                "n3+n4",
                "d",
                "c/d",
                "As,min",
                "As,req top",
                "As,req bot",
                "As",
                "ρl",
                "ØMn",
                "DCR",
            ],
            "Value": [
                self._format_longitudinal_rebar_string(self._n1_b, self._d_b1_b, self._n2_b, self._d_b2_b),
                self._format_longitudinal_rebar_string(self._n3_b, self._d_b3_b, self._n4_b, self._d_b4_b),
                round(self._d_bot.to("cm").magnitude, 2),
                self._c_d_bot,
                round(self._A_s_min_bot.to("cm**2").magnitude, 2),
                round(self._A_s_req_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_bot.to("cm**2").magnitude, 2),
                round(self._A_s_bot.to("cm**2").magnitude, 2),
                round(self._rho_l_bot.magnitude, 5),
                round(self._phi_M_n_bot.to("kN*m").magnitude, 2),
                round(self._DCRb_bot, 2),
            ],
            "Unit": [
                "",
                "",
                "cm",
                "",
                "cm²",
                "cm²",
                "cm²",
                "cm²",
                "",
                "kNm",
                check_DCR_bot,
            ],
        }
        self._flexure_all_checks = (
            self._all_flexure_checks_passed and (check_DCR_bot == "✅") and (check_DCR_top == "✅")
        )


def _compile_results_EN_1992_2004_flexure_metric(self: "RectangularBeam", force: Forces) -> Dict[str, Any]:
    # Create dictionaries for bottom and top rows
    if self._M_Ed >= 0:
        result = {
            "Label": self.label,
            "Comb.": force.label,
            "Position": "Bottom",
            "As,min": round(self._A_s_min_bot.to("cm ** 2").magnitude, 2),
            "As,req top": round(self._A_s_req_top.to("cm ** 2").magnitude, 2),
            "As,req bot": round(self._A_s_req_bot.to("cm ** 2").magnitude, 2),
            "As": round(self._A_s_bot.to("cm ** 2").magnitude, 2),
            # 'c/d': self._c_d_bot,
            "MEd": round(self._M_Ed_bot.to("kN*m").magnitude, 2),
            "MRd": round(self._M_Rd_bot.to("kN*m").magnitude, 2),
            "MEd≤MRd": self._M_Ed_bot <= self._M_Rd_bot,
            "DCR": round(self._DCRb_bot, 3),
        }
    else:
        result = {
            "Label": self.label,
            "Comb.": force.label,
            "Position": "Top",
            "As,min": round(self._A_s_min_top.to("cm ** 2").magnitude, 2),
            "As,req top": round(self._A_s_req_top.to("cm ** 2").magnitude, 2),
            "As,req bot": round(self._A_s_req_bot.to("cm ** 2").magnitude, 2),
            "As": round(self._A_s_top.to("cm ** 2").magnitude, 2),
            # 'c/d': self._c_d_top,
            "MEd": round(self._M_Ed_top.to("kN*m").magnitude, 2),
            "MRd": round(self._M_Rd_top.to("kN*m").magnitude, 2),
            "MEd≤MRd": self._M_Ed_top <= self._M_Rd_top,
            "DCR": round(self._DCRb_top, 3),
        }
    return result


def _initialize_dicts_EN_1992_2004_shear(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        """Initialize the dictionaries used in check and design methods."""
        self._materials_shear = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
                "Safety factor for concrete",
                "Safety factor for steel",
                "Coefficient for long term effects and loading effects",
            ],
            "Variable": ["", "fck", "fywk", "γc", "γs", "αcc"],
            "Value": [
                self.label,
                round(self.concrete.f_ck.to("MPa").magnitude, 2),
                round(self.steel_bar.f_y.to("MPa").magnitude, 2),
                self.concrete.gamma_c,
                self.concrete._gamma_s,
                self.concrete.alpha_cc,
            ],
            "Unit": ["", "MPa", "MPa", "", "", ""],
        }
        self._geometry_shear = {
            "Geometry": [
                "Section height",
                "Section width",
                "Clear cover",
                "Longitudinal tension rebar",
            ],
            "Variable": ["h", "b", "cc", "As"],
            "Value": [
                self.height.to("cm").magnitude,
                self.width.to("cm").magnitude,
                self.c_c.to("cm").magnitude,
                round(self._A_s_tension.to("cm**2").magnitude, 2),
            ],
            "Unit": ["cm", "cm", "cm", "cm²"],
        }
        self._forces_shear = {
            "Design forces": [
                "Axial, positive for compression",
                "Shear",
            ],
            "Variable": ["NEd", "VEd,2"],
            "Value": [
                round(self._N_Ed.to("kN").magnitude, 2),
                round(self._V_Ed_2.to("kN").magnitude, 2),
            ],
            "Unit": ["kN", "kN"],
        }
        # Min max lists
        # Same as the ACI table above: a local, not a write to the section.
        d_b_shown = self._stirrup_d_b
        if self._V_Rd_s == 0 * kN:
            d_b_shown = 0 * mm if self.concrete.unit_system == "metric" else 0 * inch
        # Min max lists
        min_values = [
            None,
            None,
            self._A_v_min,
        ]  # Use None for items without a minimum constraint
        max_values = [
            self._stirrup_s_max_l,
            self._stirrup_s_max_w,
            None,
        ]  # Use None for items without a maximum constraint
        current_values = [
            self._stirrup_s_l,
            self._stirrup_s_w,
            self._A_v,
        ]  # Current values to check

        # Generate check marks based on the range conditions
        checks = [
            "✅" if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else "❌"
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_shear_checks_passed = all(check == "✅" for check in checks)
        self._data_min_max_shear = {
            "Check": [
                "Stirrup spacing along length",
                "Stirrup spacing along width",
                "Minimum shear reinforcement",
            ],
            "Unit": ["cm", "cm", "cm²/m"],
            "Value": [
                round(self._stirrup_s_l.to("cm").magnitude, 2),
                round(self._stirrup_s_w.to("cm").magnitude, 2),
                round(self._A_v.to("cm**2/m").magnitude, 2),
            ],
            "Min.": ["", "", round(self._A_v_min.to("cm**2/m").magnitude, 2)],
            "Max.": [
                round(self._stirrup_s_max_l.to("cm").magnitude, 2),
                round(self._stirrup_s_max_w.to("cm").magnitude, 2),
                "",
            ],
            "Ok?": checks,
        }
        self._shear_reinforcement = {
            "Shear reinforcement strength": [
                "Number of stirrups",
                "Stirrup diameter",
                "Stirrup spacing",
                "Effective height",
                "Minimum shear reinforcing",
                "Required shear reinforcing",
                "Defined shear reinforcing",
                "Shear rebar strength",
            ],
            "Variable": ["ns", "db", "s", "d", "Asw,min", "Asw,req", "Asw", "VRd,s"],
            "Value": [
                self._stirrup_n,
                d_b_shown.to("mm").magnitude,
                self._stirrup_s_l.to("cm").magnitude,
                round(self._d_shear.to("cm").magnitude, 2),
                round(self._A_v_min.to("cm**2/m").magnitude, 2),
                round(self._A_v_req.to("cm**2/m").magnitude, 2),
                round(self._A_v.to("cm**2/m").magnitude, 2),
                round(self._V_Rd_s.to("kN").magnitude, 2),
            ],
            "Unit": ["", "mm", "cm", "cm", "cm²/m", "cm²/m", "cm²/m", "kN"],
        }
        check_max = "✅" if self._max_shear_ok else "❌"
        check_DCR = "✅" if self._DCRv < 1 else "❌"
        rho_l = self._rho_l_bot if self._M_Ed >= 0 * kNm else self._rho_l_top
        self._shear_concrete = {
            "Shear strength": [
                "Longitudinal reinforcement ratio",
                "k value",
                "Axial stress",
                "Concrete strut angle",
                "Concrete strength",
                "Maximum shear strength",
                "Total shear strength",
                "Max shear check",
                "Demand Capacity Ratio",
            ],
            "Variable": ["ρl", "k", "σcd", "Θ", "VRd,c", "VRd,max", "VRd", "", "DCR"],
            "Value": [
                round(rho_l.magnitude, 4),
                round(self._k_value, 2),
                round(self._sigma_cp.to("MPa").magnitude, 2),
                round(math.degrees(self._theta), 1),
                round(self._V_Rd_c.to("kN").magnitude, 2),
                round(self._V_Rd_max.to("kN").magnitude, 2),
                round(self._V_Rd.to("kN").magnitude, 2),
                check_max,
                round(self._DCRv, 3),
            ],
            "Unit": ["", "", "MPa", "deg", "kN", "kN", "kN", "", check_DCR],
        }


def _initialize_dicts_EN_1992_2004_flexure(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Update longitudinal rebar attributes
        self._update_longitudinal_rebar_attributes()
        """Initialize the dictionaries used in check and design methods."""
        self._materials_flexure = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
            ],
            "Variable": ["", "fck", "fyk"],
            "Value": [
                self.label,
                round(self.concrete.f_ck.to("MPa").magnitude, 2),
                round(self.steel_bar.f_y.to("MPa").magnitude, 2),
            ],
            "Unit": ["", "MPa", "MPa"],
        }
        self._geometry_flexure = {
            "Geometry": [
                "Section height",
                "Section width",
                "Clear cover",
                "Mechanical top cover",
                "Mechanical bottom cover",
            ],
            "Variable": ["h", "b", "cc", "cm,top", "cm,bot"],
            "Value": [
                self.height.to("cm").magnitude,
                self.width.to("cm").magnitude,
                self.c_c.to("cm").magnitude,
                round(self._c_mec_top.to("cm").magnitude, 2),
                round(self._c_mec_bot.to("cm").magnitude, 2),
            ],
            "Unit": ["cm", "cm", "cm", "cm", "cm"],
        }
        self._forces_flexure = {
            "Design forces": [
                "Top max moment",
                "Bottom max moment",
            ],
            "Variable": ["MEd,top", "MEd,bot"],
            "Value": [
                round(self._M_Ed_top.to("kN*m").magnitude, 2),
                round(self._M_Ed_bot.to("kN*m").magnitude, 2),
            ],
            "Unit": ["kNm", "kNm"],
        }
        # Min max lists
        settings = _settings(self)
        min_spacing_top: Quantity = max(settings.clear_spacing, settings.vibrator_size, self._d_b_max_top)
        min_spacing_bot: Quantity = max(settings.clear_spacing, self._d_b_max_bot)
        min_values = [
            self._A_s_min_top,
            min_spacing_top,
            self._A_s_min_bot,
            min_spacing_bot,
        ]  # Use None for items without a minimum constraint
        max_values = [
            self._A_s_max_top,
            None,
            self._A_s_max_bot,
            None,
        ]  # Use None for items without a maximum constraint
        current_values = [
            self._A_s_top,
            self._available_s_top,
            self._A_s_bot,
            self._available_s_bot,
        ]  # Current values to check

        # Generate check marks based on the range conditions
        checks = []
        for i, (curr, min_val, max_val) in enumerate(zip(current_values, min_values, max_values)):
            # --- EXCEPTION FOR DOUBLY REINFORCED SECTIONS ---
            # If doubly reinforced, ignore maximum limits for top (i=0) and bottom (i=2)
            if self._doubly_reinforced and i in (0, 2):
                # If it passes min, we give the special tag
                if min_val is None or curr >= min_val:
                    checks.append("✅ D.R.")
                    continue
                # If it fails min, let the normal logic handle it (fall through)
            # -------------------------------------------------

            passed = (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val)
            if passed:
                checks.append("✅")
            else:
                checks.append("❌")
        self._all_flexure_checks_passed = not any(check in ("❌") for check in checks)
        self._data_min_max_flexure = {
            "Check": [
                "Min/Max As rebar top",
                "Minimum spacing top",
                "Min/Max As rebar bottom",
                "Minimum spacing bottom",
            ],
            "Unit": ["cm²", "mm", "cm²", "mm"],
            "Value": [
                round(self._A_s_top.to("cm**2").magnitude, 2),
                round(self._available_s_top.to("mm").magnitude, 2),
                round(self._A_s_bot.to("cm**2").magnitude, 2),
                round(self._available_s_bot.to("mm").magnitude, 2),
            ],
            "Min.": [
                round(self._A_s_min_top.to("cm**2").magnitude, 2),
                round(min_spacing_top.to("mm").magnitude, 2),
                round(self._A_s_min_bot.to("cm**2").magnitude, 2),
                round(min_spacing_bot.to("mm").magnitude, 2),
            ],
            "Max.": [
                round(self._A_s_max_top.to("cm**2").magnitude, 2),
                "",
                round(self._A_s_max_bot.to("cm**2").magnitude, 2),
                "",
            ],
            "Ok?": checks,
        }
        check_DCR_top = "✅" if self._DCRb_top < 1 else "❌"
        check_DCR_bot = "✅" if self._DCRb_bot < 1 else "❌"
        self._flexure_capacity_top = {
            "Top reinforcement check": [
                "First layer bars",
                "Second layer bars",
                "Effective height",
                "Depth of equivalent strength block ratio",
                "Minimum rebar reinforcing",
                "Required rebar reinforcing top",
                "Required rebar reinforcing bottom",
                "Defined rebar reinforcing top",
                "Longitudinal reinforcement ratio",
                "Total flexural strength",
                "Demand Capacity Ratio",
            ],
            "Variable": [
                "n1+n2",
                "n3+n4",
                "d",
                "c/d",
                "As,min",
                "As,req top",
                "As,req bot",
                "As",
                "ρl",
                "MRd",
                "DCR",
            ],
            "Value": [
                self._format_longitudinal_rebar_string(self._n1_t, self._d_b1_t, self._n2_t, self._d_b2_t),
                self._format_longitudinal_rebar_string(self._n3_t, self._d_b3_t, self._n4_t, self._d_b4_t),
                round(self._d_top.to("cm").magnitude, 2),
                self._c_d_top,
                round(self._A_s_min_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_bot.to("cm**2").magnitude, 2),
                round(self._A_s_top.to("cm**2").magnitude, 2),
                round(self._rho_l_top.magnitude, 5),
                round(self._M_Rd_top.to("kN*m").magnitude, 2),
                round(self._DCRb_top, 2),
            ],
            "Unit": [
                "",
                "",
                "cm",
                "",
                "cm²",
                "cm²",
                "cm²",
                "cm²",
                "",
                "kNm",
                check_DCR_top,
            ],
        }
        self._flexure_capacity_bot = {
            "Bottom reinforcement check": [
                "First layer bars",
                "Second layer bars",
                "Effective height",
                "Depth of equivalent strength block ratio",
                "Minimum rebar reinforcing",
                "Required rebar reinforcing top",
                "Required rebar reinforcing bottom",
                "Defined rebar reinforcing bottom",
                "Longitudinal reinforcement ratio",
                "Total flexural strength",
                "Demand Capacity Ratio",
            ],
            "Variable": [
                "n1+n2",
                "n3+n4",
                "d",
                "c/d",
                "As,min",
                "As,req top",
                "As,req bot",
                "As",
                "ρl",
                "MRd",
                "DCR",
            ],
            "Value": [
                self._format_longitudinal_rebar_string(self._n1_b, self._d_b1_b, self._n2_b, self._d_b2_b),
                self._format_longitudinal_rebar_string(self._n3_b, self._d_b3_b, self._n4_b, self._d_b4_b),
                round(self._d_bot.to("cm").magnitude, 2),
                self._c_d_bot,
                round(self._A_s_min_bot.to("cm**2").magnitude, 2),
                round(self._A_s_req_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_bot.to("cm**2").magnitude, 2),
                round(self._A_s_bot.to("cm**2").magnitude, 2),
                round(self._rho_l_bot.magnitude, 5),
                round(self._M_Rd_bot.to("kN*m").magnitude, 2),
                round(self._DCRb_bot, 2),
            ],
            "Unit": [
                "",
                "",
                "cm",
                "",
                "cm²",
                "cm²",
                "cm²",
                "cm²",
                "",
                "kNm",
                check_DCR_bot,
            ],
        }
        self._flexure_all_checks = (
            self._all_flexure_checks_passed and (check_DCR_bot == "✅") and (check_DCR_top == "✅")
        )
