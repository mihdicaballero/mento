"""EN 1992-2004. The registry entry."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from mento.codes.check_state import apply_en_flexure_state, apply_en_shear_state
from mento.codes.EN_1992_2004_beam import (
    _check_flexure_EN_1992_2004,
    _check_shear_EN_1992_2004,
    _design_flexure_EN_1992_2004,
    _design_shear_EN_1992_2004,
)
from mento.codes.registry import DesignCode, register
from mento.material import Concrete_EN_1992_2004
from mento.units import cm, kN, kNm, mm, MPa

if TYPE_CHECKING:
    from mento.beam import RectangularBeam
    from mento.rebar import Rebar

#: EN names the demand M_Ed and the capacity M_Rd, where ACI says M_u / phi*Mn.
_FLEXURE_SYMBOLS = {
    "demand": "MEd",
    "demand_top": "MEd,top",
    "demand_bot": "MEd,bot",
    "md_demand": "M_{Ed}",
    "md_capacity": "M_{Rd}",
}

_UNITS_ROW_SHEAR = {
    "Label": "",
    "Comb.": "",
    "Av,min": "cm²/m",
    "Av,req": "cm²/m",
    "Av": "cm²/m",
    "NEd": "kN",
    "VEd,1": "kN",
    "VEd,2": "kN",
    "VRd,c": "kN",
    "VRd,s": "kN",
    "VRd": "kN",
    "VRd,max": "kN",
    "VEd,1≤VRd,max": "",
    "VEd,2≤VRd": "",
    "DCR": "",
}

_UNITS_ROW_FLEXURE = {
    "Label": "",
    "Comb.": "",
    "Position": "",
    "As,min": "cm²",
    "As,req top": "cm²",
    "As,req bot": "cm²",
    "As": "cm²",
    "MEd": "kNm",
    "MRd": "kNm",
    "MEd≤MRd": "",
    "DCR": "",
}


_SHEAR_SYMBOLS = {
    "reinforcement": "A_{sw}",
    "demand": "V_{Ed,2}",
    "capacity": "V_{Rd}",
    # Row the design resistance sits on in this code's shear detail table.
    "capacity_row": 6,
}

#: Columns the all-beams summaries leave out. The ticks repeat what the DCR
#: column beside them already says, and the resistances are right there for a
#: reader who wants the margin; VEd,1 is the demand at the support edge, which
#: the detailed output carries. They all stay in `check_shear` and
#: `check_flexure`, which is where a caller reads them programmatically.
_SUMMARY_DROP_COLUMNS = ("VEd,1", "VEd,1≤VRd,max", "VEd,2≤VRd", "MEd≤MRd")


_SUMMARY_COLUMNS = {
    "moment_demand": "MEd",
    "shear_demand": "VEd",
    # EN's shear table reports the demand at the edge and at d; the summary
    # takes the one at d.
    "shear_demand_source": "VEd,2",
    "axial_demand": "NEd",
    "moment_capacity_top": "MRd,top",
    "moment_capacity_bot": "MRd,bot",
    "shear_capacity": "VRd",
}


def _capacity_columns(section: "RectangularBeam") -> dict:
    """EN reports a design resistance M_Rd on each face."""
    return {
        "MRd,top": round(section._M_Rd_top.to("kN*m").magnitude, 1),
        "MRd,bot": round(section._M_Rd_bot.to("kN*m").magnitude, 1),
    }


def _initialize_attributes(section: "RectangularBeam") -> None:
    """The zeroed result attributes the EN report tables read off the beam.

    Moved here from ``RectangularBeam`` verbatim; see the ACI counterpart.
    """
    section._V_Ed_1 = 0 * kN
    section._V_Ed_2 = 0 * kN
    section._N_Ed = 0 * kN
    section._M_Ed = 0 * kNm
    section._sigma_cd = 0 * MPa
    section._V_Rd_c = 0 * kN
    section._V_Rd_s = 0 * kN
    section._V_Rd_max = 0 * kN
    section._V_Rd = 0 * kN
    section._k_value = 0
    section._f_ywk = section.steel_bar.f_y
    section._f_ywd = 0 * MPa
    section._f_cd = 0 * MPa
    section._f_cd_shear = 0 * MPa
    section._A_p = 0 * cm**2  # No prestressed for now
    section._sigma_cp = 0 * MPa
    section._theta = 0
    section._cot_theta = 0
    section._z = 0 * cm
    section._M_Rd_bot = 0 * kNm
    section._M_Rd_top = 0 * kNm
    section._M_Ed_bot = 0 * kNm
    section._M_Ed_top = 0 * kNm
    # Shared with the ACI branch: the flexural design driver and the
    # results tables read these before any check has run.
    section._A_s_min_bot = 0 * cm**2
    section._A_s_min_top = 0 * cm**2
    section._A_s_max_bot = 0 * cm**2
    section._A_s_max_top = 0 * cm**2
    section.flexure_design_results_bot = None
    section.flexure_design_results_top = None


def _transverse_rebar(rebar: "Rebar", V_s_req: Any, alpha: float) -> Any:
    return rebar.transverse_rebar_EN_1992_2004(alpha)


def _longitudinal_rebar(rebar: "Rebar", A_s_req: Any, A_s_max: Any, mech_cover: Any) -> Any:
    return rebar.longitudinal_rebar_EN_1992_2004(A_s_req, A_s_max, mech_cover)


def _max_bar_spacing_slab(section: "RectangularBeam") -> Any:
    """EN 1992-1-1 9.3.1.1(3): 3h, and no more than 400 mm.

    The general provision for the principal reinforcement of a slab. The tighter
    one it states for areas of concentrated load or of maximum moment (2h, and
    no more than 250 mm) is not applied here: a section is checked against a set
    of forces, with no position along the span for mento to tell that it is in
    one of those areas.
    """
    return min(3 * section.height, 400 * mm)


EN_1992_2004 = register(
    DesignCode(
        title="EN 1992-2004",
        year=2004,
        materials=(Concrete_EN_1992_2004,),
        check_shear=_check_shear_EN_1992_2004,
        check_flexure=_check_flexure_EN_1992_2004,
        apply_shear_state=apply_en_shear_state,
        apply_flexure_state=apply_en_flexure_state,
        design_shear=_design_shear_EN_1992_2004,
        design_flexure=_design_flexure_EN_1992_2004,
        transverse_rebar=_transverse_rebar,
        longitudinal_rebar=_longitudinal_rebar,
        initialize_attributes=_initialize_attributes,
        # EN shear walls are not implemented; requires() names the code.
        flexure_symbols=_FLEXURE_SYMBOLS,
        units_row_shear=_UNITS_ROW_SHEAR,
        units_row_flexure=_UNITS_ROW_FLEXURE,
        summary_columns=_SUMMARY_COLUMNS,
        capacity_columns=_capacity_columns,
        shear_symbols=_SHEAR_SYMBOLS,
        summary_drop_columns=_SUMMARY_DROP_COLUMNS,
        max_bar_spacing_slab=_max_bar_spacing_slab,
    )
)
