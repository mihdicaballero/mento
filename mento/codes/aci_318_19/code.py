"""ACI 318-19, and CIRSOC 201-25 which shares its formulas.

The registry entries. CIRSOC reuses every ACI hook -- the two codes publish the
same equations -- and differs only in the bar catalogue its transverse rebar
selection draws from, which is exactly what the registry lets a code override
without a fork.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from mento.codes.ACI_318_19_beam import (
    _check_flexure_ACI_318_19,
    _check_shear_ACI_318_19,
    _design_flexure_ACI_318_19,
    _design_shear_ACI_318_19,
)
from mento.codes.ACI_318_19_wall import _check_shear_ACI_318_19_wall, _design_shear_ACI_318_19_wall
from mento.codes.check_state import (
    apply_flexure_state,
    apply_shear_state,
    apply_wall_shear_state,
)
from mento.codes.registry import DesignCode, register
from mento.material import Concrete_ACI_318_19, Concrete_CIRSOC_201_25
from mento.units import cm, dimensionless, inch, kN, kNm, mm, MPa

if TYPE_CHECKING:
    from mento.beam import RectangularBeam
    from mento.rebar import Rebar

#: ACI names the demand M_u and the capacity phi*Mn.
_FLEXURE_SYMBOLS = {
    "demand": "Mu",
    "demand_top": "Mu,top",
    "demand_bot": "Mu,bot",
    "md_demand": "M_u",
    "md_capacity": r"\phi M_n",
}

_UNITS_ROW_SHEAR = {
    "Label": "",
    "Comb.": "",
    "Av,min": "cm²/m",
    "Av,req": "cm²/m",
    "Av": "cm²/m",
    "Vu": "kN",
    "Nu": "kN",
    "ØVc": "kN",
    "ØVs": "kN",
    "ØVn": "kN",
    "ØVmax": "kN",
    "Vu≤ØVmax": "",
    "Vu≤ØVn": "",
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
    "Mu": "kNm",
    "ØMn": "kNm",
    "Mu≤ØMn": "",
    "DCR": "",
}


_SHEAR_SYMBOLS = {
    "reinforcement": "A_v",
    "demand": "V_u",
    "capacity": r"\phi V_n",
    # Row the reduced capacity sits on in this code's shear detail table.
    "capacity_row": 7,
}


def _min_stirrup_diameter(concrete: Any) -> Any:
    """ACI's smallest stirrup: a #3 bar, or 10 mm in metric practice."""
    return 10 * mm if concrete.unit_system == "metric" else 3 / 8 * inch


def _min_stirrup_diameter_cirsoc(concrete: Any) -> Any:
    """CIRSOC's catalogue starts one size below ACI's."""
    return 6 * mm


_SUMMARY_COLUMNS = {
    "moment_demand": "Mu",
    "shear_demand": "Vu",
    "shear_demand_source": "Vu",
    "axial_demand": "Nu",
    "moment_capacity_top": "ØMn,top",
    "moment_capacity_bot": "ØMn,bot",
    "shear_capacity": "ØVn",
}


def _capacity_columns(section: "RectangularBeam") -> dict:
    """ACI reports a reduced capacity phi*Mn on each face."""
    return {
        "ØMn,top": round(section._phi_M_n_top.to("kN*m").magnitude, 1),
        "ØMn,bot": round(section._phi_M_n_bot.to("kN*m").magnitude, 1),
    }


def _initialize_attributes(section: "RectangularBeam") -> None:
    """The zeroed result attributes the ACI report tables read off the beam.

    Moved here from ``RectangularBeam`` verbatim: they are this code's
    compatibility layer (ADR-0001), so the code owns them and a new code brings
    its own set without editing the element.
    """
    section._phi_V_n = 0 * kN
    section._phi_V_s = 0 * kN
    section._phi_V_c = 0 * kN
    section._phi_V_max = 0 * kN
    section._V_u = 0 * kN
    section._M_u = 0 * kNm
    section._M_u_bot = 0 * kNm
    section._M_u_top = 0 * kNm
    section._N_u = 0 * kN
    section._A_cv = 0 * cm**2
    section._k_c_min = 0 * MPa
    section._sigma_Nu = 0 * MPa
    section.V_c = 0 * kN
    section._rho_w = 0 * dimensionless
    section._lambda_s = 0
    section.f_yt = 0 * MPa
    section._max_shear_ok = False
    section._A_s_min_bot = 0 * cm**2
    section._A_s_min_top = 0 * cm**2
    section._A_s_max_bot = 0 * cm**2
    section._A_s_max_top = 0 * cm**2
    section._phi_M_n_bot = 0 * kNm
    section._phi_M_n_top = 0 * kNm
    section._d_b_max_bot = 0 * mm
    section._d_b_max_top = 0 * mm
    section.flexure_design_results_bot = None
    section.flexure_design_results_top = None
    section._A_s_bool_bot = False
    section._A_s_bool_top = False


def _transverse_rebar_aci(rebar: "Rebar", V_s_req: Any, alpha: float) -> Any:
    return rebar.transverse_rebar_ACI_318_19(V_s_req)


def _transverse_rebar_cirsoc(rebar: "Rebar", V_s_req: Any, alpha: float) -> Any:
    return rebar.transverse_rebar_CIRSOC_201_25(V_s_req)


def _longitudinal_rebar_aci(rebar: "Rebar", A_s_req: Any, A_s_max: Any, mech_cover: Any) -> Any:
    return rebar.longitudinal_rebar_ACI_318_19(A_s_req, A_s_max, mech_cover)


_COMMON = dict(
    check_shear=_check_shear_ACI_318_19,
    check_flexure=_check_flexure_ACI_318_19,
    apply_shear_state=apply_shear_state,
    apply_flexure_state=apply_flexure_state,
    design_shear=_design_shear_ACI_318_19,
    design_flexure=_design_flexure_ACI_318_19,
    longitudinal_rebar=_longitudinal_rebar_aci,
    initialize_attributes=_initialize_attributes,
    check_shear_wall=_check_shear_ACI_318_19_wall,
    # CIRSOC 201-25 reuses the ACI wall design; its bar catalogue is
    # selected inside the design itself.
    design_shear_wall=_design_shear_ACI_318_19_wall,
    apply_wall_shear_state=apply_wall_shear_state,
    flexure_symbols=_FLEXURE_SYMBOLS,
    units_row_shear=_UNITS_ROW_SHEAR,
    units_row_flexure=_UNITS_ROW_FLEXURE,
    summary_columns=_SUMMARY_COLUMNS,
    capacity_columns=_capacity_columns,
    shear_symbols=_SHEAR_SYMBOLS,
)

ACI_318_19 = register(
    DesignCode(
        title="ACI 318-19",
        year=2019,
        materials=(Concrete_ACI_318_19,),
        transverse_rebar=_transverse_rebar_aci,
        min_stirrup_diameter=_min_stirrup_diameter,
        **_COMMON,  # type: ignore[arg-type]
    )
)

CIRSOC_201_25 = register(
    DesignCode(
        title="CIRSOC 201-25",
        year=2025,
        materials=(Concrete_CIRSOC_201_25,),
        # The one thing CIRSOC does differently: its own bar catalogue.
        transverse_rebar=_transverse_rebar_cirsoc,
        min_stirrup_diameter=_min_stirrup_diameter_cirsoc,
        **_COMMON,  # type: ignore[arg-type]
    )
)
