"""The values one check produces, held off the section.

A check used to write these onto the element, which made checking a section a
mutation of it and left the caller reading whichever combination happened to
run last. They live here instead: the design-code helpers fill a state, the
check returns it, and only the reporting path copies it back onto the element
as the compatibility layer ADR-0001 keeps for one release cycle.

The fields are pre-zeroed in the section's own unit system by
:func:`new_shear_state`, so nothing is ``Optional`` and no reader has to narrow.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from pint import Quantity

from mento.units import cm, inch, kip, kN, kNm, m, mm, psi, MPa, ft, dimensionless

if TYPE_CHECKING:
    from mento.beam import RectangularBeam
    from mento.shear_wall import ShearWall

#: The attributes a state maps onto, in the order the dataclass declares them.
#: The compatibility layer and the tests that pin it both read this.
BEAM_ATTRIBUTES = {
    "V_u": "_V_u",
    "N_u": "_N_u",
    "f_yt": "f_yt",
    "A_s_tension": "_A_s_tension",
    "A_cv": "_A_cv",
    "rho_w": "_rho_w",
    "lambda_s": "_lambda_s",
    "sigma_Nu": "_sigma_Nu",
    "k_c_min": "_k_c_min",
    "V_c": "V_c",
    "phi_V_c": "_phi_V_c",
    "phi_V_s": "_phi_V_s",
    "phi_V_n": "_phi_V_n",
    "phi_V_max": "_phi_V_max",
    "V_s_req": "_V_s_req",
    "A_v_min": "_A_v_min",
    "A_v_req": "_A_v_req",
    "stirrup_s_w": "_stirrup_s_w",
    "stirrup_s_max_l": "_stirrup_s_max_l",
    "stirrup_s_max_w": "_stirrup_s_max_w",
    "max_shear_ok": "_max_shear_ok",
    "DCR": "_DCRv",
}


@dataclass
class ShearCheckState:
    """One combination's shear result. Mutable while the check fills it in."""

    V_u: Quantity
    N_u: Quantity
    f_yt: Quantity
    A_s_tension: Quantity
    A_cv: Quantity
    rho_w: Quantity
    lambda_s: float
    sigma_Nu: Quantity
    k_c_min: Quantity
    V_c: Quantity
    phi_V_c: Quantity
    phi_V_s: Quantity
    phi_V_n: Quantity
    phi_V_max: Quantity
    V_s_req: Quantity
    A_v_min: Quantity
    A_v_req: Quantity
    stirrup_s_w: Quantity
    stirrup_s_max_l: Quantity
    stirrup_s_max_w: Quantity
    max_shear_ok: bool
    DCR: float


def new_shear_state(section: "RectangularBeam") -> ShearCheckState:
    """A zeroed state carrying the section's unit system."""
    imperial = section.concrete.is_imperial
    force = 0 * (kip if imperial else kN)
    stress = 0 * (psi if imperial else MPa)
    length = 0 * (inch if imperial else cm)
    per_length = 0 * (inch**2 / ft if imperial else cm**2 / m)
    area = 0 * (inch**2 if imperial else cm**2)

    return ShearCheckState(
        V_u=force,
        N_u=force,
        f_yt=stress,
        A_s_tension=area,
        A_cv=area,
        rho_w=0 * dimensionless,
        lambda_s=0.0,
        sigma_Nu=stress,
        k_c_min=stress,
        V_c=force,
        phi_V_c=force,
        phi_V_s=force,
        phi_V_n=force,
        phi_V_max=force,
        V_s_req=force,
        A_v_min=per_length,
        A_v_req=per_length,
        stirrup_s_w=length,
        stirrup_s_max_l=length,
        stirrup_s_max_w=length,
        max_shear_ok=False,
        DCR=0.0,
    )


def apply_shear_state(section: "RectangularBeam", state: ShearCheckState) -> None:
    """Copy a state onto the section.

    The compatibility layer: the report tables and a handful of tests still read
    these attributes off the element. The values-only entry points skip this, so
    a loop over many sections leaves every one of them untouched.
    """
    for field_name, attribute in BEAM_ATTRIBUTES.items():
        setattr(section, attribute, getattr(state, field_name))


#: EN 1992-1-1 works with a different set of quantities than ACI, so it gets its
#: own state rather than a union of both.
EN_BEAM_ATTRIBUTES = {
    "N_Ed": "_N_Ed",
    "V_Ed_1": "_V_Ed_1",
    "V_Ed_2": "_V_Ed_2",
    "f_cd": "_f_cd",
    "f_cd_shear": "_f_cd_shear",
    "f_ywd": "_f_ywd",
    "A_s_tension": "_A_s_tension",
    "A_v_min": "_A_v_min",
    "A_v_req": "_A_v_req",
    "rho_l_bot": "_rho_l_bot",
    "rho_l_top": "_rho_l_top",
    "k_value": "_k_value",
    "sigma_cp": "_sigma_cp",
    "theta": "_theta",
    "cot_theta": "_cot_theta",
    "z": "_z",
    "V_Rd_c": "_V_Rd_c",
    "V_Rd_s": "_V_Rd_s",
    "V_Rd": "_V_Rd",
    "V_Rd_max": "_V_Rd_max",
    "V_s_req": "_V_s_req",
    "stirrup_s_w": "_stirrup_s_w",
    "stirrup_s_max_l": "_stirrup_s_max_l",
    "stirrup_s_max_w": "_stirrup_s_max_w",
    "max_shear_ok": "_max_shear_ok",
    "DCR": "_DCRv",
}


@dataclass
class ENShearCheckState:
    """One combination's EN 1992-1-1 shear result."""

    N_Ed: Quantity
    V_Ed_1: Quantity
    V_Ed_2: Quantity
    f_cd: Quantity
    f_cd_shear: Quantity
    f_ywd: Quantity
    A_s_tension: Quantity
    A_v_min: Quantity
    A_v_req: Quantity
    rho_l_bot: Quantity
    rho_l_top: Quantity
    k_value: float
    sigma_cp: Quantity
    theta: float
    cot_theta: float
    z: Quantity
    V_Rd_c: Quantity
    V_Rd_s: Quantity
    V_Rd: Quantity
    V_Rd_max: Quantity
    V_s_req: Quantity
    stirrup_s_w: Quantity
    stirrup_s_max_l: Quantity
    stirrup_s_max_w: Quantity
    max_shear_ok: bool
    DCR: float


def new_en_shear_state(section: "RectangularBeam") -> ENShearCheckState:
    """A zeroed EN state. Eurocode sections are metric only."""
    force = 0 * kN
    stress = 0 * MPa
    length = 0 * cm
    per_length = 0 * cm**2 / m
    area = 0 * cm**2
    return ENShearCheckState(
        N_Ed=force,
        V_Ed_1=force,
        V_Ed_2=force,
        f_cd=stress,
        f_cd_shear=stress,
        f_ywd=stress,
        A_s_tension=area,
        A_v_min=per_length,
        A_v_req=per_length,
        rho_l_bot=0 * dimensionless,
        rho_l_top=0 * dimensionless,
        k_value=0.0,
        sigma_cp=stress,
        theta=0.0,
        cot_theta=0.0,
        z=length,
        V_Rd_c=force,
        V_Rd_s=force,
        V_Rd=force,
        V_Rd_max=force,
        V_s_req=force,
        stirrup_s_w=length,
        stirrup_s_max_l=length,
        stirrup_s_max_w=length,
        max_shear_ok=False,
        DCR=0.0,
    )


def apply_en_shear_state(section: "RectangularBeam", state: ENShearCheckState) -> None:
    """Copy an EN state onto the section — the same compatibility layer."""
    for field_name, attribute in EN_BEAM_ATTRIBUTES.items():
        setattr(section, attribute, getattr(state, field_name))


#: A structural wall reports different quantities again — rho_t/rho_l rather
#: than A_v, and one capacity per curtain direction.
WALL_BEAM_ATTRIBUTES = {
    "V_u": "_V_u",
    "N_u": "_N_u",
    "Acv": "_Acv",
    "alpha_c": "_alpha_c",
    "hw_lw": "_hw_lw",
    "f_yt_wall": "_f_yt_wall",
    "V_c_wall": "_V_c_wall",
    "V_s_wall": "_V_s_wall",
    "V_n_wall": "_V_n_wall",
    "V_n_max": "_V_n_max",
    "phi_V_n_wall": "_phi_V_n_wall",
    "phi_V_n_max_wall": "_phi_V_n_max_wall",
    "rho_t_min": "_rho_t_min",
    "rho_l_min": "_rho_l_min",
    "rho_t_req": "_rho_t_req",
    "s_h_max": "_s_h_max",
    "s_v_max": "_s_v_max",
    "DCR": "_DCRv_wall",
}


@dataclass
class WallShearCheckState:
    """One combination's ACI 318-19 structural wall shear result."""

    V_u: Quantity
    N_u: Quantity
    Acv: Quantity
    alpha_c: float
    hw_lw: float
    f_yt_wall: Quantity
    V_c_wall: Quantity
    V_s_wall: Quantity
    V_n_wall: Quantity
    V_n_max: Quantity
    phi_V_n_wall: Quantity
    phi_V_n_max_wall: Quantity
    rho_t_min: Quantity
    rho_l_min: Quantity
    rho_t_req: Quantity
    s_h_max: Quantity
    s_v_max: Quantity
    DCR: float


def new_wall_shear_state(section: "ShearWall") -> WallShearCheckState:
    """A zeroed wall state carrying the section's unit system."""
    imperial = section.concrete.is_imperial
    force = 0 * (kip if imperial else kN)
    stress = 0 * (psi if imperial else MPa)
    length = 0 * (inch if imperial else mm)
    area = 0 * (inch**2 if imperial else cm**2)
    return WallShearCheckState(
        V_u=force,
        N_u=force,
        Acv=area,
        alpha_c=0.0,
        hw_lw=0.0,
        f_yt_wall=stress,
        V_c_wall=force,
        V_s_wall=force,
        V_n_wall=force,
        V_n_max=force,
        phi_V_n_wall=force,
        phi_V_n_max_wall=force,
        rho_t_min=0 * dimensionless,
        rho_l_min=0 * dimensionless,
        rho_t_req=0 * dimensionless,
        s_h_max=length,
        s_v_max=length,
        DCR=0.0,
    )


def apply_wall_shear_state(section: "ShearWall", state: WallShearCheckState) -> None:
    """Copy a wall state onto the section — the same compatibility layer."""
    for field_name, attribute in WALL_BEAM_ATTRIBUTES.items():
        setattr(section, attribute, getattr(state, field_name))


#: Flexure reports per face, so its state nests two of them.
FLEXURE_BEAM_ATTRIBUTES = {
    "M_u": "_M_u",
    "M_u_bot": "_M_u_bot",
    "M_u_top": "_M_u_top",
    "f_yt": "f_yt",
    "A_s_tension": "_A_s_tension",
    "A_s_min_bot": "_A_s_min_bot",
    "A_s_min_top": "_A_s_min_top",
    "A_s_max_bot": "_A_s_max_bot",
    "A_s_max_top": "_A_s_max_top",
    "A_s_req_bot": "_A_s_req_bot",
    "A_s_req_top": "_A_s_req_top",
    "phi_M_n_bot": "_phi_M_n_bot",
    "phi_M_n_top": "_phi_M_n_top",
    "c_d_bot": "_c_d_bot",
    "c_d_top": "_c_d_top",
    "d_b_max_bot": "_d_b_max_bot",
    "d_b_max_top": "_d_b_max_top",
    "rho_l_bot": "_rho_l_bot",
    "rho_l_top": "_rho_l_top",
    "DCR_bot": "_DCRb_bot",
    "DCR_top": "_DCRb_top",
    "A_s_bool_bot": "_A_s_bool_bot",
    "A_s_bool_top": "_A_s_bool_top",
}


@dataclass
class FlexureCheckState:
    """One combination's flexure result, both faces."""

    M_u: Quantity
    M_u_bot: Quantity
    M_u_top: Quantity
    f_yt: Quantity
    A_s_tension: Quantity
    A_s_min_bot: Quantity
    A_s_min_top: Quantity
    A_s_max_bot: Quantity
    A_s_max_top: Quantity
    A_s_req_bot: Quantity
    A_s_req_top: Quantity
    phi_M_n_bot: Quantity
    phi_M_n_top: Quantity
    c_d_bot: float
    c_d_top: float
    d_b_max_bot: Quantity
    d_b_max_top: Quantity
    rho_l_bot: Quantity
    rho_l_top: Quantity
    DCR_bot: float
    DCR_top: float
    A_s_bool_bot: bool
    A_s_bool_top: bool
    doubly_reinforced: bool


def new_flexure_state(section: "RectangularBeam") -> FlexureCheckState:
    """A zeroed flexure state carrying the section's unit system."""
    imperial = section.concrete.is_imperial
    stress = 0 * (psi if imperial else MPa)
    length = 0 * (inch if imperial else mm)
    area = 0 * (inch**2 if imperial else cm**2)
    moment = 0 * (kip * ft if imperial else kNm)
    return FlexureCheckState(
        M_u=moment,
        M_u_bot=moment,
        M_u_top=moment,
        f_yt=stress,
        A_s_tension=area,
        A_s_min_bot=area,
        A_s_min_top=area,
        A_s_max_bot=area,
        A_s_max_top=area,
        A_s_req_bot=area,
        A_s_req_top=area,
        phi_M_n_bot=moment,
        phi_M_n_top=moment,
        c_d_bot=0.0,
        c_d_top=0.0,
        d_b_max_bot=length,
        d_b_max_top=length,
        rho_l_bot=0 * dimensionless,
        rho_l_top=0 * dimensionless,
        DCR_bot=0.0,
        DCR_top=0.0,
        A_s_bool_bot=False,
        A_s_bool_top=False,
        doubly_reinforced=False,
    )


def apply_flexure_state(section: "RectangularBeam", state: FlexureCheckState) -> None:
    """Copy a flexure state onto the section — the same compatibility layer."""
    for field_name, attribute in FLEXURE_BEAM_ATTRIBUTES.items():
        setattr(section, attribute, getattr(state, field_name))
    # Sticky rather than copied: the flag means "some combination needed
    # compression steel", so the last one checked must not clear it.
    section._doubly_reinforced = section._doubly_reinforced or state.doubly_reinforced


#: EN names its demand M_Ed and its capacity M_Rd, so flexure gets a second
#: state rather than a union of both codes' vocabularies.
EN_FLEXURE_BEAM_ATTRIBUTES = {
    "M_Ed": "_M_Ed",
    "M_Ed_bot": "_M_Ed_bot",
    "M_Ed_top": "_M_Ed_top",
    "M_Rd_bot": "_M_Rd_bot",
    "M_Rd_top": "_M_Rd_top",
    "f_cd": "_f_cd",
    "f_cd_shear": "_f_cd_shear",
    "f_ywd": "_f_ywd",
    "A_s_min_bot": "_A_s_min_bot",
    "A_s_min_top": "_A_s_min_top",
    "A_s_max_bot": "_A_s_max_bot",
    "A_s_max_top": "_A_s_max_top",
    "A_s_req_bot": "_A_s_req_bot",
    "A_s_req_top": "_A_s_req_top",
    "c_d_bot": "_c_d_bot",
    "c_d_top": "_c_d_top",
    "d_b_max_bot": "_d_b_max_bot",
    "d_b_max_top": "_d_b_max_top",
    "rho_l_bot": "_rho_l_bot",
    "rho_l_top": "_rho_l_top",
    "DCR_bot": "_DCRb_bot",
    "DCR_top": "_DCRb_top",
}


@dataclass
class ENFlexureCheckState:
    """One combination's EN 1992-1-1 flexure result, both faces."""

    M_Ed: Quantity
    M_Ed_bot: Quantity
    M_Ed_top: Quantity
    M_Rd_bot: Quantity
    M_Rd_top: Quantity
    f_cd: Quantity
    f_cd_shear: Quantity
    f_ywd: Quantity
    A_s_min_bot: Quantity
    A_s_min_top: Quantity
    A_s_max_bot: Quantity
    A_s_max_top: Quantity
    A_s_req_bot: Quantity
    A_s_req_top: Quantity
    c_d_bot: float
    c_d_top: float
    d_b_max_bot: Quantity
    d_b_max_top: Quantity
    rho_l_bot: Quantity
    rho_l_top: Quantity
    DCR_bot: float
    DCR_top: float


def new_en_flexure_state(section: "RectangularBeam") -> ENFlexureCheckState:
    """A zeroed EN flexure state. Eurocode sections are metric only."""
    stress, length, area, moment = 0 * MPa, 0 * mm, 0 * cm**2, 0 * kNm
    return ENFlexureCheckState(
        M_Ed=moment,
        M_Ed_bot=moment,
        M_Ed_top=moment,
        M_Rd_bot=moment,
        M_Rd_top=moment,
        f_cd=stress,
        f_cd_shear=stress,
        f_ywd=stress,
        A_s_min_bot=area,
        A_s_min_top=area,
        A_s_max_bot=area,
        A_s_max_top=area,
        A_s_req_bot=area,
        A_s_req_top=area,
        c_d_bot=0.0,
        c_d_top=0.0,
        d_b_max_bot=length,
        d_b_max_top=length,
        rho_l_bot=0 * dimensionless,
        rho_l_top=0 * dimensionless,
        DCR_bot=0.0,
        DCR_top=0.0,
    )


def apply_en_flexure_state(section: "RectangularBeam", state: ENFlexureCheckState) -> None:
    """Copy an EN flexure state onto the section — the same compatibility layer."""
    for field_name, attribute in EN_FLEXURE_BEAM_ATTRIBUTES.items():
        setattr(section, attribute, getattr(state, field_name))
