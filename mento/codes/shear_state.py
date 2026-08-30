"""The values one shear check produces, held off the section.

A check used to write these onto the beam, which made checking a section a
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

from mento.units import cm, inch, kip, kN, m, psi, MPa, ft, dimensionless

if TYPE_CHECKING:
    from mento.beam import RectangularBeam

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
