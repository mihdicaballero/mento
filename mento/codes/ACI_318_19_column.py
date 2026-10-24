"""One-way shear of a column section — ACI 318-19 Chapter 10 with §22.5, and CIRSOC 201-25 with it.

The checker of ADR-0002 for a column: it reads the section's float view
(:class:`~mento.column_section.ColumnShearSection`) and the combination once,
calls the float equations of :mod:`mento.codes.aci_318_19.equations.shear`
and returns a :class:`~mento.column_results.ColumnShearCheck`. It writes
nothing back onto the column.

The resistance is the beam's -- V_c of Table 22.5.5.1 with N_u/(6·A_g), V_s of
§22.5.8.5.3, the section limit of §22.5.1.2 -- on the column's own section,
and it gives the same numbers as :mod:`mento.codes.ACI_318_19_beam` on the
same rectangle (``tests/elements/test_column.py``). What the column reads
from its own chapter, the same clauses in both codes:

* §10.6.2.1 -- shear reinforcement is required where V_u > 0.5·φ·V_c, with
  the axial load in V_c, in place of the beam threshold of §9.6.3.1.
* §10.6.2.2 -- where it is, at least A_v,min of Table 10.6.2.2, the
  expressions of Table 9.6.3.4.
* Table 10.7.6.5.2 -- the spacing of the legs along and across, the rows and
  the threshold of Table 9.7.6.2.2, with the caps each code prints there
  (``stirrup_spacing_caps``).
* §22.5.3.2 lifts the 8.3 MPa (100 psi) cap on √f'c of §22.5.3.1 for beams
  and joists only, so a column keeps it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Tuple, cast

from mento.codes.aci_318_19.equations import shear as shear_eq
from mento.codes.check_state import _scaled_to_display
from mento.codes.registry import design_code
from mento.column_results import ColumnShearCheck
from mento.precompute import CANONICAL

if TYPE_CHECKING:
    from mento.column_section import ColumnSection, ColumnShearSection, Direction
    from mento.forces import Forces
    from mento.material import Concrete_ACI_318_19


def _concrete_shear_strength(
    sec: "ColumnShearSection", lambda_factor: float, sigma_Nu: float, *, has_min_rebar: bool
) -> float:
    """V_c of Table 22.5.5.1 for the row ``has_min_rebar`` selects, capped by §22.5.5.1.1 and never negative.

    Rows (a)/(b) with A_v >= A_v,min, row (c) with the size effect without it.
    ``allow_high_strength=False``: §22.5.3.2 does not reach a column.
    """
    A_cv = sec.b_w * sec.d
    rho_w = sec.A_s / A_cv
    lambda_s = shear_eq.size_effect_factor(sec.d, is_imperial=sec.is_imperial)
    v_c = shear_eq.concrete_shear_stress(
        sec.f_c,
        lambda_factor,
        rho_w,
        sigma_Nu,
        lambda_s,
        has_min_rebar=has_min_rebar,
        allow_high_strength=False,
        is_imperial=sec.is_imperial,
    )
    v_c_max = shear_eq.max_concrete_shear_stress(
        sec.f_c, lambda_factor, has_min_rebar=has_min_rebar, allow_high_strength=False, is_imperial=sec.is_imperial
    )
    return min(v_c_max * A_cv, max(0.0, v_c * A_cv))


def _spacing_limits(column: "ColumnSection", sec: "ColumnShearSection", V_s_req: float) -> Tuple[float, float, bool]:
    """``(s_max_l, s_max_w, halved)`` of Table 10.7.6.5.2, in the canonical length."""
    length = CANONICAL[sec.is_imperial]["length"]
    caps = design_code(column.concrete).requires("stirrup_spacing_caps")(column.concrete)
    cap_low, cap_high = (float(cap.to(length).magnitude) for cap in caps)
    A_cv = sec.b_w * sec.d
    s_max_l, s_max_w = shear_eq.max_stirrup_spacing(
        V_s_req, sec.f_c, A_cv, sec.d, cap_low, cap_high, is_imperial=sec.is_imperial
    )
    halved = shear_eq.stirrup_spacing_halved(V_s_req, sec.f_c, A_cv, is_imperial=sec.is_imperial)
    return s_max_l, s_max_w, halved


def check_column_shear_ACI_318_19(column: "ColumnSection", force: "Forces", direction: "Direction") -> ColumnShearCheck:
    """The shear of one combination in one direction of a column.

    φ·V_n = φ·(V_c + V_s) >= V_u, ACI 318-19 / CIRSOC 201-25 §22.5.1.1 with
    §10.5.1.1, never above φ·V_max of §22.5.1.2. φ = 0.75, Table 21.2.1(b).
    The DCR is V_u over the smaller of the two, infinite where both are zero.
    """
    sec = column.shear_section(direction)
    imperial = sec.is_imperial
    force_unit = CANONICAL[imperial]["force"]
    concrete = cast("Concrete_ACI_318_19", column.concrete)
    phi, lambda_factor = concrete.phi_v, concrete.lambda_factor

    N_u = float(force._N_x.to(force_unit).magnitude)
    V = force._V_y if direction == "y" else force._V_z
    V_u = abs(float(V.to(force_unit).magnitude))
    f_yt = shear_eq.max_yield_strength_for_shear(sec.f_y, is_imperial=imperial)
    A_cv = sec.b_w * sec.d

    # Table 22.5.5.1 picks its row on the A_v,min it defines, whether or not
    # the demand asks for it (R22.5.5.1).
    A_v_min_table = shear_eq.min_shear_reinforcement_ratio(sec.f_c, f_yt, sec.b_w, is_imperial=imperial)
    has_min_rebar = sec.A_v >= A_v_min_table
    sigma_Nu = shear_eq.axial_stress_influence(N_u, sec.A_g, sec.f_c)
    V_c = _concrete_shear_strength(sec, lambda_factor, sigma_Nu, has_min_rebar=has_min_rebar)
    phi_V_c = phi * V_c

    required = V_u > shear_eq.column_shear_reinforcement_threshold(phi_V_c)
    A_v_min = A_v_min_table if required else 0.0
    V_s_req = max((V_u - phi_V_c) / phi, 0.0)
    A_v_req = max(V_s_req / (f_yt * sec.d), A_v_min) if required else 0.0

    increment = shear_eq.shear_stress_capacity_increment(sec.f_c, lambda_factor, is_imperial=imperial) * A_cv
    phi_V_max = phi * (V_c + increment)
    V_c_reinforced = (
        V_c if has_min_rebar else _concrete_shear_strength(sec, lambda_factor, sigma_Nu, has_min_rebar=True)
    )
    phi_V_n = phi * (V_c + shear_eq.shear_strength_of_reinforcement(sec.A_v, f_yt, sec.d))
    capacity = min(phi_V_n, phi_V_max)
    DCR = V_u / capacity if capacity > 0 else float("inf")

    s_max_l: Optional[float] = None
    s_max_w: Optional[float] = None
    halved: Optional[bool] = None
    if required:
        s_max_l, s_max_w, halved = _spacing_limits(column, sec, V_s_req)

    return ColumnShearCheck(
        label=force.label,
        direction=direction,
        N_demand=_q(N_u, "force", imperial),
        V_demand=_q(V_u, "force", imperial),
        b_w=_q(sec.b_w, "length", imperial),
        d=_q(sec.d, "length", imperial),
        V_concrete=_q(phi_V_c, "force", imperial),
        reinforcement_required=required,
        A_v=_q(sec.A_v, "per_length", imperial),
        A_v_req=_q(A_v_req, "per_length", imperial),
        A_v_min=_q(A_v_min, "per_length", imperial),
        V_capacity=_q(capacity, "force", imperial),
        V_max=_q(phi * (V_c_reinforced + increment), "force", imperial),
        DCR=DCR,
        s_l=_length_or_none(sec.s_l, imperial),
        s_w=_length_or_none(sec.s_w, imperial),
        s_max_l=_length_or_none(s_max_l, imperial),
        # No limit across where there is no second leg to hold to it: a circle, or no ties.
        s_max_w=_length_or_none(s_max_w, imperial) if sec.s_w else None,
        spacing_halved=halved,
    )


def _q(value: float, kind: str, imperial: bool) -> Any:
    """A canonical float as the quantity the result carries."""
    return _scaled_to_display(value, kind, imperial)


def _length_or_none(value: Optional[float], imperial: bool) -> Any:
    """A spacing as a quantity, or ``None`` where there is none."""
    return _scaled_to_display(value, "length", imperial) if value else None
