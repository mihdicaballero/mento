"""Structural wall equations — ACI 318-19, Chapter 11.

CIRSOC 201-25 reprints this chapter under the same numbering, so every function
below carries both citations. The divisor of Eq. (11.5.4.4) differs:
3.45*A_g in ACI 318-19 against 3.5*A_g in CIRSOC 201-25. It is supplied by
the code registry to :func:`alpha_c_in_tension`.

Pure functions of floats; see the package docstring for the unit convention.

The f_yt cap is not repeated here: walls and beams share that clause, so the
wall checker calls
:func:`mento.codes.aci_318_19.equations.shear.max_yield_strength_for_shear`.
Each code reaches the same 420 MPa by a different route, which that function
documents — ACI 318-19 §22.5.3.3 → Table 20.2.2.4(a), CIRSOC 201-25 §22.5.3.3 →
§20.2.1.3 with C 22.5.3.3.

Scope: ordinary structural walls. Both codes send a special structural wall
elsewhere — ACI 318-19 §11.1.2 to Chapter 18 (§18.10), CIRSOC 201-25 §11.1.2 to
INPRES-CIRSOC 103 Parte II-2026 — and mento implements none of that.
"""

from __future__ import annotations

import math

__all__ = [
    "alpha_c",
    "alpha_c_in_tension",
    "concrete_shear_stress",
    "max_shear_stress",
    "reinforcement_shear_stress",
    "min_vertical_reinforcement_ratio",
    "max_horizontal_spacing",
    "max_vertical_spacing",
]

# ACI 318-19 §11.6.2 / CIRSOC 201-25 §11.6.2: the floor both directions share —
# (a) for rho_l, (b) for rho_t. The same 0.0025 in both codes.
MIN_REINFORCEMENT_RATIO = 0.0025


def alpha_c(hw_lw: float, *, is_imperial: bool = False) -> float:
    """Coefficient alpha_c for wall shear — ACI 318-19 §11.5.4.3 / CIRSOC 201-25 §11.5.4.3.

    Not an equation of its own: both codes define alpha_c in the three lines
    printed under Eq. (11.5.4.3), with the same numbers — 0.25 for hw/lw <= 1.5,
    0.17 for hw/lw >= 2.0 (3 and 2 in psi), varying linearly between. Squat walls
    carry more shear in the concrete than slender ones, so the coefficient falls
    across the transition and is flat outside it.

    CIRSOC 201-25 §11.5.4.4 differs from ACI 318-19 §11.5.4.4 in the alpha_c a
    wall in net axial tension takes: 0.17*(1 + N_u/(3.5*A_g)) against ACI's
    0.17*(1 + N_u/(3.45*A_g)), both floored at zero, N_u negative in tension
    (2*(1 + N_u/(500*A_g)) in ACI's in-lb edition). That branch is evaluated
    by :func:`alpha_c_in_tension`; this function covers walls without net tension.

    Args:
        hw_lw: Wall height-to-length ratio h_w/l_w. h_w is the height of the
            whole wall, or the clear height of the wall segment or pier
            considered (Chapter 2 of both codes), not the storey height of a
            multi-storey wall.

    Returns:
        alpha_c: 0.25 to 0.17 in SI, 3.0 to 2.0 in US customary.
    """
    alpha_hi, alpha_lo = (3.0, 2.0) if is_imperial else (0.25, 0.17)

    if hw_lw <= 1.5:
        return alpha_hi
    if hw_lw >= 2.0:
        return alpha_lo
    return alpha_hi + (hw_lw - 1.5) / 0.5 * (alpha_lo - alpha_hi)


def alpha_c_in_tension(N_u: float, A_g: float, divisor: float, *, is_imperial: bool = False) -> float:
    """ACI 318-19 / CIRSOC 201-25 Eq. (11.5.4.4): alpha_c under net tension.

    N_u is negative in tension (N or lb); A_g is the gross area (mm² or in²).
    The divisor is 3.45 MPa for ACI SI, 3.5 MPa for CIRSOC, or 500 psi.
    """
    coefficient = 2.0 if is_imperial else 0.17
    return max(coefficient * (1.0 + N_u / (divisor * A_g)), 0.0)


def concrete_shear_stress(f_c: float, alpha_c_value: float, lambda_factor: float) -> float:
    """Concrete shear stress carried by a wall — ACI 318-19 Eq. (11.5.4.3) / CIRSOC 201-25 ec. (11.5.4.3).

    The concrete term of ``V_n = (alpha_c*lambda*sqrt(f_c) + rho_t*f_yt)*A_cv``,
    printed the same way in both codes, lambda included; this returns the
    stress, so the caller multiplies by A_cv.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        alpha_c_value: From :func:`alpha_c`.
        lambda_factor: Lightweight concrete factor lambda.

    Returns:
        The stress (MPa, or psi). Both unit systems share this expression — the
        system-dependent constants live in ``alpha_c``.
    """
    return alpha_c_value * lambda_factor * math.sqrt(f_c)


def max_shear_stress(f_c: float, lambda_factor: float, *, is_imperial: bool = False) -> float:
    """Upper limit on total wall shear stress — ACI 318-19 §11.5.4.2 / CIRSOC 201-25 §11.5.4.2.

    Both codes write it as V_n <= 0.66*sqrt(f_c)*A_cv (8*sqrt(f_c)*A_cv in psi),
    with no lambda — unlike Eq. (11.5.4.3) one clause below, which does carry
    one. The ``lambda_factor`` below is mento's, and it changes nothing while
    lambda = 1, the only value mento supports today.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        lambda_factor: Lightweight concrete factor lambda.

    Returns:
        The stress that, times A_cv, gives V_n,max (MPa, or psi).
    """
    coeff = 8.0 if is_imperial else 0.66
    return coeff * lambda_factor * math.sqrt(f_c)


def reinforcement_shear_stress(rho_t: float, f_yt: float) -> float:
    """Stress carried by the horizontal reinforcement — ACI 318-19 Eq. (11.5.4.3) / CIRSOC 201-25 ec. (11.5.4.3).

    The reinforcement term of the same equation, ``rho_t*f_yt*A_cv``, identical
    in both codes; this returns the stress.

    Args:
        rho_t: Transverse (horizontal) reinforcement ratio.
        f_yt: Yield strength used for shear (MPa, or psi).

    Returns:
        The stress (MPa, or psi).
    """
    return rho_t * f_yt


def min_vertical_reinforcement_ratio(hw_lw: float, rho_t: float, rho_t_req: float) -> float:
    """Minimum vertical reinforcement ratio — ACI 318-19 §11.6.2(a) / CIRSOC 201-25 §11.6.2(a).

        rho_l >= 0.0025 + 0.5*(2.5 - hw/lw)*(rho_t - 0.0025)        Eq. (11.6.2)

    §11.6.2(a) of both codes, in three parts: rho_l shall be at least the
    greater of Eq. (11.6.2) and 0.0025, but need not exceed the rho_t required
    for strength by §11.5.4.3. The rho_t of the equation is the plain ratio of
    Chapter 2 — the one the wall provides. The clause names the required one
    apart, as a ceiling, and that ceiling can only bind if the equation is fed
    the provided ratio: 0.5*(2.5 - hw/lw) is at most 1 over the clamped range,
    so with the required ratio as its argument the equation could never return
    more than it, and the "need not exceed" would be dead text. Read literally,
    then, a wall that carries more horizontal steel than its shear needs has
    to carry more vertical steel too, up to what the shear needed:

        rho_l,min = max(0.0025, min(eq(rho_t), rho_t_req))

    That is the literal reading and, for a wall whose mesh meets its shear, the
    conservative one — eq(rho_t) >= eq(rho_t_req) whenever rho_t >= rho_t_req.
    It is a reading of the clause and not a printed formula, and the wall
    docstrings say so. The 0.0025 floor is never lifted: a ceiling below it
    leaves the floor.

    Not implemented, in either code: the §11.6.1 branch with Table 11.6.1, which
    a wall whose in-plane V_u stays below 0.04*phi*alpha_c*lambda*sqrt(f_c)*A_cv
    may use instead and which allows less. Always applying §11.6.2 is the
    conservative reading of the two.

    Args:
        hw_lw: Wall height-to-length ratio, clamped to [0.5, 2.5] by the clause.
            Above 2.5 only the floor applies; at or below 0.5 the equation
            returns the horizontal ratio itself.
        rho_t: Transverse (horizontal) reinforcement ratio the wall provides —
            the rho_t printed in Eq. (11.6.2).
        rho_t_req: Transverse reinforcement ratio required for strength by
            §11.5.4.3, the ceiling of §11.6.2(a). The caller passes it already
            floored at the 0.0025 of §11.6.2(b); an unfloored one changes
            nothing, since the floor below wins.

    Returns:
        rho_l,min: the equation with the provided rho_t, capped by rho_t_req,
        never below the 0.0025 floor of §11.6.2(a) in both codes.
    """
    r = max(0.5, min(hw_lw, 2.5))
    rho_l_eq = MIN_REINFORCEMENT_RATIO + 0.5 * (2.5 - r) * (rho_t - MIN_REINFORCEMENT_RATIO)
    return max(MIN_REINFORCEMENT_RATIO, min(rho_l_eq, rho_t_req))


def max_horizontal_spacing(l_w: float, thickness: float, *, is_imperial: bool = False) -> float:
    """Maximum horizontal bar spacing — ACI 318-19 §11.7.3.1 / CIRSOC 201-25 §11.7.3.1.

    The spacing of the transverse reinforcement in a cast-in-place wall: both
    codes cap it at the lesser of 3h and 450 mm (18 in.), and add that s shall
    not exceed l_w/5 *where shear reinforcement is required for in-plane
    strength*. mento applies the l_w/5 always, which is the conservative
    reading. Precast walls (§11.7.3.2 in both) are out of scope.

    Args:
        l_w: Wall length (mm, or in).
        thickness: Wall thickness — the h of the clause (mm, or in).

    Returns:
        min(l_w/5, 3*t, 450 mm or 18 in).
    """
    absolute_cap = 18.0 if is_imperial else 450.0
    return min(l_w / 5, 3 * thickness, absolute_cap)


def max_vertical_spacing(l_w: float, thickness: float, *, is_imperial: bool = False) -> float:
    """Maximum vertical bar spacing — ACI 318-19 §11.7.2.1 / CIRSOC 201-25 §11.7.2.1.

    The spacing of the longitudinal reinforcement in a cast-in-place wall: both
    codes cap it at the lesser of 3h and 450 mm (18 in.), and add that s shall
    not exceed l_w/3 *where shear reinforcement is required for in-plane
    strength*. As in :func:`max_horizontal_spacing`, mento applies the l_w/3
    always. Precast walls (§11.7.2.2 in both) are out of scope.

    Args:
        l_w: Wall length (mm, or in).
        thickness: Wall thickness — the h of the clause (mm, or in).

    Returns:
        min(l_w/3, 3*t, 450 mm or 18 in). Looser than the horizontal limit
        because the horizontal bars are the ones resisting in-plane shear.
    """
    absolute_cap = 18.0 if is_imperial else 450.0
    return min(l_w / 3, 3 * thickness, absolute_cap)
