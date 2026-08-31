"""Structural wall equations — ACI 318-19, Chapter 11.

Pure functions of floats; see the package docstring for the unit convention.

The f_yt cap of §20.2.2.4 is not repeated here: walls and beams share that
clause, so the wall checker calls
:func:`mento.codes.aci_318_19.equations.shear.max_yield_strength_for_shear`.
"""

from __future__ import annotations

import math

__all__ = [
    "alpha_c",
    "concrete_shear_stress",
    "max_shear_stress",
    "reinforcement_shear_stress",
    "min_vertical_reinforcement_ratio",
    "max_horizontal_spacing",
    "max_vertical_spacing",
]

# §11.6.1: the floor both directions share.
MIN_REINFORCEMENT_RATIO = 0.0025


def alpha_c(hw_lw: float, *, is_imperial: bool = False) -> float:
    """Coefficient alpha_c for wall shear — ACI 318-19 §11.5.4.6.

    Squat walls carry more shear in the concrete than slender ones, so the
    coefficient falls linearly across the transition 1.5 <= hw/lw <= 2.0 and is
    flat outside it.

    Args:
        hw_lw: Wall height-to-length ratio.

    Returns:
        alpha_c: 0.25 to 0.17 in SI, 3.0 to 2.0 in US customary.
    """
    alpha_hi, alpha_lo = (3.0, 2.0) if is_imperial else (0.25, 0.17)

    if hw_lw <= 1.5:
        return alpha_hi
    if hw_lw >= 2.0:
        return alpha_lo
    return alpha_hi + (hw_lw - 1.5) / 0.5 * (alpha_lo - alpha_hi)


def concrete_shear_stress(f_c: float, alpha_c_value: float, lambda_factor: float) -> float:
    """Concrete shear stress carried by a wall — ACI 318-19 §11.5.4.6.

    ``V_c = alpha_c * lambda * sqrt(f_c) * A_cv``; this returns the stress, so
    the caller multiplies by A_cv.

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
    """Upper limit on total wall shear stress — ACI 318-19 §11.5.4.6.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        lambda_factor: Lightweight concrete factor lambda.

    Returns:
        The stress that, times A_cv, gives V_n,max (MPa, or psi).
    """
    coeff = 8.0 if is_imperial else 0.66
    return coeff * lambda_factor * math.sqrt(f_c)


def reinforcement_shear_stress(rho_t: float, f_yt: float) -> float:
    """Shear stress carried by the horizontal reinforcement — ACI 318-19 §11.5.4.3.

    ``V_s = rho_t * f_yt * A_cv``; this returns the stress.

    Args:
        rho_t: Transverse (horizontal) reinforcement ratio.
        f_yt: Yield strength used for shear (MPa, or psi).

    Returns:
        The stress (MPa, or psi).
    """
    return rho_t * f_yt


def min_vertical_reinforcement_ratio(hw_lw: float, rho_t_req: float) -> float:
    """Minimum vertical reinforcement ratio — ACI 318-19 Eq. (11.6.2).

        rho_l >= 0.0025 + 0.5*(2.5 - hw/lw)*(rho_t - 0.0025)

    Args:
        hw_lw: Wall height-to-length ratio, clamped to [0.5, 2.5] by the clause.
            Above 2.5 only the floor applies; at or below 0.5 the vertical ratio
            matches the horizontal one.
        rho_t_req: Transverse reinforcement ratio required for strength.

    Returns:
        rho_l,min, never below the 0.0025 floor of §11.6.1.
    """
    r = max(0.5, min(hw_lw, 2.5))
    rho_l_eq = MIN_REINFORCEMENT_RATIO + 0.5 * (2.5 - r) * (rho_t_req - MIN_REINFORCEMENT_RATIO)
    return max(MIN_REINFORCEMENT_RATIO, rho_l_eq)


def max_horizontal_spacing(l_w: float, thickness: float, *, is_imperial: bool = False) -> float:
    """Maximum horizontal bar spacing — ACI 318-19 §11.7.3.

    Args:
        l_w: Wall length (mm, or in).
        thickness: Wall thickness (mm, or in).

    Returns:
        min(l_w/5, 3*t, 450 mm or 18 in).
    """
    absolute_cap = 18.0 if is_imperial else 450.0
    return min(l_w / 5, 3 * thickness, absolute_cap)


def max_vertical_spacing(l_w: float, thickness: float, *, is_imperial: bool = False) -> float:
    """Maximum vertical bar spacing — ACI 318-19 §11.7.3.

    Args:
        l_w: Wall length (mm, or in).
        thickness: Wall thickness (mm, or in).

    Returns:
        min(l_w/3, 3*t, 450 mm or 18 in). Looser than the horizontal limit
        because the horizontal bars are the ones resisting in-plane shear.
    """
    absolute_cap = 18.0 if is_imperial else 450.0
    return min(l_w / 3, 3 * thickness, absolute_cap)
