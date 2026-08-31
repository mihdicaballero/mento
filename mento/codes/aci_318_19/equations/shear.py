"""One-way shear equations — ACI 318-19, Chapter 22.5 and Table 9.6.3.4.

Pure functions of floats; see the package docstring for the unit convention and
for why the SI and US customary coefficients are both spelled out instead of
one being converted into the other.
"""

from __future__ import annotations

import math

__all__ = [
    "size_effect_factor",
    "axial_stress_influence",
    "concrete_shear_stress",
    "max_concrete_shear_stress",
    "shear_stress_capacity_increment",
    "min_shear_reinforcement_threshold_stress",
    "max_yield_strength_for_shear",
    "min_shear_reinforcement_ratio",
    "shear_strength_of_reinforcement",
    "max_stirrup_spacing",
]


def size_effect_factor(d: float, *, is_imperial: bool = False) -> float:
    """Size effect factor lambda_s — ACI 318-19 Eq. 22.5.5.1.3.

    Args:
        d: Effective depth for shear (mm, or in when ``is_imperial``).

    Returns:
        lambda_s, capped at 1.0.

    The cap is part of the equation, not a guard: the factor exists to reduce
    V_c in deep members, so without it a shallow member (d below 250 mm / 10 in)
    would get lambda_s > 1 and an inflated V_c — the opposite of its purpose.
    """
    if is_imperial:
        return min(math.sqrt(2 / (1 + d / 10)), 1.0)
    return min(math.sqrt(2 / (1 + 0.004 * d)), 1.0)


def axial_stress_influence(N_u: float, A_g: float, f_c: float) -> float:
    """Axial stress term sigma_Nu = N_u / (6*A_g) — ACI 318-19 §22.5.5.1.

    Args:
        N_u: Factored axial force, positive in compression (N, or lb).
        A_g: Gross section area (mm², or in²).
        f_c: Specified concrete compressive strength (MPa, or psi).

    Returns:
        The axial contribution to shear stress, capped at 0.05*f_c by the
        clause. Same expression in both unit systems.
    """
    return min(N_u / (6 * A_g), 0.05 * f_c)


def concrete_shear_stress(
    f_c: float,
    lambda_factor: float,
    rho_w: float,
    sigma_Nu: float,
    lambda_s: float,
    *,
    has_min_rebar: bool,
    is_imperial: bool = False,
) -> float:
    """Concrete shear stress v_c — ACI 318-19 Table 22.5.5.1.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        lambda_factor: Lightweight concrete factor lambda.
        rho_w: Longitudinal tension reinforcement ratio A_s / (b_w*d).
        sigma_Nu: Axial stress term from :func:`axial_stress_influence`.
        lambda_s: Size effect factor from :func:`size_effect_factor`.
        has_min_rebar: True when A_v >= A_v_min, which selects rows (a)/(b) of
            the table; False selects row (c), where the size effect factor
            applies because there are no stirrups to control the crack.

    Returns:
        The shear stress carried by the concrete (MPa, or psi).
    """
    sqrt_f_c = math.sqrt(f_c)

    if not has_min_rebar:
        # Table 22.5.5.1(c)
        coeff = 8 if is_imperial else 0.66
        return coeff * lambda_s * lambda_factor * rho_w ** (1 / 3) * sqrt_f_c + sigma_Nu

    # Rows (a) and (b): the member may take whichever is larger.
    coeff_a = 2 if is_imperial else 0.17
    coeff_b = 8 if is_imperial else 0.66
    return max(
        coeff_a * lambda_factor * sqrt_f_c + sigma_Nu,
        coeff_b * lambda_factor * rho_w ** (1 / 3) * sqrt_f_c + sigma_Nu,
    )


def max_concrete_shear_stress(f_c: float, lambda_factor: float, *, is_imperial: bool = False) -> float:
    """Upper limit on the concrete shear stress — ACI 318-19 §22.5.5.1.1.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        lambda_factor: Lightweight concrete factor lambda.

    Returns:
        The largest v_c the table may produce (MPa, or psi).
    """
    coeff = 5 if is_imperial else 0.42
    return coeff * lambda_factor * math.sqrt(f_c)


def shear_stress_capacity_increment(f_c: float, lambda_factor: float, *, is_imperial: bool = False) -> float:
    """Largest stress the stirrups may add to V_c — ACI 318-19 §22.5.1.2.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        lambda_factor: Lightweight concrete factor lambda.

    Returns:
        The stress that, times A_cv, bounds V_s and therefore caps the total
        shear the section may be designed for (MPa, or psi).
    """
    coeff = 8 if is_imperial else 0.66
    return coeff * lambda_factor * math.sqrt(f_c)


def min_shear_reinforcement_threshold_stress(f_c: float, lambda_factor: float, *, is_imperial: bool = False) -> float:
    """Stress below which no shear reinforcement is required — ACI 318-19 §9.6.3.1.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        lambda_factor: Lightweight concrete factor lambda.

    Returns:
        The stress to compare V_u/(phi_v*A_cv) against (MPa, or psi). Below it,
        the clause waives both A_v_min and A_v_req.
    """
    coeff = 1 if is_imperial else 0.083
    return coeff * lambda_factor * math.sqrt(f_c)


def max_yield_strength_for_shear(f_y: float, *, is_imperial: bool = False) -> float:
    """Yield strength usable for shear reinforcement — ACI 318-19 §20.2.2.4.

    Args:
        f_y: Specified yield strength of the reinforcement (MPa, or psi).

    Returns:
        f_yt, capped at 420 MPa / 60,000 psi.
    """
    cap = 60_000.0 if is_imperial else 420.0
    return min(f_y, cap)


def min_shear_reinforcement_ratio(f_c: float, f_yt: float, b_w: float, *, is_imperial: bool = False) -> float:
    """Minimum shear reinforcement A_v,min / s — ACI 318-19 Table 9.6.3.4.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        f_yt: Yield strength used for shear, from
            :func:`max_yield_strength_for_shear` (MPa, or psi).
        b_w: Web width (mm, or in).

    Returns:
        Reinforcement area per unit length along the member — mm²/mm, which is
        a length, or in²/in.
    """
    if is_imperial:
        return max(0.75 * math.sqrt(f_c) / f_yt, 50 / f_yt) * b_w
    return max(0.062 * math.sqrt(f_c) / f_yt, 0.35 / f_yt) * b_w


def max_stirrup_spacing(
    V_s_req: float,
    f_c: float,
    lambda_factor: float,
    A_cv: float,
    d: float,
    *,
    is_imperial: bool = False,
) -> tuple[float, float]:
    """Maximum stirrup spacing along and across the member — ACI 318-19 Table 9.7.6.2.2.

    Once the demand on the stirrups passes the table's threshold the spacing
    limits halve, because a wider crack needs more legs crossing it.

    Args:
        V_s_req: Shear the stirrups must carry (N, or lb).
        f_c: Specified concrete compressive strength (MPa, or psi).
        lambda_factor: Lightweight concrete factor lambda.
        A_cv: Effective shear area (mm², or in²).
        d: Effective depth for shear (mm, or in).

    Returns:
        ``(s_max_l, s_max_w)`` — the limits along the member and across its
        width (mm, or in).
    """
    threshold_coeff = 4 if is_imperial else 0.083
    cap_low, cap_high = (24.0, 12.0) if is_imperial else (600.0, 300.0)

    if V_s_req <= threshold_coeff * lambda_factor * math.sqrt(f_c) * A_cv:
        return min(d / 2, cap_low), min(d, cap_low)
    return min(d / 4, cap_high), min(d / 2, cap_high)


def shear_strength_of_reinforcement(A_v: float, f_yt: float, d: float) -> float:
    """Shear carried by the stirrups V_s = A_v*f_yt*d/s — ACI 318-19 §22.5.8.5.3.

    Args:
        A_v: Stirrup area per unit length A_v/s (mm²/mm, or in²/in).
        f_yt: Yield strength used for shear (MPa, or psi).
        d: Effective depth for shear (mm, or in).

    Returns:
        V_s (N, or lb). Same expression in both unit systems.
    """
    return A_v * f_yt * d
