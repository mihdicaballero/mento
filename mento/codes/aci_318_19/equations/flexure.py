"""Flexure equations — ACI 318-19, Chapters 9, 21 and 22.

Pure functions of floats; see the package docstring for the unit convention.
Areas are mm² (or in²) and moments N·mm (or lb·in), consistent with the stress
and length units of the system in use.
"""

from __future__ import annotations

import math

__all__ = [
    "max_reinforcement_ratio",
    "min_reinforcement_ratio",
    "neutral_axis_at_ductility_limit",
    "compression_steel_net_stress",
    "flexural_resistance_factor",
    "singly_reinforced_discriminant",
    "tension_steel_for_moment",
    "neutral_axis_depth",
    "nominal_moment_singly_reinforced",
    "nominal_moment_doubly_reinforced",
]


def max_reinforcement_ratio(f_c: float, f_y: float, beta_1: float, epsilon_c: float, epsilon_y: float) -> float:
    """Maximum tension reinforcement ratio rho_max — ACI 318-19 §21.2.2.

    The ductility limit: the largest ratio that still leaves the section
    tension-controlled, from strain compatibility at eps_t = eps_y + 2*eps_c.
    Form taken from the CRSI Design Guide, Beam Theory p. 6-3.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        f_y: Specified steel yield strength, same unit as ``f_c``.
        beta_1: Stress block factor, ACI Table 22.2.2.4.3.
        epsilon_c: Concrete crushing strain (0.003 in ACI).
        epsilon_y: Steel yield strain.

    Returns:
        rho_max, dimensionless.
    """
    return 0.85 * beta_1 * f_c / f_y * (epsilon_c / (epsilon_y + epsilon_c * 2))


def min_reinforcement_ratio(f_c: float, f_y: float, *, is_imperial: bool = False) -> float:
    """Minimum flexural reinforcement ratio rho_min — ACI 318-19 §9.6.1.2.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        f_y: Specified steel yield strength, same unit as ``f_c``.

    Returns:
        rho_min, dimensionless. The clause is the larger of a sqrt(f_c) term
        and a floor that governs at low concrete strengths.
    """
    if is_imperial:
        return max(3 * math.sqrt(f_c) / f_y, 200 / f_y)
    return max(0.25 * math.sqrt(f_c) / f_y, 1.4 / f_y)


def neutral_axis_at_ductility_limit(d: float, epsilon_y: float) -> float:
    """Neutral axis depth at the tension-controlled boundary — ACI 318-19 §21.2.2.

    At the limit eps_t = eps_y + 0.006, so strain compatibility
    (c/d = eps_cu / (eps_cu + eps_t)) gives c_t = 0.003*d / (eps_y + 0.006).

    Args:
        d: Effective depth of the tension reinforcement (mm, or in).
        epsilon_y: Steel yield strain.

    Returns:
        c_t (mm, or in). A section with c < c_t is tension-controlled.
    """
    return 0.003 * d / (epsilon_y + 0.006)


def compression_steel_net_stress(d_prime: float, c_t: float, E_s: float, f_y: float, f_c: float) -> float:
    """Compression steel stress at the ductility limit, net of displaced concrete.

    ACI 318-19 §22.2 (strain compatibility) with the §22.2.2.4.1 stress block.

        eps_s' = (c_t - d')/c_t * 0.003
        f_s'   = min(eps_s'*E_s, f_y)
        net    = f_s' - 0.85*f_c

    The subtraction accounts for the concrete the bar displaces: that volume is
    already carrying 0.85*f_c inside the rectangular block, so counting the full
    bar stress would count it twice.

    Args:
        d_prime: Depth to the compression reinforcement (mm, or in).
        c_t: Neutral axis at the ductility limit, from
            :func:`neutral_axis_at_ductility_limit` (mm, or in).
        E_s: Steel modulus of elasticity (MPa, or psi).
        f_y: Steel yield strength (MPa, or psi).
        f_c: Concrete compressive strength (MPa, or psi).

    Returns:
        The net compression steel stress (MPa, or psi).
    """
    f_s_prime = min(0.003 * E_s * (1 - d_prime / c_t), f_y)
    return f_s_prime - 0.85 * f_c


def flexural_resistance_factor(M_u: float, phi: float, b: float, d: float) -> float:
    """Required flexural resistance R_n = M_u/(phi*b*d²) — ACI 318-19 §22.3.

    Args:
        M_u: Factored moment, positive (N·mm, or lb·in).
        phi: Strength reduction factor for flexure.
        b: Section width (mm, or in).
        d: Effective depth (mm, or in).

    Returns:
        R_n as a stress (MPa, or psi).
    """
    return M_u / (phi * b * d**2)


def singly_reinforced_discriminant(R_n: float, f_c: float) -> float:
    """Term under the root of the singly reinforced solution — ACI 318-19 §22.2.

    Args:
        R_n: Required flexural resistance from :func:`flexural_resistance_factor`.
        f_c: Concrete compressive strength, same unit as ``R_n``.

    Returns:
        ``1 - 2*R_n/(0.85*f_c)``. Negative means the demand exceeds what the
        section can carry with tension steel alone, and the caller must fall
        back to a doubly reinforced solution or report the section as
        insufficient — that decision is orchestration, not this equation.
    """
    return 1 - 2 * R_n / (0.85 * f_c)


def tension_steel_for_moment(R_n: float, f_c: float, f_y: float, b: float, d: float) -> float:
    """Tension steel required for a moment — ACI 318-19 §22.2, singly reinforced.

    Inverts the rectangular stress block for A_s. Only valid when
    :func:`singly_reinforced_discriminant` is non-negative.

    Args:
        R_n: Required flexural resistance (MPa, or psi).
        f_c: Concrete compressive strength (MPa, or psi).
        f_y: Steel yield strength (MPa, or psi).
        b: Section width (mm, or in).
        d: Effective depth (mm, or in).

    Returns:
        A_s (mm², or in²).
    """
    return 0.85 * f_c * b * d / f_y * (1 - math.sqrt(singly_reinforced_discriminant(R_n, f_c)))


def neutral_axis_depth(A_s: float, f_y: float, f_c: float, b: float, beta_1: float) -> float:
    """Neutral axis depth from equilibrium — ACI 318-19 §22.2.2.4.1.

    From 0.85*f_c*(beta_1*c)*b = A_s*f_y.

    Args:
        A_s: Tension reinforcement area (mm², or in²).
        f_y: Steel yield strength (MPa, or psi).
        f_c: Concrete compressive strength (MPa, or psi).
        b: Section width (mm, or in).
        beta_1: Stress block factor.

    Returns:
        c (mm, or in).
    """
    return A_s * f_y / (0.85 * f_c * b * beta_1)


def nominal_moment_singly_reinforced(A_s: float, f_y: float, f_c: float, b: float, d: float) -> float:
    """Nominal moment of a singly reinforced section — ACI 318-19 §22.3.

    Force equilibrium 0.85*f_c*a*b = A_s*f_y fixes the stress block depth
    a = A_s*f_y/(0.85*f_c*b); the lever arm is then d - a/2. Valid only while
    A_s does not exceed the maximum of §21.2.2.

    Args:
        A_s: Tension reinforcement area (mm², or in²).
        f_y: Steel yield strength (MPa, or psi).
        f_c: Concrete compressive strength (MPa, or psi).
        b: Section width (mm, or in).
        d: Effective depth (mm, or in).

    Returns:
        M_n (N·mm, or lb·in).
    """
    a = A_s * f_y / (0.85 * f_c * b)
    return A_s * f_y * (d - a / 2)


def nominal_moment_doubly_reinforced(
    A_s: float,
    A_s_prime: float,
    f_y: float,
    f_c: float,
    b: float,
    d: float,
    d_prime: float,
    beta_1: float,
    epsilon_c: float,
    epsilon_y: float,
    E_s: float,
) -> float:
    """Nominal moment with compression reinforcement — ACI 318-19 §22.3.

    Two cases, decided by whether the compression steel yields.

    Equilibrium carries the displaced-concrete correction throughout: the
    compression bar sits inside the 0.85*f_c*a*b block, so its effective
    contribution is (f_y - 0.85*f_c), not f_y.

    First the neutral axis is solved assuming the compression steel yields::

        c = (A_s*f_y - A_s'*(f_y - 0.85*f_c)) / (0.85*f_c*b*beta_1)

    If the resulting strain reaches eps_y the assumption holds. Otherwise c is
    re-solved from the quadratic that keeps the compression steel elastic,
    A*c² + B*c + C = 0 with A = 0.85*f_c*b*beta_1,
    B = A_s'*(eps_c*E_s - 0.85*f_c) - A_s*f_y and C = -d'*A_s'*eps_c*E_s.

    Args:
        A_s: Tension reinforcement area (mm², or in²).
        A_s_prime: Compression reinforcement area (mm², or in²).
        f_y: Steel yield strength (MPa, or psi).
        f_c: Concrete compressive strength (MPa, or psi).
        b: Section width (mm, or in).
        d: Effective depth of the tension steel (mm, or in).
        d_prime: Depth to the compression steel (mm, or in).
        beta_1: Stress block factor.
        epsilon_c: Concrete crushing strain.
        epsilon_y: Steel yield strain.
        E_s: Steel modulus of elasticity (MPa, or psi).

    Returns:
        M_n (N·mm, or lb·in).
    """
    # Step 1: assume the compression steel yields and solve for the neutral axis.
    c_assumed = (A_s * f_y - A_s_prime * (f_y - 0.85 * f_c)) / (0.85 * f_c * b * beta_1)

    # A zero neutral axis happens when top and bottom steel are equal; the
    # compression steel is then certainly not yielding.
    epsilon_s = (c_assumed - d_prime) / c_assumed * epsilon_c if c_assumed > 0 else 0

    # Step 2: if the assumed strain reaches yield, the assumption stands.
    if epsilon_s >= epsilon_y:
        a_assumed = c_assumed * beta_1
        return 0.85 * f_c * a_assumed * b * (d - a_assumed / 2) + A_s_prime * (f_y - 0.85 * f_c) * (d - d_prime)

    # Otherwise re-solve with the compression steel still elastic.
    A = 0.85 * f_c * b * beta_1
    B = A_s_prime * (epsilon_c * E_s - 0.85 * f_c) - A_s * f_y
    C = -d_prime * A_s_prime * epsilon_c * E_s

    c = (-B + math.sqrt(B**2 - 4 * A * C)) / (2 * A)
    a = c * beta_1

    f_s_prime_net = (c - d_prime) / c * epsilon_c * E_s - 0.85 * f_c

    # Part of the tension steel is balanced by the net compression steel force.
    A_s_2 = A_s_prime * f_s_prime_net / f_y
    A_s_1 = A_s - A_s_2

    return A_s_1 * f_y * (d - a / 2) + A_s_prime * f_s_prime_net * (d - d_prime)
