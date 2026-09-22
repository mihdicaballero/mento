"""Flexure equations — ACI 318-19, Chapters 9, 21 and 22.

CIRSOC 201-25 reprints these chapters under the same numbering, so every
function below carries both citations. One thing the two codes print
differently reaches this module, and it arrives as an argument rather than
written in: the cap on f_y inside §9.6.1.2 (ACI 550 MPa / 80,000 psi against
CIRSOC 500 MPa), which the caller reads off the code's registry entry.

Pure functions of floats; see the package docstring for the unit convention.
Areas are mm² (or in²) and moments N·mm (or lb·in), consistent with the stress
and length units of the system in use.
"""

from __future__ import annotations

import math

__all__ = [
    "max_reinforcement_ratio",
    "min_reinforcement_ratio",
    "shrinkage_and_temperature_ratio",
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
    """Maximum tension reinforcement ratio rho_max — ACI 318-19 §9.3.3.1 / CIRSOC 201-25 §9.3.3.1.

    The ductility limit: the largest ratio that still leaves the section
    tension-controlled. A nonprestressed beam has to be tension-controlled per
    Table 21.2.2 of both codes (a slab too, by ACI 318-19 §7.3.3.1 / CIRSOC
    201-25 §7.3.3), and that table puts the limit at eps_t = eps_ty + 0.003.
    With eps_cu = 0.003 (§22.2.2.1 in both) and linear strains (§22.2.1.2 in
    both), c/d = eps_cu/(eps_cu + eps_t), so the denominator below is
    eps_y + 2*eps_c — which is the sum, not eps_t itself. Form taken from the
    CRSI Design Guide, Beam Theory p. 6-3.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        f_y: Specified steel yield strength, same unit as ``f_c``.
        beta_1: Stress block factor, ACI 318-19 Table 22.2.2.4.3 / CIRSOC
            201-25 Tabla 22.2.2.4.3.
        epsilon_c: Concrete crushing strain (0.003 in ACI).
        epsilon_y: Steel yield strain.

    Returns:
        rho_max, dimensionless.
    """
    return 0.85 * beta_1 * f_c / f_y * (epsilon_c / (epsilon_y + epsilon_c * 2))


def min_reinforcement_ratio(
    f_c: float, f_y: float, f_y_cap: float | None = None, *, is_imperial: bool = False
) -> float:
    """Minimum flexural reinforcement ratio rho_min — ACI 318-19 §9.6.1.2 / CIRSOC 201-25 §9.6.1.2.

    Expressions (a) and (b) are the same in both codes: 0.25*sqrt(f_c)/f_y and
    1.4/f_y in SI, 3*sqrt(f_c)/f_y and 200/f_y in psi.

    Both clauses also cap the f_y that may be put into them, and that cap is
    where they differ: "The value of fy shall be limited to a maximum of
    550 MPa" (ACI 318-19 §9.6.1.2; "80,000 psi" in the in-lb edition) against
    "El valor de fy debe limitarse a un máximo de 500 MPa" (CIRSOC 201-25
    §9.6.1.2). What it caps is the value entering the formula, not the steel
    the section is designed with, and since f_y is in the denominator of both
    expressions the cap only ever raises rho_min — ignoring it is
    unconservative. Which of the two numbers applies is a datum of the design
    code, so it arrives as ``f_y_cap``, the way the absolute spacing caps of
    Table 9.7.6.2.2 reach :func:`shear.max_stirrup_spacing`. Neither cap is
    reached by the ADN 420 of ordinary practice.

    Args:
        f_c: Specified concrete compressive strength (MPa, or psi).
        f_y: Specified steel yield strength, same unit as ``f_c``.
        f_y_cap: Largest f_y this code lets into the formula, same unit as
            ``f_y``: 550 MPa (80,000 psi) under ACI 318-19, 500 MPa under
            CIRSOC 201-25. Omitted, it defaults to the ACI cap in the selected
            unit system, preserving the two-argument equation interface.

    Returns:
        rho_min, dimensionless. The clause is the larger of a sqrt(f_c) term
        and a floor that governs at low concrete strengths.
    """
    cap = (80000.0 if is_imperial else 550.0) if f_y_cap is None else f_y_cap
    f_y_eff = min(f_y, cap)
    if is_imperial:
        return max(3 * math.sqrt(f_c) / f_y_eff, 200 / f_y_eff)
    return max(0.25 * math.sqrt(f_c) / f_y_eff, 1.4 / f_y_eff)


def shrinkage_and_temperature_ratio(f_y: float, *, is_imperial: bool = False) -> float:
    """Shrinkage and temperature reinforcement ratio — ACI 318-19 §24.4.3.2 / CIRSOC 201-25 §24.4.3.2.

    The minimum that governs a member on the ground, written on the GROSS
    section ``b*h`` rather than on the effective depth the flexural minimum of
    §9.6.1.2 uses. The path to it is the same in both codes: ACI 318-19
    §13.3.2.1 sends a one-way shallow foundation to Chapters 7 and 9, and
    §7.6.1.1 asks there for A_s,min = 0.0018*A_g, the same ratio §24.4.3.2
    gives (two-way: §13.3.3.1 → §8.6.1.1). CIRSOC 201-25 §13.3.2.1 → §7.6.1
    (an unnumbered paragraph, "As,min de 0,0018Ag") and §13.3.3.1 → §8.6.1.1.

    Both codes now state one flat ratio for every f_y: 0.0018. The scaling by
    f_y and the 0.0014 floor below are the Table 24.4.3.2 of earlier editions,
    which neither code carries any more — ACI 318-19 R24.4.3.2 records that the
    reduction for f_y over 420 MPa was withdrawn because increased yield
    strength gives no benefit for crack control. Correcting the returned value
    is a behavioural change and is not made here.

    Args:
        f_y: Specified steel yield strength (MPa, or psi).
        is_imperial: Selects the reference yield strength the withdrawn table
            was anchored at: 60000 psi in US customary, 420 MPa in SI.

    Returns:
        A_s,min/(b*h), dimensionless. 0.0018 at the reference yield strength,
        scaled inversely for a stronger steel, and never below the 0.0014
        floor.
    """
    f_y_reference = 60000.0 if is_imperial else 420.0
    return max(0.0018 * f_y_reference / f_y, 0.0014)


def neutral_axis_at_ductility_limit(d: float, epsilon_y: float) -> float:
    """Neutral axis depth at the tension-controlled boundary — ACI 318-19 Table 21.2.2 / CIRSOC 201-25 Tabla 21.2.2.

    The limit is eps_t = eps_ty + 0.003 in both codes, so with eps_cu = 0.003
    (§22.2.2.1 in both) strain compatibility (§22.2.1.2 in both,
    c/d = eps_cu / (eps_cu + eps_t)) gives
    c_t = 0.003*d / (0.003 + eps_y + 0.003) = 0.003*d / (eps_y + 0.006).

    Args:
        d: Effective depth of the tension reinforcement (mm, or in).
        epsilon_y: Steel yield strain.

    Returns:
        c_t (mm, or in). A section with c < c_t is tension-controlled.
    """
    return 0.003 * d / (epsilon_y + 0.006)


def compression_steel_net_stress(d_prime: float, c_t: float, E_s: float, f_y: float, f_c: float) -> float:
    """Compression steel stress at the ductility limit, net of displaced concrete.

    ACI 318-19 §22.2.1.2 / CIRSOC 201-25 §22.2.1.2 (strains proportional to the
    distance from the neutral axis) with the equivalent rectangular stress block
    of ACI 318-19 §22.2.2.4.1 / CIRSOC 201-25 §22.2.2.4.1, and the bar stress
    capped at f_y by ACI 318-19 §20.2.2.1 / CIRSOC 201-25 §20.2.2.1. The same
    three clauses in both codes.

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
    """Required flexural resistance R_n = M_u/(phi*b*d²) — ACI 318-19 §22.3.1.1 / CIRSOC 201-25 §22.3.1.1.

    Not a published formula: the textbook rearrangement of phi*M_n >= M_u under
    the §22.2 assumptions, which both codes share. Cited so the derivation can
    be traced, not because either code prints it.

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
    """Term under the root of the singly reinforced solution — ACI 318-19 §22.2 / CIRSOC 201-25 §22.2.

    Part of the same textbook inversion as :func:`tension_steel_for_moment`,
    resting on the §22.2 assumptions, which both codes print alike.

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
    """Tension steel required for a moment — ACI 318-19 §22.2 / CIRSOC 201-25 §22.2, singly reinforced.

    Inverts the rectangular stress block of ACI 318-19 §22.2.2.4.1 / CIRSOC
    201-25 §22.2.2.4.1 for A_s — a derivation, not a published formula, and the
    same one under either code. Only valid when
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
    """Neutral axis depth from equilibrium — ACI 318-19 §22.2.2.4.1 / CIRSOC 201-25 §22.2.2.4.1.

    From 0.85*f_c*(beta_1*c)*b = A_s*f_y, with a = beta_1*c. Both codes print
    the same stress block: 0.85*f_c uniform over a depth a.

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
    """Nominal moment of a singly reinforced section — ACI 318-19 §22.3.1.1 / CIRSOC 201-25 §22.3.1.1.

    Force equilibrium 0.85*f_c*a*b = A_s*f_y over the rectangular stress block
    of §22.2.2.4.1 in both codes fixes the depth a = A_s*f_y/(0.85*f_c*b); the
    lever arm is then d - a/2. Valid only while A_s does not exceed the maximum
    that keeps the section tension-controlled (Table 21.2.2 in both codes), for
    which see :func:`max_reinforcement_ratio`.

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
    """Nominal moment with compression reinforcement — ACI 318-19 §22.3.1.1 / CIRSOC 201-25 §22.3.1.1.

    Two cases, decided by whether the compression steel yields. The assumptions
    behind both are §22.2 of either code: linear strains (§22.2.1.2), the
    rectangular stress block (§22.2.2.4.1) and f_s = E_s*eps_s capped at f_y
    (§20.2.2.1). Nothing here differs between the two codes.

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
