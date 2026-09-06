"""Two-way (punching) shear equations of ACI 318-19.

**Not implemented yet.** This module is the shape of Phase 2, not its content:
every function below is the signature its clause asks for, with the units each
argument is expected in, and a body that raises. The formulas are being
recreated from a validated Calcpad sheet and land one at a time, each with the
worked example that checks it.

Per ADR-0005 these are free functions over plain floats. Two unit systems, one
per set of published coefficients:

- metric — N, mm, MPa (= N/mm²), N·mm
- imperial — lb, in, psi, lb·in

ACI publishes separately rounded coefficients for the two (``0.33√f'c`` in MPa
against ``4√f'c`` in psi), so a function whose coefficients differ takes an
``imperial: bool`` and works in either. The unit *system* is a parameter; the
units within it are fixed and stated per argument.

Geometry — where the critical section runs, how long it is, and where its
centroid sits once an opening or a free edge has clipped it — is not here. It
belongs to the element, which hands this module a length and a section modulus
already computed; see :mod:`mento.punching`.
"""

from __future__ import annotations

_NOT_YET = "ACI 318-19 punching equations are not implemented yet (Phase 2 of docs/architecture/punching-roadmap.md)."


def critical_section_offset(d: float) -> float:
    """Distance from the column face to the critical section, ACI 318-19 §22.6.4.1.

    Args:
        d: effective depth [mm | in].

    Returns:
        The offset [mm | in]. Half the effective depth, in either unit system.
    """
    raise NotImplementedError(_NOT_YET)


def size_effect_factor(d: float, imperial: bool = False) -> float:
    """λ_s, the size effect factor, ACI 318-19 §22.5.5.1.3.

    Args:
        d: effective depth [mm | in].
        imperial: work in (lb, in, psi) instead of (N, mm, MPa).

    Returns:
        λ_s, dimensionless and not greater than 1.0.
    """
    raise NotImplementedError(_NOT_YET)


def location_factor(position: str) -> float:
    """α_s, the factor for the column's position in the slab, ACI 318-19 §22.6.5.3.

    Args:
        position: ``"interior"``, ``"edge"`` or ``"corner"``.

    Returns:
        α_s, dimensionless.
    """
    raise NotImplementedError(_NOT_YET)


def concrete_shear_stress(
    f_c: float,
    beta: float,
    alpha_s: float,
    d: float,
    b_0: float,
    lambda_s: float,
    lambda_factor: float = 1.0,
    imperial: bool = False,
) -> float:
    """v_c, the two-way concrete shear stress, ACI 318-19 §22.6.5.2 (table).

    The least of the three expressions in the table governs.

    Args:
        f_c: concrete compressive strength [MPa | psi].
        beta: ratio of the long to the short side of the column, dimensionless.
        alpha_s: the location factor of :func:`location_factor`.
        d: effective depth [mm | in].
        b_0: length of the critical perimeter [mm | in].
        lambda_s: the size effect factor of :func:`size_effect_factor`.
        lambda_factor: λ, the lightweight concrete factor of §19.2.4.
        imperial: work in (lb, in, psi) instead of (N, mm, MPa).

    Returns:
        v_c [MPa | psi].
    """
    raise NotImplementedError(_NOT_YET)


def moment_fraction_by_flexure(b_1: float, b_2: float) -> float:
    """γ_f, the fraction of the unbalanced moment carried in flexure, ACI 318-19 §8.4.2.3.4.

    Args:
        b_1: width of the critical section measured in the direction of the span
            the moment acts in [mm | in].
        b_2: width of the critical section measured perpendicular to it [mm | in].

    Returns:
        γ_f, dimensionless.
    """
    raise NotImplementedError(_NOT_YET)


def moment_fraction_by_shear(gamma_f: float) -> float:
    """γ_v, the fraction of the unbalanced moment carried in shear, ACI 318-19 §8.4.4.2.2.

    Args:
        gamma_f: the flexural fraction of :func:`moment_fraction_by_flexure`.

    Returns:
        γ_v = 1 − γ_f, dimensionless.
    """
    raise NotImplementedError(_NOT_YET)


def shear_stress(
    V_u: float,
    b_0: float,
    d: float,
    gamma_v_x: float,
    M_sc_x: float,
    c_x: float,
    J_c_x: float,
    gamma_v_y: float,
    M_sc_y: float,
    c_y: float,
    J_c_y: float,
) -> float:
    """v_u, the factored shear stress at the critical point, ACI 318-19 §8.4.4.2.3.

    The direct stress plus the contribution of each unbalanced moment, evaluated
    at the point of the critical section farthest from its centroid.

    Args:
        V_u: factored shear transferred to the slab [N | lb] — the ``V_z`` of
            the load combination, see ``docs/source/user_guide/local_axes.rst``.
        b_0: length of the critical perimeter [mm | in].
        d: effective depth [mm | in].
        gamma_v_x: γ_v for the moment about x, from :func:`moment_fraction_by_shear`.
        M_sc_x: unbalanced moment transferred about x [N·mm | lb·in].
        c_x: distance from the centroid of the critical section to the point
            being checked, measured along x [mm | in].
        J_c_x: the section property analogous to a polar moment of inertia,
            about x [mm⁴ | in⁴].
        gamma_v_y: γ_v for the moment about y.
        M_sc_y: unbalanced moment transferred about y [N·mm | lb·in].
        c_y: distance to the point being checked, along y [mm | in].
        J_c_y: the section property about y [mm⁴ | in⁴].

    Returns:
        v_u [MPa | psi].
    """
    raise NotImplementedError(_NOT_YET)
