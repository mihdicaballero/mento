"""Two-way shear check of a slab-column connection to ACI 318-19.

The checker of ADR-0002 for punching: it converts the node's pint quantities
once, calls the float equations in
:mod:`mento.codes.aci_318_19.equations.punching`, and returns a
:class:`~mento.punching_results.PunchingCheck`. It writes nothing back onto the
node and builds no tables.

**Phase 2 is not finished.** The preconditions below are real and enforced; the
calculation between them is what the Calcpad validation still has to produce.
CIRSOC 201-25 shares this checker, as it shares the rest of the ACI formulas.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mento.punching_results import PunchingCheck

if TYPE_CHECKING:
    from mento.forces import Forces
    from mento.punching import PunchingNode


def punching_load(force: "Forces") -> "Forces":
    """The design punching load, refused if the combination carries none.

    ``V_z`` is ``V_u`` in ACI 318-19 §22.6 and ``V_Ed`` in EN 1992-1-1 §6.4 —
    the vertical force transferred at the connection. A combination that carries
    only moments would otherwise be checked at zero load and report a DCR of
    0.00, which reads as a pass.
    """
    if force.V_z.magnitude == 0:
        raise ValueError(
            f"Forces '{force.label}' carries no vertical load, so there is nothing to punch the slab. "
            "The punching load is V_z — Vu under ACI 318-19, VEd under EN 1992-1-1 — "
            "see docs/source/user_guide/local_axes.rst."
        )
    return force


def supported_geometry(node: "PunchingNode") -> "PunchingNode":
    """Refuse the geometry whose perimeter rules are not written yet.

    A capital and an opening both change the critical section rather than only
    the numbers on it, so ignoring one would not be conservative — it would be
    wrong in the unsafe direction for the opening and the safe one for the
    capital, with nothing on the result to say so.
    """
    if node.capital is not None:
        raise NotImplementedError(
            "A capital changes which critical perimeters have to be checked; that arrives in Phase 3 "
            "of docs/architecture/punching-roadmap.md."
        )
    if node.openings:
        raise NotImplementedError(
            "An opening near the column removes part of the critical perimeter; that arrives in Phase 3 "
            "of docs/architecture/punching-roadmap.md."
        )
    return node


def check_punching_ACI_318_19(node: "PunchingNode", force: "Forces") -> PunchingCheck:
    """Two-way shear at the connection, for one load combination.

    ACI 318-19 §22.6 for the resistance and §8.4.4.2 for the share of the
    unbalanced moment that reaches the critical section as shear.
    """
    supported_geometry(node)
    punching_load(force)
    raise NotImplementedError(
        "The ACI 318-19 punching check is not implemented yet: Phase 2 of "
        "docs/architecture/punching-roadmap.md. The geometry, the reinforcement and the result type are "
        "in place; the equations of mento.codes.aci_318_19.equations.punching are not."
    )
