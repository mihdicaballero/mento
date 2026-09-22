"""Two-way shear check of a slab-column connection to ACI 318-19 / CIRSOC 201-25.

The checker of ADR-0002 for punching: it converts the node's pint quantities
once, calls the float equations in
:mod:`mento.codes.aci_318_19.equations.punching`, and returns a
:class:`~mento.punching_results.PunchingCheck`. It writes nothing back onto the
node and builds no tables.

**Phase 2 is not finished.** The preconditions below are real and enforced; the
calculation between them is what the Calcpad validation still has to produce.

CIRSOC 201-25 shares this checker, and every clause cited in this module carries
the same number in both codes. Two things do not match, though.

First, the cap on fyt: §22.6.3.2 sends it to ACI 318-19 §20.2.2.4 and to
CIRSOC 201-25 §20.2.1.3 — the same split as §22.5.3.3 for one-way shear, and
the same 420 MPa in practice (ACI Table 20.2.2.4(a) prints it; CIRSOC C 22.6.3.2
explains it). A reference difference, not a numeric one.

Second, and this one is numeric, Table 22.6.5.2:

* ACI 318-19 (SI) Table 22.6.5.2 prints ``(0.17 + 0.33/beta)*lambda_s*lambda*
  sqrt(f'c)`` for (b) and ``(0.17 + 0.083*alpha_s*d/b_o)*lambda_s*lambda*
  sqrt(f'c)`` for (c).
* CIRSOC 201-25 Table 22.6.5.2 prints ``0.17*(1 + 2/beta)*...`` for (b), which
  expands to ``0.17 + 0.34/beta``, and ``0.083*(2 + alpha_s*d/b_o)*...`` for
  (c), which expands to ``0.166 + 0.083*alpha_s*d/b_o``.
* ACI 318-19 (in-lb) Table 22.6.5.2 prints ``(2 + 4/beta)`` and
  ``(2 + alpha_s*d/b_o)``.

Row (a), ``0.33*lambda_s*lambda*sqrt(f'c)``, is identical, and so is the whole
of Table 22.6.6.1 — where ACI 318-19 SI itself writes rows (c) and (d) in the
same factored form CIRSOC uses in Table 22.6.5.2. The gap is below 2.5 percent,
but it is real: when Phase 2 lands, the Table 22.6.5.2 coefficients have to
reach the equations as a per-code registry datum, on the pattern of
``stirrup_spacing_caps``, and never as a ``design_code`` string comparison. No
other punching clause needs a hook, Table 22.6.6.1 least of all.

Everything else in §22.6 was read side by side and matches: §22.6.1.2
(``vn = vc``), §22.6.1.3 (``vn = vc + vs``), the 8.3 MPa ceiling on
``sqrt(f'c)`` in §22.6.3.1, the critical sections of §22.6.4, the
alpha_s of §22.6.5.3 (40 interior, 30 edge, 20 corner) and the strength
condition of §8.5.1.1(d).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mento.punching_results import PunchingCheck

if TYPE_CHECKING:
    from mento.forces import Forces
    from mento.punching import PunchingNode


def punching_load(force: "Forces") -> "Forces":
    """The design punching load, refused if the combination carries none.

    ``V_z`` is ``V_u`` in ACI 318-19 / CIRSOC 201-25 §22.6.1.4 — two-way shear
    is resisted by a section of depth d on the critical perimeter b_o of §22.6.4
    — and ``V_Ed`` in EN 1992-1-1 §6.4. It is the vertical force transferred at
    the connection; the unbalanced moment adds to it through §8.4.4.2.1 in both
    codes, where v_u combines v_uv with the stress produced by gamma_v*M_sc.

    A combination that carries only moments would otherwise be checked at zero
    load and report a DCR of 0.00, which reads as a pass.
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

    Both rules are the same in ACI 318-19 and CIRSOC 201-25:

    * §22.6.4.1(b) — b_o need not be closer than d/2 to a change in slab or
      footing thickness, such as the edge of a capital, drop panel or shear cap,
      so a capital adds critical sections instead of moving one.
    * §22.6.4.3 — for an opening closer than 4h to the periphery of the column,
      the portion of b_o enclosed by the lines projected from the column
      centroid tangent to the opening is ineffective.
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

    ACI 318-19 / CIRSOC 201-25 §22.6 for the resistance — §22.6.4 for the
    critical perimeter, §22.6.5.2 with Table 22.6.5.2 for v_c without shear
    reinforcement, §22.6.6 with Table 22.6.6.1 for v_c with it, §22.6.1.2 and
    §22.6.1.3 for ``v_n``, and §22.6.3.1 for the 8.3 MPa (100 psi) ceiling on
    ``sqrt(f'c)`` — and §8.4.4.2 for the share of the unbalanced moment that
    reaches the critical section as shear (§8.4.4.2.2 for gamma_v, §8.4.4.2.3
    for the linear variation). What has to hold is §8.5.1.1(d) in both:
    ``phi*v_n >= v_u`` at the critical sections of §8.4.4.1.

    Table 22.6.5.2 is the one place where the two codes print different
    coefficients; see the module docstring. Everything else cited here is
    word-for-word the same clause in both, which is why one checker serves them.
    """
    supported_geometry(node)
    punching_load(force)
    raise NotImplementedError(
        "The ACI 318-19 punching check is not implemented yet: Phase 2 of "
        "docs/architecture/punching-roadmap.md. The geometry, the reinforcement and the result type are "
        "in place; the equations of mento.codes.aci_318_19.equations.punching are not."
    )
