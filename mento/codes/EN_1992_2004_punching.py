"""Punching shear check of a slab-column connection to EN 1992-1-1:2004 §6.4.

The EN counterpart of :mod:`mento.codes.ACI_318_19_punching`, and the same
contract: convert once, call the float equations of
:mod:`mento.codes.en_1992_2004.equations.punching`, return a
:class:`~mento.punching_results.PunchingCheck`.

**Phase 5 is not finished.** The preconditions below are real and enforced; the
calculation between them is what the Calcpad validation still has to produce.

The precondition EN has and ACI does not is ρ: ``v_Rd,c`` reads the flexural
reinforcement ratio, so a slab whose reinforcement was never declared cannot be
checked to this code.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mento.codes.ACI_318_19_punching import punching_load, supported_geometry
from mento.punching_results import PunchingCheck

if TYPE_CHECKING:
    from mento.forces import Forces
    from mento.punching import PunchingNode


def declared_reinforcement(node: "PunchingNode") -> "PunchingNode":
    """Refuse a slab with no ρ, which EN needs and ACI does not.

    ``v_Rd,c`` scales with ``(100·ρ_l·f_ck)^(1/3)``, so an undeclared ρ read as
    zero would drop the resistance to ``v_min`` — a quiet under-report, in a
    check whose whole output is one ratio.
    """
    if not node.slab.has_rebar:
        raise ValueError(
            "EN 1992-1-1 §6.4.4 computes v_Rd,c from the flexural reinforcement ratio, and this slab has "
            "none declared in both directions. Declare it with set_rebar_x() and set_rebar_y()."
        )
    return node


def check_punching_EN_1992_2004(node: "PunchingNode", force: "Forces") -> PunchingCheck:
    """Punching at the connection, for one load combination.

    EN 1992-1-1 §6.4: the basic control perimeter at 2d, ``v_Rd,c`` from
    eq. (6.47), and β for the eccentricity of the transfer from §6.4.3.
    """
    supported_geometry(node)
    declared_reinforcement(node)
    punching_load(force)
    raise NotImplementedError(
        "The EN 1992-1-1 punching check is not implemented yet: Phase 5 of "
        "docs/architecture/punching-roadmap.md. The geometry, the reinforcement and the result type are "
        "in place; the equations of mento.codes.en_1992_2004.equations.punching are not."
    )
