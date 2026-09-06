"""Public, read-only result of a punching shear check.

The counterpart of :mod:`mento.design_results` for a slab-column connection.
A check returns one of these per load combination and an envelope across them,
rather than a DataFrame: the tables are presentation and are built in
:mod:`mento.reports.punching` from these values (ADR-0001, ADR-0004)::

    governing = node.check()

    governing.DCR
    governing.v_u        # design shear stress on the critical perimeter
    governing.v_c        # the resistance it is compared against
    governing.b_0        # length of that perimeter

Every quantity is a pint ``Quantity`` in the unit system of the slab's concrete.

The two codes name these quantities differently and mento keeps its own names
here, once, rather than a union of both spellings:

===========  ============================  ==============================
This module  ACI 318-19 / CIRSOC 201-25    EN 1992-1-1:2004
===========  ============================  ==============================
``b_0``      ``b_o``, at ``d/2`` (§22.6.4)  ``u_1``, at ``2d`` (§6.4.2)
``d``        ``d`` (§22.6.4.1)              ``d_eff`` (eq. 6.32)
``v_u``      ``v_u`` (§8.4.4.2.3)           ``v_Ed`` (eq. 6.38)
``v_c``      ``φv_c`` (§22.6.5)             ``v_Rd,c`` (eq. 6.47)
===========  ============================  ==============================

The field set is deliberately the intersection that both codes report and that
an engineer reads off a check. Phase 2 may add to it as the equations land --
mento is pre-1.0 and ADR-0003 allows that -- but nothing here is expected to
change meaning.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

from pint import Quantity


class PunchingCheckNotRunError(RuntimeError):
    """Raised when a punching result is read before ``check()`` has been run."""


@dataclass(frozen=True)
class PunchingCheck:
    """The punching result of one load combination.

    ``v_c`` is the resistance the ``DCR`` was formed from, so ``v_u / DCR``
    gives it back — the same contract the beam results keep (see
    :class:`mento.design_results.ShearCheck`).
    """

    label: str
    b_0: Quantity
    d: Quantity
    v_u: Quantity
    v_c: Quantity
    DCR: float


def _governing(checks: Sequence[PunchingCheck]) -> Optional[PunchingCheck]:
    """The combination with the worst DCR; ties break to the lowest resistance.

    A punching resistance is not the connection's alone: it moves with the
    combination through the eccentricity of the load transfer, so the envelope
    has to carry the ``v_c`` that its own ``DCR`` was formed from rather than
    the largest or the smallest seen. Among combinations tied on DCR — every
    one of them when nothing is demanded — the lowest resistance is the safe
    reading.
    """
    if not checks:
        return None
    return min(checks, key=lambda check: (-check.DCR, check.v_c))


def envelope_punching(checks: Sequence[PunchingCheck]) -> PunchingCheck:
    """The governing combination, relabelled.

    A pure function of the results. Unlike flexure, where each face envelopes
    its quantities independently, a punching check is a single stress against a
    single resistance on one perimeter — so the envelope *is* one of the
    combinations, not a synthetic mixture of several.
    """
    governing = _governing(checks)
    if governing is None:
        raise ValueError("envelope_punching needs at least one checked combination.")
    return PunchingCheck(
        label="envelope",
        b_0=governing.b_0,
        d=governing.d,
        v_u=governing.v_u,
        v_c=governing.v_c,
        DCR=governing.DCR,
    )
