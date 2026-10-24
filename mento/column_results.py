"""Public, read-only result of the shear check of a column section.

The counterpart of :mod:`mento.design_results` for :class:`~mento.column_section.RectangularColumn`
and :class:`~mento.column_section.CircularColumn`. A check returns one of these
per load combination and direction, as values (ADR-0001)::

    for check in column.shear_check_results(forces):
        check.direction, check.DCR, check.reinforcement_required, check.A_v_req

Every quantity is a pint ``Quantity`` in the unit system of the column's
concrete. The plan this belongs to is
``docs/architecture/column-shear-roadmap.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

from mento.units import Quantity


@dataclass(frozen=True)
class ColumnShearCheck:
    """The shear result of one load combination in one direction of a column.

    ``direction`` is the shear checked: ``"z"`` reads ``V_z``, which acts
    across ``height`` with ``b_w = width``, and ``"y"`` reads ``V_y``, across
    ``width`` with ``b_w = height`` (``docs/source/user_guide/local_axes.rst``).
    A circle reads the same section both ways.

    ``V_concrete`` is the share of the concrete under the combination's axial
    load: ``ØVc`` under ACI 318-19 and CIRSOC 201-25. ``reinforcement_required``
    is whether the demand passes the share a column has to reinforce past:
    ``Vu > 0.5·ØVc`` of §10.6.2.1. Only there do ``A_v_min`` (Table 10.6.2.2),
    ``A_v_req`` and the spacing limits of Table 10.7.6.5.2 (``s_max_l``,
    ``s_max_w``) apply; elsewhere ``A_v_req`` and ``A_v_min`` are zero and the
    limits ``None``. A column is held to the ties of its detailing everywhere,
    which this check does not read.

    ``V_capacity`` is the design resistance the ``DCR`` was formed from, so
    ``V_demand / DCR`` gives it back, and ``V_max`` the most the section can
    carry however it is reinforced (ACI 318-19 §22.5.1.2, with the row of
    Table 22.5.5.1 that A_v,min selects). ``A_v`` is what the ties give in this
    direction, ``s_l`` their spacing along the column and ``s_w`` the largest
    distance between adjacent legs across it (``None`` in a circle and with no
    ties).
    """

    label: Optional[str]
    direction: Literal["y", "z"]
    N_demand: Quantity
    V_demand: Quantity
    b_w: Quantity
    d: Quantity
    V_concrete: Quantity
    reinforcement_required: bool
    A_v: Quantity
    A_v_req: Quantity
    A_v_min: Quantity
    V_capacity: Quantity
    V_max: Quantity
    DCR: float
    s_l: Optional[Quantity] = None
    s_w: Optional[Quantity] = None
    s_max_l: Optional[Quantity] = None
    s_max_w: Optional[Quantity] = None
    spacing_halved: Optional[bool] = None
