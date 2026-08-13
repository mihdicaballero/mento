"""Public, read-only view of the reinforcement produced by a check or a design.

The element classes keep their results in private attributes (``_A_s_bot``,
``_stirrup_s_l`` and so on). Those are implementation details and their names
and units can change. The dataclasses in this module are the supported way to
read a result from code::

    node.design()

    beam.flexure_design.bottom.A_s        # provided bottom steel area
    beam.flexure_design.bottom.A_s_req    # what the design required
    beam.shear_design.s_l                 # stirrup spacing

Every quantity is returned as a pint ``Quantity`` in the unit system of the
section's concrete, so a value can be converted with ``.to()`` as usual.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, Optional, Tuple

from pint import Quantity

if TYPE_CHECKING:
    from mento.beam import RectangularBeam


class DesignNotRunError(RuntimeError):
    """Raised when results are read before a check or design has been run."""


@dataclass(frozen=True)
class RebarLayer:
    """One layer of longitudinal bars: ``n`` bars of diameter ``d_b``."""

    n: int
    d_b: Quantity

    @property
    def A_s(self) -> Quantity:
        """Steel area of this layer."""
        return self.n * (self.d_b**2) * math.pi / 4

    def __str__(self) -> str:
        return f"{self.n}Ø{self.d_b:.4g~P}"


@dataclass(frozen=True)
class FlexureFaceDesign:
    """Longitudinal reinforcement on one face of the section.

    ``layers`` only contains the layers that actually carry bars, so an
    empty tuple means no reinforcement was placed on this face.

    ``A_s`` is what the section carries. ``A_s_req``, ``A_s_min``, ``A_s_max``
    and ``DCR`` are the envelope over every load combination that was checked,
    so ``DCR`` is the one of the combination that governs this face.
    """

    layers: Tuple[RebarLayer, ...]
    A_s: Quantity
    A_s_req: Quantity
    A_s_min: Quantity
    A_s_max: Quantity
    DCR: float

    @property
    def n_bars(self) -> int:
        """Total number of bars across every layer of this face."""
        return sum(layer.n for layer in self.layers)

    def __str__(self) -> str:
        if not self.layers:
            return "no reinforcement"
        return " + ".join(str(layer) for layer in self.layers)


@dataclass(frozen=True)
class FlexureDesign:
    """Longitudinal reinforcement of the whole section."""

    bottom: FlexureFaceDesign
    top: FlexureFaceDesign

    @property
    def DCR(self) -> float:
        """Governing demand-to-capacity ratio of the two faces."""
        return max(self.bottom.DCR, self.top.DCR)

    def __str__(self) -> str:
        return f"bottom: {self.bottom} / top: {self.top}"


@dataclass(frozen=True)
class ShearDesign:
    """Transverse reinforcement of the section.

    ``n_stirrups`` counts the stirrups; ``n_legs`` counts the legs crossing
    the shear plane, which is what enters the ``A_v`` calculation.

    ``A_v_req``, ``A_v_min`` and ``DCR`` are the envelope over every load
    combination that was checked, so ``DCR`` is the governing one.
    """

    n_stirrups: int
    d_b: Quantity
    s_l: Quantity
    A_v: Quantity
    A_v_req: Quantity
    A_v_min: Quantity
    DCR: float

    @property
    def n_legs(self) -> int:
        """Number of stirrup legs crossing the shear plane."""
        return self.n_stirrups * 2

    def __str__(self) -> str:
        if self.n_stirrups == 0:
            return "no stirrups"
        return f"{self.n_stirrups}eØ{self.d_b:.4g~P}/{self.s_l:.4g~P}"


def _layers(beam: RectangularBeam, face: str) -> Tuple[RebarLayer, ...]:
    """Collect the non-empty rebar layers of one face ("b" for bottom, "t" for top)."""
    layers = []
    for index in (1, 2, 3, 4):
        n = getattr(beam, f"_n{index}_{face}", 0)
        d_b: Optional[Quantity] = getattr(beam, f"_d_b{index}_{face}", None)
        if n and d_b is not None and d_b.magnitude > 0:
            layers.append(RebarLayer(n=int(n), d_b=d_b))
    return tuple(layers)


def _face(beam: RectangularBeam, face: str) -> FlexureFaceDesign:
    suffix = "bot" if face == "b" else "top"
    # _A_s_bot is set when the section is built, so it always carries the area
    # unit of this beam and can seed the values a check may not have produced.
    zero: Quantity = 0 * beam._A_s_bot.units
    # The envelope over every combination that was checked. The private attributes
    # only describe the combination that ran last, which is not necessarily the one
    # that governs this face.
    envelope: Dict[str, Any] = getattr(beam, "_flexure_envelope", {}).get(suffix, {})

    def area(name: str) -> Quantity:
        """Envelope value of a per-combination area, or the last one if there is none."""
        value: Optional[Quantity] = envelope.get(name, getattr(beam, f"_{name}_{suffix}", None))
        return zero if value is None else value

    # A_s is the reinforcement the section carries, not a per-combination result.
    A_s: Optional[Quantity] = getattr(beam, f"_A_s_{suffix}", None)

    return FlexureFaceDesign(
        layers=_layers(beam, face),
        A_s=zero if A_s is None else A_s,
        A_s_req=area("A_s_req"),
        A_s_min=area("A_s_min"),
        A_s_max=area("A_s_max"),
        DCR=float(envelope.get("DCR", getattr(beam, f"_DCRb_{suffix}", 0.0))),
    )


def build_flexure_design(beam: RectangularBeam) -> FlexureDesign:
    """Build the public flexure result of ``beam``.

    Raises:
        DesignNotRunError: if no flexure check or design has been run yet.
    """
    if not getattr(beam, "_flexure_checked", False):
        raise DesignNotRunError(
            "No flexure results yet. Run node.design() or node.check_flexure() before reading flexure_design."
        )
    return FlexureDesign(bottom=_face(beam, "b"), top=_face(beam, "t"))


def build_shear_design(beam: RectangularBeam) -> ShearDesign:
    """Build the public shear result of ``beam``.

    Raises:
        DesignNotRunError: if no shear check or design has been run yet.
    """
    if not getattr(beam, "_shear_checked", False):
        raise DesignNotRunError(
            "No shear results yet. Run node.design() or node.check_shear() before reading shear_design."
        )
    A_v: Quantity = beam._A_v
    zero: Quantity = 0 * A_v.units
    # As in _face: the private attributes describe the combination that ran last, the
    # envelope describes every combination that was checked.
    envelope: Dict[str, Any] = getattr(beam, "_shear_envelope", {})

    def area(name: str) -> Quantity:
        """Envelope value of a per-combination area, or the last one if there is none."""
        value: Optional[Quantity] = envelope.get(name, getattr(beam, f"_{name}", None))
        return zero if value is None else value

    return ShearDesign(
        n_stirrups=int(beam._stirrup_n),
        d_b=beam._stirrup_d_b,
        s_l=beam._stirrup_s_l,
        A_v=A_v,
        A_v_req=area("A_v_req"),
        A_v_min=area("A_v_min"),
        DCR=float(envelope.get("DCR", getattr(beam, "_DCRv", 0.0))),
    )
