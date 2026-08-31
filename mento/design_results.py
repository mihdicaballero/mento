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
from typing import TYPE_CHECKING, Any, Optional, Sequence, Tuple

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
class FlexureFaceCheck:
    """What a single load combination demanded of one face.

    A field is ``None`` when the design code did not set it for this
    combination; enveloping skips those rather than treating them as zero.
    """

    A_s_req: Optional[Quantity]
    A_s_min: Optional[Quantity]
    A_s_max: Optional[Quantity]
    DCR: float


@dataclass(frozen=True)
class FlexureCheck:
    """The flexure result of one load combination, on both faces."""

    label: str
    bottom: FlexureFaceCheck
    top: FlexureFaceCheck


@dataclass(frozen=True)
class ShearCheck:
    """The shear result of one load combination."""

    label: str
    A_v_req: Optional[Quantity]
    A_v_min: Optional[Quantity]
    DCR: float


def _worst(values: Sequence[Optional[Quantity]]) -> Optional[Quantity]:
    """Largest of the values that are present, or None if none of them are."""
    present = [value for value in values if value is not None]
    return max(present) if present else None


def envelope_flexure_face(checks: Sequence[FlexureCheck], face: str) -> FlexureFaceCheck:
    """Worst demand on one face across every combination checked.

    A pure function of the results: the governing combination differs per
    quantity and per face, so each one is enveloped independently. ``face`` is
    ``"bottom"`` or ``"top"``.
    """
    faces = [getattr(check, face) for check in checks]
    return FlexureFaceCheck(
        A_s_req=_worst([f.A_s_req for f in faces]),
        A_s_min=_worst([f.A_s_min for f in faces]),
        A_s_max=_worst([f.A_s_max for f in faces]),
        DCR=max([f.DCR for f in faces], default=0.0),
    )


def envelope_shear(checks: Sequence[ShearCheck]) -> ShearCheck:
    """Worst shear demand across every combination checked."""
    return ShearCheck(
        label="envelope",
        A_v_req=_worst([c.A_v_req for c in checks]),
        A_v_min=_worst([c.A_v_min for c in checks]),
        DCR=max([c.DCR for c in checks], default=0.0),
    )


def capture_flexure_check(beam: RectangularBeam, label: str, state: Any = None) -> FlexureCheck:
    """The flexure result of the combination just run.

    Reads ``state`` when the design code returned one, so nothing had to be
    written to the beam; falls back to its attributes for codes still on the
    older path.
    """
    if state is not None:
        return FlexureCheck(
            label=label,
            bottom=FlexureFaceCheck(
                A_s_req=state.A_s_req_bot,
                A_s_min=state.A_s_min_bot,
                A_s_max=state.A_s_max_bot,
                DCR=float(state.DCR_bot),
            ),
            top=FlexureFaceCheck(
                A_s_req=state.A_s_req_top,
                A_s_min=state.A_s_min_top,
                A_s_max=state.A_s_max_top,
                DCR=float(state.DCR_top),
            ),
        )

    def face(suffix: str) -> FlexureFaceCheck:
        return FlexureFaceCheck(
            A_s_req=getattr(beam, f"_A_s_req_{suffix}", None),
            A_s_min=getattr(beam, f"_A_s_min_{suffix}", None),
            A_s_max=getattr(beam, f"_A_s_max_{suffix}", None),
            DCR=float(getattr(beam, f"_DCRb_{suffix}", 0.0)),
        )

    return FlexureCheck(label=label, bottom=face("bot"), top=face("top"))


def capture_shear_check(beam: RectangularBeam, label: str, state: Any = None) -> ShearCheck:
    """The shear result of the combination just run.

    Reads ``state`` when the design code returned one — then nothing had to be
    written to the beam for this to work. Falls back to the beam's attributes
    for the codes still on the older path.
    """
    if state is not None:
        A_v_req, A_v_min = state.shear_reinforcement_quantities(beam.concrete.is_imperial)
        return ShearCheck(label=label, A_v_req=A_v_req, A_v_min=A_v_min, DCR=float(state.DCR))
    return ShearCheck(
        label=label,
        A_v_req=getattr(beam, "_A_v_req", None),
        A_v_min=getattr(beam, "_A_v_min", None),
        DCR=float(getattr(beam, "_DCRv", 0.0)),
    )


@dataclass(frozen=True)
class FaceReinforcement:
    """The longitudinal bars a section carries on one face.

    This is *configuration*, not a result: it says what is detailed on the
    section right now, whether that came from a design, from
    ``set_longitudinal_rebar_bot``, or from the constructor defaults. It is
    therefore readable at any time, unlike :class:`FlexureFaceDesign`, which
    also carries what a check demanded and so needs a check to have run.
    """

    layers: Tuple[RebarLayer, ...]
    A_s: Quantity

    @property
    def n_bars(self) -> int:
        """Total number of bars across every layer of this face."""
        return sum(layer.n for layer in self.layers)

    def __str__(self) -> str:
        if not self.layers:
            return "no reinforcement"
        return " + ".join(str(layer) for layer in self.layers)


#: How a section carries its transverse reinforcement. A beam carries closed
#: stirrups; a slab strip carries a grid of legs, with no cage to count.
STIRRUPS = "stirrups"
GRID = "grid"


def format_transverse_rebar(layout: str, n_stirrups: int, d_b: str, s_l: str, s_w: str) -> str:
    """Label the transverse reinforcement in the notation of its element.

    A beam is a number of closed stirrups of one diameter at one spacing along
    the length, so the count leads: ``2eØ10/15cm``. A slab strip has no cage --
    the same bar sits on a grid -- so what identifies it is the diameter once
    and a spacing each way, longitudinal first: ``Ø10/15cm×20cm``. The diameter
    is not repeated: both directions are the same bar.

    Takes the numbers already formatted, so each caller keeps its own precision
    and units while the shape of the label is decided in one place.
    """
    if n_stirrups == 0:
        return "no stirrups"
    if layout == GRID:
        return f"Ø{d_b}/{s_l}×{s_w}"
    return f"{n_stirrups}eØ{d_b}/{s_l}"


@dataclass(frozen=True)
class TransverseReinforcement:
    """The transverse reinforcement a section carries, as configured.

    ``s_l`` is the spacing along the length and ``s_w`` the spacing across the
    width. On a beam the second is not detailed but implied -- it is where the
    legs of the cage fall -- while on a slab strip both are chosen, which is
    what ``layout`` distinguishes.
    """

    n_stirrups: int
    d_b: Quantity
    s_l: Quantity
    A_v: Quantity
    s_w: Quantity
    layout: str = STIRRUPS

    @property
    def n_legs(self) -> int:
        """Number of stirrup legs crossing the shear plane."""
        return self.n_stirrups * 2

    def __str__(self) -> str:
        return format_transverse_rebar(
            self.layout,
            self.n_stirrups,
            f"{self.d_b:.4g~P}",
            f"{self.s_l:.4g~P}",
            f"{self.s_w:.4g~P}",
        )


@dataclass(frozen=True)
class SectionReinforcement:
    """Every bar a section carries: both faces plus the stirrups."""

    bottom: FaceReinforcement
    top: FaceReinforcement
    transverse: TransverseReinforcement

    def __str__(self) -> str:
        return f"bottom: {self.bottom} / top: {self.top} / stirrups: {self.transverse}"


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
    s_w: Quantity
    layout: str = STIRRUPS

    @property
    def n_legs(self) -> int:
        """Number of stirrup legs crossing the shear plane."""
        return self.n_stirrups * 2

    def __str__(self) -> str:
        return format_transverse_rebar(
            self.layout,
            self.n_stirrups,
            f"{self.d_b:.4g~P}",
            f"{self.s_l:.4g~P}",
            f"{self.s_w:.4g~P}",
        )


def transverse_layout(beam: RectangularBeam) -> str:
    """Which notation the section's transverse reinforcement is written in."""
    return GRID if getattr(beam, "mode", "beam") == "slab" else STIRRUPS


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
    # Envelope over every combination that was checked, computed here rather than
    # accumulated during the loop: the beam's private attributes only describe the
    # combination that ran last, which is not necessarily the one that governs.
    worst = envelope_flexure_face(getattr(beam, "_flexure_checks", ()), "bottom" if face == "b" else "top")

    # A_s is the reinforcement the section carries, not a per-combination result.
    A_s: Optional[Quantity] = getattr(beam, f"_A_s_{suffix}", None)

    return FlexureFaceDesign(
        layers=_layers(beam, face),
        A_s=zero if A_s is None else A_s,
        A_s_req=zero if worst.A_s_req is None else worst.A_s_req,
        A_s_min=zero if worst.A_s_min is None else worst.A_s_min,
        A_s_max=zero if worst.A_s_max is None else worst.A_s_max,
        DCR=worst.DCR,
    )


def build_reinforcement(beam: RectangularBeam) -> SectionReinforcement:
    """Build the public view of the reinforcement ``beam`` currently carries.

    Never raises: a section always has *some* reinforcement state, even if that
    state is "none". This is what separates it from the design results, which
    only exist once a check has run.
    """
    # These are set when the section is built, so they always exist.
    return SectionReinforcement(
        bottom=FaceReinforcement(layers=_layers(beam, "b"), A_s=beam._A_s_bot),
        top=FaceReinforcement(layers=_layers(beam, "t"), A_s=beam._A_s_top),
        transverse=TransverseReinforcement(
            n_stirrups=int(beam._stirrup_n),
            d_b=beam._stirrup_d_b,
            s_l=beam._stirrup_s_l,
            A_v=beam._A_v,
            s_w=beam._leg_spacing_across_width(),
            layout=transverse_layout(beam),
        ),
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
    # As in _face: the private attributes describe the combination that ran last,
    # so the envelope is taken over the results of every combination checked.
    worst = envelope_shear(getattr(beam, "_shear_checks", ()))

    return ShearDesign(
        n_stirrups=int(beam._stirrup_n),
        d_b=beam._stirrup_d_b,
        s_l=beam._stirrup_s_l,
        A_v=A_v,
        A_v_req=zero if worst.A_v_req is None else worst.A_v_req,
        A_v_min=zero if worst.A_v_min is None else worst.A_v_min,
        DCR=worst.DCR,
        s_w=beam._leg_spacing_across_width(),
        layout=transverse_layout(beam),
    )
