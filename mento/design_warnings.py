"""Structured warnings: the detailing limits a section does not meet.

The detailed reports mark these with a ❌ in their limit tables and a line of
text. A program reading a design needs them as data, so every check and design
also records them here::

    node.design()
    for warning in node.warnings:
        print(warning.code, warning.message, warning.values)

``code`` is stable -- compare against it, not against the message. ``message``
is written in the language set with :func:`mento.set_language` at the time the
warnings are read. ``values`` carries the numbers the message quotes, as pint
quantities.

A warning does not change a DCR. The DCR measures strength; these are the
limits around it -- minimum and maximum steel, spacing, bar fit -- which a
section can miss while still carrying its load, or meet while failing to.

Codes
-----
``As_below_min``
    A face carries less steel than its minimum, and the 4/3 relief of
    ACI 318-19 §9.6.1.3 / CIRSOC 201-25 §9.6.1.3 does not cover it. The
    minimum it quotes is the one left after that relief, ``A_s_min_eff``.
``As_above_max``
    A face carries more steel than its maximum. Under ACI 318-19 and CIRSOC
    201-25 that is the tension face past the tension-controlled limit of
    §9.3.3.1 with the compression steel it has, ``A_s_max_eff``: the section
    is no longer tension-controlled, and its capacity already carries the
    lower phi of the strain it reaches. Under EN 1992-1-1 it is either face
    past the 4 % of §9.2.1.1(3).
``clear_spacing_below_min``
    The clear distance between the bars of a beam face is below the minimum
    its settings ask for (bar diameter, 25 mm / 1 in., vibrator on top).
``bar_spacing_below_min`` / ``bar_spacing_exceeds_max``
    The same for a slab, which is detailed centre to centre: the spacing of
    the layer nearest the face against the code's limits.
``bars_do_not_fit``
    The bars on a face leave no clear space between them, or a design found no
    layout that fits the width -- which holds until the face is given bars by
    hand, since those are the spacing check's to judge.
``As_below_required``
    A design found no layout that fits the section and carries the moment --
    or, under ACI 318-19 / CIRSOC 201-25, carries it tension-controlled -- so
    it left the closest it found, and the section has to grow. It quotes the
    area the face was asked for and the one it was given, and stays while the
    face carries what the design left: bars set by hand afterwards are the
    check's to judge.
``stirrups_required``
    The section has no stirrups, and a combination asks for shear
    reinforcement -- beyond what the concrete carries, or the code minimum.
``Av_below_min``
    The stirrups provide less than the minimum shear reinforcement.
``stirrup_spacing_exceeds_max``
    The stirrups are further apart than the code allows, along the member or
    across its width.
``stirrup_diameter_below_min``
    The stirrup bar is thinner than the code's minimum.
``shear_exceeds_section_limit``
    The shear exceeds the most the section can carry however it is
    reinforced (ACI 318-19 §22.5.1.2, EN 1992-1-1 V_Rd,max): the section has
    to grow. A wall reports it against ØVn,max of §11.5.4.3.
``mesh_ratio_below_min``
    A wall mesh gives less than its direction asks for: the horizontal one
    below the ρt the shear needs (never below its minimum), the vertical one
    below ρl,min of ACI 318-19 / CIRSOC 201-25 §11.6.2. ``values`` carries
    ``direction``, ``"h"`` or ``"v"``.
``mesh_spacing_exceeds_max``
    The bars of a wall mesh are further apart than §11.7 allows.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Optional, Tuple

from mento.codes.check_state import to_display
from mento.codes.registry import design_code
from mento.i18n import translate
from mento.precompute import DISPLAY
from mento.units import Quantity, inch, mm

if TYPE_CHECKING:
    from mento.beam import RectangularBeam
    from mento.wall_results import WallMesh, WallShearCheck


@dataclass(frozen=True)
class DesignWarning:
    """One detailing limit the section does not meet.

    ``code`` is the stable identifier (``"stirrup_spacing_exceeds_max"``);
    ``message`` the text, in the language set with :func:`mento.set_language`;
    ``values`` the numbers the message quotes, by name (``s``, ``s_max``,
    ``A_s``, ``A_s_min`` ...), as quantities. ``face`` is ``"bottom"`` or
    ``"top"`` for a longitudinal warning and ``None`` otherwise, and
    ``combinations`` holds the labels of the load combinations the limit is
    missed under -- empty for a limit of the section alone, such as the bar
    spacing or the stirrup diameter.
    """

    code: str
    message: str
    values: Mapping[str, Any] = field(default_factory=dict)
    face: Optional[str] = None
    combinations: Tuple[str, ...] = ()

    def __str__(self) -> str:
        return self.message


@dataclass
class _Raw:
    """A warning as the check found it, before it is worded.

    Kept apart from :class:`DesignWarning` so the message is written when it
    is read, in the language in force then, and so the same limit missed under
    several combinations collapses into one warning.
    """

    code: str
    values: Dict[str, Any]
    face: Optional[str] = None
    combination: Optional[str] = None
    #: How far past its limit this occurrence is; of the same limit missed under
    #: several combinations, the furthest is the one reported.
    severity: float = 0.0


#: The English wording of each code; the text is also the key of the Spanish
#: catalog in :mod:`mento.i18n`. ``{face}`` is filled with the translated face.
_MESSAGES: Dict[str, str] = {
    "As_below_min": "Steel on the {face}: A_s = {A_s} is below the minimum A_s,min = {A_s_min}.",
    "As_above_max": "Steel on the {face}: A_s = {A_s} exceeds the maximum A_s,max = {A_s_max}.",
    "clear_spacing_below_min": "Clear spacing between the bars on the {face}: {s} is below the minimum {s_min}.",
    "bar_spacing_below_min": "Bar spacing on the {face}: {s} is below the minimum {s_min}.",
    "bar_spacing_exceeds_max": "Bar spacing on the {face}: {s} exceeds the maximum {s_max}.",
    "bars_do_not_fit": "The bars on the {face} do not fit in the width of the section.",
    "As_below_required": (
        "Steel on the {face}: no layout that fits the section carries the moment "
        "(A_s,req = {A_s_req}); the design left A_s = {A_s}. Enlarge the section."
    ),
    "stirrups_required": "The section has no stirrups and requires shear reinforcement A_v = {A_v_req}.",
    "Av_below_min": "The stirrups provide A_v = {A_v}, below the minimum A_v,min = {A_v_min}.",
    "stirrup_spacing_exceeds_max_l": "Stirrup spacing along the member: {s} exceeds the maximum {s_max}.",
    "stirrup_spacing_exceeds_max_w": "Stirrup leg spacing across the width: {s} exceeds the maximum {s_max}.",
    "stirrup_diameter_below_min": "Stirrup diameter {d_b} is below the minimum {d_b_min}.",
    "shear_exceeds_section_limit": "Shear V = {V} exceeds the most the section can carry, {V_max}: enlarge the section.",
    "mesh_ratio_below_min_h": "Horizontal wall mesh: ρt = {rho} is below the required ρt = {rho_min}.",
    "mesh_ratio_below_min_v": "Vertical wall mesh: ρl = {rho} is below the minimum ρl,min = {rho_min}.",
    "mesh_spacing_exceeds_max_h": "Horizontal wall mesh spacing: {s} exceeds the maximum {s_max}.",
    "mesh_spacing_exceeds_max_v": "Vertical wall mesh spacing: {s} exceeds the maximum {s_max}.",
}

_FACES = {"bottom": "bottom face", "top": "top face"}


def _format(value: Any) -> str:
    """A value as the messages print it: three significant figures.

    Everything a warning quotes is a quantity -- an area, a spacing, a force --
    so there is one format, and pint's ``~P`` writes the unit with it. A
    reinforcement ratio is the exception, a bare number.
    """
    if not isinstance(value, Quantity):
        return f"{value:.3g}"
    return f"{value:.3g~P}"


def _q(value: float, kind: str, beam: "RectangularBeam") -> Quantity:
    """A float of a check state back as a quantity in the section's display units."""
    return to_display(value, kind, beam.concrete.is_imperial)


def _shown(values: Dict[str, Any], beam: "RectangularBeam") -> Dict[str, Any]:
    """The values in the units the detailed reports print them in.

    Chosen by name: ``A_s*`` is an area, ``A_v*`` an area per length, ``s*`` a
    spacing, ``d_b*`` a bar diameter, ``V*`` a force. Anything else is left as
    it is.
    """
    units = DISPLAY[beam.concrete.is_imperial]
    kinds: Tuple[Tuple[str, Any], ...] = (
        ("A_s", units["area"]),
        ("A_v", units["per_length"]),
        ("d_b", inch if beam.concrete.is_imperial else mm),
        ("s", units["length"]),
        ("V", units["force"]),
    )
    shown = {}
    for name, value in values.items():
        unit = next((unit for prefix, unit in kinds if name.startswith(prefix)), None)
        shown[name] = value.to(unit) if isinstance(value, Quantity) and unit is not None else value
    return shown


def _with_units(raw: _Raw, beam: "RectangularBeam") -> _Raw:
    raw.values = _shown(raw.values, beam)
    return raw


def _face_name(suffix: str) -> str:
    return "bottom" if suffix in ("bot", "b") else "top"


# ---------------------------------------------------------------------------
# What a flexure check leaves
# ---------------------------------------------------------------------------


def flexure_warnings(beam: "RectangularBeam", label: str, state: Any) -> List[_Raw]:
    """The steel-area limits one flexure combination finds on each face.

    Mirrors the limit rows of the detailed report: a face below its minimum
    warns unless the 4/3 relief of ACI 318-19 §9.6.1.3 covers it, and a face
    above its maximum warns.

    Which maximum, and on which face, is the code's. Where it is the
    ductility limit of the tension steel (ACI 318-19 / CIRSOC 201-25
    §9.3.3.1) only the face the combination puts in tension is held to it --
    the bars a negative moment asks for on the bottom are compression steel,
    and a combination with no moment pulls neither face -- and the limit is
    ``A_s_max_eff``, which the compression steel opposite extends: a doubly
    reinforced face is judged by the strain it reaches, not excused. Where it
    caps any bars (EN 1992-1-1 §9.2.1.1(3)) both faces are read against it.
    """
    found: List[_Raw] = []
    ductility_limit = design_code(beam.concrete).max_steel_is_ductility_limit
    doubly = bool(getattr(state, "doubly_reinforced", False))
    M = getattr(state, "M_u", getattr(state, "M_Ed", 0.0))
    tension_face = "bot" if M > 0 else "top" if M < 0 else None
    for suffix in ("bot", "top"):
        A_s: Quantity = getattr(beam, f"_A_s_{suffix}")
        # The minimum the face has to meet, relief of §9.6.1.3 included; a
        # code with no relief (EN 1992-1-1) leaves it at A_s_min.
        A_s_min_raw = getattr(state, f"A_s_min_{suffix}")
        A_s_min = _q(getattr(state, f"A_s_min_eff_{suffix}", A_s_min_raw), "area", beam).to(A_s.units)
        if ductility_limit:
            limit = getattr(state, f"A_s_max_eff_{suffix}")
            applies = suffix == tension_face
        else:
            limit = getattr(state, f"A_s_max_{suffix}")
            applies = not doubly
        A_s_max = _q(limit, "area", beam).to(A_s.units)
        if A_s < A_s_min and not math.isclose(A_s.magnitude, A_s_min.magnitude):
            found.append(
                _Raw(
                    "As_below_min",
                    {"A_s": A_s, "A_s_min": A_s_min},
                    _face_name(suffix),
                    label,
                    severity=float((A_s_min - A_s).magnitude),
                )
            )
        if applies and A_s_max.magnitude > 0 and A_s > A_s_max and not math.isclose(A_s.magnitude, A_s_max.magnitude):
            found.append(
                _Raw(
                    "As_above_max",
                    {"A_s": A_s, "A_s_max": A_s_max},
                    _face_name(suffix),
                    label,
                    severity=float((A_s - A_s_max).magnitude),
                )
            )
    return [_with_units(raw, beam) for raw in found]


def spacing_warnings(beam: "RectangularBeam") -> List[_Raw]:
    """The bar-spacing limits of the section itself, whatever the forces.

    Uses the same rows as the detailed report
    (:func:`mento.reports.tables._bar_spacing_row`): the clear distance of a
    beam face against its minimum, the centre-to-centre spacing of a slab
    against both of its limits.
    """
    from mento.reports.tables import _bar_spacing_row

    settings = beam.settings
    assert settings is not None
    found: List[_Raw] = []
    for face, suffix in (("t", "top"), ("b", "bot")):
        name = _face_name(suffix)
        layers = beam.reinforcement.top.layers if face == "t" else beam.reinforcement.bottom.layers
        if not layers:
            continue
        d_b_max = max(layer.d_b for layer in layers)
        minimum = max(settings.clear_spacing, d_b_max)
        if face == "t":
            minimum = max(minimum, settings.vibrator_size)
        _label, value, s_min, s_max = _bar_spacing_row(beam, face, minimum)
        is_slab = getattr(beam, f"_s_b1_{face}", None) is not None
        if not is_slab and value.magnitude <= 0 and sum(layer.n for layer in layers) > 1:
            found.append(_Raw("bars_do_not_fit", {"s": value}, name))
            continue
        if s_min is not None and value < s_min and not math.isclose(value.magnitude, s_min.to(value.units).magnitude):
            code = "bar_spacing_below_min" if is_slab else "clear_spacing_below_min"
            found.append(_Raw(code, {"s": value, "s_min": s_min}, name))
        if s_max is not None and value > s_max and not math.isclose(value.magnitude, s_max.to(value.units).magnitude):
            found.append(_Raw("bar_spacing_exceeds_max", {"s": value, "s_max": s_max}, name))
    for face in sorted(getattr(beam, "_infeasible_faces", ())):
        name = _face_name(face)
        if not any(raw.code == "bars_do_not_fit" and raw.face == name for raw in found):
            found.append(_Raw("bars_do_not_fit", {}, name))
    return [_with_units(raw, beam) for raw in found]


def shortfall_warnings(beam: "RectangularBeam") -> List[_Raw]:
    """The faces the last design could not bring up to what they need.

    The design records, per face, the area it asked for and the one it left.
    The warning holds while the face still carries that one: bars changed by
    hand afterwards are a different section, which the check judges.
    """
    found: List[_Raw] = []
    for suffix, (needed, placed) in sorted(getattr(beam, "_short_faces", {}).items()):
        A_s: Quantity = getattr(beam, f"_A_s_{suffix}")
        if not math.isclose(A_s.magnitude, placed.to(A_s.units).magnitude):
            continue
        found.append(_Raw("As_below_required", {"A_s": A_s, "A_s_req": needed.to(A_s.units)}, _face_name(suffix)))
    return [_with_units(raw, beam) for raw in found]


# ---------------------------------------------------------------------------
# What a shear check leaves
# ---------------------------------------------------------------------------


def shear_warnings(beam: "RectangularBeam", label: str, state: Any) -> List[_Raw]:
    """The transverse-reinforcement limits one shear combination finds.

    Mirrors the limit rows of the detailed report -- spacing along and across,
    the minimum area and, where the code states one, the minimum diameter --
    and adds the two a section without stirrups or with too small a web runs
    into.
    """
    found: List[_Raw] = []
    A_v: Quantity = beam._A_v
    A_v_min = _q(state.A_v_min, "per_length", beam).to(A_v.units)
    A_v_req = _q(state.A_v_req, "per_length", beam).to(A_v.units)

    if not state.max_shear_ok:
        V = _q(state.V_u if hasattr(state, "V_u") else state.V_Ed_1, "force", beam)
        V_max = _q(state.phi_V_max if hasattr(state, "phi_V_max") else state.V_Rd_max, "force", beam)
        found.append(
            _Raw(
                "shear_exceeds_section_limit",
                {"V": V, "V_max": V_max},
                None,
                label,
                severity=float((V - V_max).magnitude),
            )
        )

    if int(beam._stirrup_n) == 0 or A_v.magnitude == 0:
        if A_v_req.magnitude > 0:
            found.append(
                _Raw(
                    "stirrups_required",
                    {"A_v_req": A_v_req, "A_v_min": A_v_min},
                    None,
                    label,
                    severity=float(A_v_req.magnitude),
                )
            )
        return [_with_units(raw, beam) for raw in found]

    if A_v < A_v_min and not math.isclose(A_v.magnitude, A_v_min.magnitude):
        found.append(
            _Raw(
                "Av_below_min",
                {"A_v": A_v, "A_v_min": A_v_min},
                None,
                label,
                severity=float((A_v_min - A_v).magnitude),
            )
        )

    s_l: Quantity = beam._stirrup_s_l
    s_w: Quantity = beam._leg_spacing_across_width()
    s_max_l = _q(state.stirrup_s_max_l, "length", beam).to(s_l.units)
    s_max_w = _q(state.stirrup_s_max_w, "length", beam).to(s_l.units)
    for direction, s, s_max in (("l", s_l, s_max_l), ("w", s_w.to(s_l.units), s_max_w)):
        if s_max.magnitude > 0 and s > s_max and not math.isclose(s.magnitude, s_max.magnitude):
            found.append(
                _Raw(
                    "stirrup_spacing_exceeds_max",
                    {"s": s, "s_max": s_max, "direction": direction},
                    None,
                    label,
                    severity=float((s - s_max).magnitude),
                )
            )

    minimum = design_code(beam.concrete).min_stirrup_diameter
    if minimum is not None:
        d_b_min: Quantity = minimum(beam.concrete)
        d_b: Quantity = beam._stirrup_d_b
        if d_b < d_b_min:
            found.append(_Raw("stirrup_diameter_below_min", {"d_b": d_b, "d_b_min": d_b_min.to(d_b.units)}, None))
    return [_with_units(raw, beam) for raw in found]


# ---------------------------------------------------------------------------
# What a wall shear check leaves
# ---------------------------------------------------------------------------


def wall_warnings(wall: "RectangularBeam", mesh: "WallMesh", checks: Tuple["WallShearCheck", ...]) -> List[_Raw]:
    """The mesh limits a wall misses, over every combination checked.

    Mirrors the limit rows of the wall report: each ratio against what its
    direction asks for, each spacing against its maximum, and the shear
    against the most the section can carry. A mesh with a zero spacing has no
    bars, so its spacing is not a limit it misses; its ratio is.
    """
    found: List[_Raw] = []
    for check in checks:
        for direction, provided, required in (
            ("h", mesh.horizontal, check.rho_t_req),
            ("v", mesh.vertical, check.rho_l_min),
        ):
            if provided.rho < required and not math.isclose(provided.rho, required, rel_tol=1e-9):
                values = {"direction": direction, "rho": round(provided.rho, 5), "rho_min": round(required, 5)}
                found.append(_Raw("mesh_ratio_below_min", values, None, check.label, required - provided.rho))
        for direction, provided, s_max in (
            ("h", mesh.horizontal, check.s_h_max),
            ("v", mesh.vertical, check.s_v_max),
        ):
            s_max = s_max.to(provided.s.units)
            if provided.has_bars and provided.s > s_max and not math.isclose(provided.s.magnitude, s_max.magnitude):
                values = {"direction": direction, "s": provided.s, "s_max": s_max}
                found.append(
                    _Raw("mesh_spacing_exceeds_max", values, None, None, float((provided.s - s_max).magnitude))
                )
        if check.V_u > check.V_max:
            values = {"V": check.V_u, "V_max": check.V_max}
            severity = float((check.V_u - check.V_max).magnitude)
            found.append(_Raw("shear_exceeds_section_limit", values, None, check.label, severity))
    return [_with_units(raw, wall) for raw in found]


# ---------------------------------------------------------------------------
# Wording
# ---------------------------------------------------------------------------


def collect(raws: List[_Raw]) -> Tuple[DesignWarning, ...]:
    """Collapse the raw findings into one worded warning per limit and face.

    The same limit missed under several combinations is one warning, with the
    values of the combination that misses it by most and the labels of all of
    them.
    """
    groups: Dict[Tuple[str, Optional[str], Optional[str]], List[_Raw]] = {}
    for raw in raws:
        key = (raw.code, raw.face, raw.values.get("direction"))
        groups.setdefault(key, []).append(raw)

    warnings: List[DesignWarning] = []
    for (code, face, direction), group in groups.items():
        worst = max(group, key=lambda raw: raw.severity)
        labels = tuple(dict.fromkeys(raw.combination for raw in group if raw.combination is not None))
        values = {name: value for name, value in worst.values.items() if name != "direction"}
        template = _MESSAGES[f"{code}_{direction}" if direction else code]
        fields = {name: _format(value) for name, value in values.items()}
        if face is not None:
            fields["face"] = translate(_FACES[face])
        warnings.append(
            DesignWarning(
                code=code,
                message=translate(template, **fields),
                values=MappingProxyType(values),
                face=face,
                combinations=labels,
            )
        )
    return tuple(warnings)
