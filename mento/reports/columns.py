"""Notebook views, console tables and Word reports of a column's shear.

Built from the column's results, :class:`~mento.column_results.ColumnShearCheck`,
and its inputs; nothing here computes a resistance or writes to the column
(ADR-0004). One combination at a time, both directions: the one that governs
by default, or the one asked for.

The labels are English and go through the catalog of :mod:`mento.i18n` when
printed or written; the symbols are the design code's own
(``DesignCode.column_shear_symbols``).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, cast

import pandas as pd
from IPython.display import Markdown, display

from mento._version import __version__ as MENTO_VERSION
from mento.bar_sizes import bar_designation
from mento.codes.registry import design_code
from mento.column_section import CircularColumn, RectangularColumn
from mento.i18n import get_language, translate
from mento.precompute import shown, unit_label
from mento.results import DocumentBuilder, Formatter, TablePrinter
from mento.units import Quantity, mm

if TYPE_CHECKING:
    from mento.column_results import ColumnShearCheck
    from mento.column_section import ColumnSection
    from mento.forces import Forces

Table = Dict[str, List[Any]]


def _show(markdown: str) -> None:
    """Render Markdown in a notebook; IPython ships no type information."""
    display(Markdown(markdown))  # type: ignore[no-untyped-call]


def _symbols(column: "ColumnSection") -> Dict[str, Any]:
    return design_code(column.concrete).column_shear_symbols


def _bar(d_mm: float, imperial: bool) -> str:
    """A bar as a drawing names it: ``Ø20`` in SI, ``#6`` in US customary."""
    return bar_designation(d_mm * mm) if imperial else f"Ø{d_mm:g}"


def _length(value: Optional[Quantity], imperial: bool) -> Any:
    return "" if value is None else shown(value, "length", imperial, 2)


def _ok(passes: bool) -> str:
    return "✅" if passes else "❌"


def _within(value: Quantity, limit: Quantity) -> bool:
    """``value`` at or under ``limit``, rounding aside."""
    limit = limit.to(value.units)
    return bool(value <= limit) or math.isclose(value.magnitude, limit.magnitude)


# ---------------------------------------------------------------------------
# Which combination
# ---------------------------------------------------------------------------


def combination(
    column: "ColumnSection", force: Optional["Forces"] = None
) -> Tuple["Forces", Tuple["ColumnShearCheck", ...]]:
    """The combination a report describes and its checks, ``y`` then ``z``.

    ``force`` picks one of those checked; without it, the one whose worse
    direction has the largest DCR. Raises ``DesignNotRunError`` before a
    check, and ``ValueError`` for a force that was not checked.
    """
    checks = column.shear_checks
    forces = column._checked_forces
    if force is None:
        worst = max(range(len(checks)), key=lambda k: checks[k].DCR)
        index = worst // 2
    else:
        found = [k for k, checked in enumerate(forces) if checked is force or checked.id == force.id]
        if not found:
            raise ValueError(f"No results found for Forces object with ID {force.id}.")
        index = found[0]
    return forces[index], checks[2 * index : 2 * index + 2]


# ---------------------------------------------------------------------------
# The tables
# ---------------------------------------------------------------------------


def materials_table(column: "ColumnSection") -> Table:
    imperial = column.concrete.is_imperial
    rows: List[Tuple[str, str, Any, str]] = [
        ("Section Label", "", column.label, ""),
        ("Concrete strength", "fc", shown(column.concrete.f_c, "stress", imperial, 2), unit_label("stress", imperial)),
        (
            "Steel reinforcement yield strength",
            "fy",
            shown(column.steel_bar.f_y, "steel_stress", imperial, 2),
            unit_label("steel_stress", imperial),
        ),
    ]
    # The factors each code applies: φ for shear and λ under ACI 318-19 and
    # CIRSOC 201-25, the partial factors under EN 1992-1-1.
    concrete: Any = column.concrete
    if hasattr(concrete, "phi_v"):
        rows += [
            ("Normalweight concrete", "λ", concrete.lambda_factor, ""),
            ("Safety factor for shear", "Øv", concrete.phi_v, ""),
        ]
    if hasattr(concrete, "gamma_c"):
        rows += [
            ("Partial factor for concrete", "γc", concrete.gamma_c, ""),
            ("Partial factor for reinforcement", "γs", column.steel_bar.gamma_s, ""),
        ]
    return _table("Materials", rows)


def geometry_table(column: "ColumnSection") -> Table:
    imperial = column.concrete.is_imperial
    length = unit_label("length", imperial)
    rows: List[Tuple[str, str, Any, str]] = []
    if isinstance(column, CircularColumn):
        circle = column
        rows += [
            ("Diameter", "D", shown(circle.diameter, "length", imperial, 2), length),
            ("Clear cover", "cc", shown(circle.c_c, "length", imperial, 2), length),
            ("Longitudinal rebar", "", f"{circle._n_bars}{_bar(circle._d_b, imperial)}", ""),
            ("Spiral" if circle._spiral else "Hoops", "", _ties(circle._tie_d_b, circle._tie_s, imperial), ""),
        ]
        if circle._lever_arm is not None:
            d, z = (value * mm for value in circle._lever_arm)
            rows += [
                ("Depth to the tension chord", "d", _length(d, imperial), length),
                ("Lever arm", "z", _length(z, imperial), length),
            ]
        return _table("Geometry", rows)
    rect = cast(RectangularColumn, column)
    corner = _bar(rect._d_b_corner, imperial)
    rows += [
        ("Section width", "b", shown(rect.width, "length", imperial, 2), length),
        ("Section height", "h", shown(rect.height, "length", imperial, 2), length),
        ("Clear cover", "cc", shown(rect.c_c, "length", imperial, 2), length),
        ("Corner bars", "", f"4{corner}", ""),
        ("Bars between corners, each face of the width", "", _bars(rect._n_width, rect._d_b_width, imperial), ""),
        ("Bars between corners, each face of the height", "", _bars(rect._n_height, rect._d_b_height, imperial), ""),
        ("Ties", "", _ties(rect._tie_d_b, rect._tie_s, imperial), ""),
        ("Crossties, faces of the width", "", rect._crossties_width, ""),
        ("Crossties, faces of the height", "", rect._crossties_height, ""),
    ]
    if rect._crossties_width or rect._crossties_height:
        rows.append(("Crosstie bar", "", _bar(rect._d_b_crosstie, imperial), ""))
    return _table("Geometry", rows)


def _bars(n: int, d_mm: float, imperial: bool) -> str:
    return f"{n}{_bar(d_mm, imperial)}" if n else "-"


def _ties(d_mm: float, s_mm: float, imperial: bool) -> str:
    """``Ø10 / 15 cm``: the bar and its spacing, or ``-`` with none."""
    if not s_mm:
        return "-"
    spacing = shown(s_mm * mm, "length", imperial, 2)
    return f"{_bar(d_mm, imperial)} / {spacing:g} {unit_label('length', imperial)}"


def forces_table(column: "ColumnSection", force: "Forces") -> Table:
    imperial = column.concrete.is_imperial
    symbols = _symbols(column)
    unit = unit_label("force", imperial)
    rows = [
        ("Axial, positive for compression", symbols["axial"], shown(force.N_x, "force", imperial, 2), unit),
        ("Shear along y", f"{symbols['demand']},y", shown(force.V_y, "force", imperial, 2), unit),
        ("Shear along z", f"{symbols['demand']},z", shown(force.V_z, "force", imperial, 2), unit),
    ]
    return _table("Design forces", rows)


def limits_table(column: "ColumnSection", check: "ColumnShearCheck") -> Table:
    """The ties against what this direction asks of them: area, spacings, and the section's limit."""
    imperial = column.concrete.is_imperial
    per_length = unit_label("per_length", imperial)
    length = unit_label("length", imperial)
    force = unit_label("force", imperial)
    digits = 3 if imperial else 2
    rows: List[List[Any]] = [
        [
            "Shear reinforcement",
            per_length,
            shown(check.A_v, "per_length", imperial, digits),
            shown(check.A_v_req, "per_length", imperial, digits),
            "",
            _ok(_within(check.A_v_req, check.A_v)),
        ]
    ]
    for label, value, limit in (
        ("Tie spacing along the column", check.s_l, check.s_max_l),
        ("Leg spacing across the section", check.s_w, check.s_max_w),
    ):
        if value is not None and limit is not None:
            rows.append(
                [label, length, _length(value, imperial), "", _length(limit, imperial), _ok(_within(value, limit))]
            )
    rows.append(
        [
            "Maximum shear capacity",
            force,
            shown(check.V_demand, "force", imperial, 2),
            "",
            shown(check.V_max, "force", imperial, 2),
            _ok(_within(check.V_demand, check.V_max)),
        ]
    )
    columns = ["Check", "Unit", "Value", "Min.", "Max.", "Ok?"]
    return {name: [row[k] for row in rows] for k, name in enumerate(columns)}


def strength_table(column: "ColumnSection", check: "ColumnShearCheck") -> Table:
    """The resistances of this direction, ending on its DCR."""
    imperial = column.concrete.is_imperial
    symbols = _symbols(column)
    force = unit_label("force", imperial)
    length = unit_label("length", imperial)
    threshold = check.V_concrete * symbols["threshold_share"]
    rows: List[Tuple[str, str, Any, str]] = [
        ("Shear width", "bw", shown(check.b_w, "length", imperial, 2), length),
        ("Effective height", "d", shown(check.d, "length", imperial, 2), length),
        ("Concrete shear strength", symbols["concrete"], shown(check.V_concrete, "force", imperial, 2), force),
        (
            "Shear reinforcement required past",
            f"{symbols['threshold']} ({symbols['required_clause']})",
            shown(threshold, "force", imperial, 2),
            force,
        ),
        ("Shear reinforcement required", "", translate("Yes" if check.reinforcement_required else "No"), ""),
        ("Total shear strength", symbols["capacity"], shown(check.V_capacity, "force", imperial, 2), force),
        ("Maximum shear strength", symbols["max"], shown(check.V_max, "force", imperial, 2), force),
        ("Demand Capacity Ratio", "DCR", round(check.DCR, 3), _ok(check.DCR <= 1)),
    ]
    return _table("Shear strength", rows)


def _table(title: str, rows: List[Tuple[str, str, Any, str]]) -> Table:
    return {
        title: [row[0] for row in rows],
        "Variable": [row[1] for row in rows],
        "Value": [row[2] for row in rows],
        "Unit": [row[3] for row in rows],
    }


# ---------------------------------------------------------------------------
# What the element shows
# ---------------------------------------------------------------------------


def column_data(column: "ColumnSection") -> None:
    """The column's inputs in one line of Markdown: its geometry, ties and materials."""
    geometry = geometry_table(column)
    parts = [
        f"{label} {value}{' ' + unit if unit else ''}"
        for label, value, unit in zip(geometry["Geometry"], geometry["Value"], geometry["Unit"])
    ]
    _show(
        f"Column {column.label}: "
        + ", ".join(parts)
        + f". Concrete {column.concrete.name}, rebar {column.steel_bar.name}."
    )


def column_shear_results(column: "ColumnSection") -> None:
    """The governing combination in one line of Markdown per direction."""
    force, checks = combination(column, None)
    imperial = column.concrete.is_imperial
    symbols = _symbols(column)
    unit = unit_label("force", imperial)
    formatter = Formatter()
    lines = [f"Column {column.label}, combination {force.label}:"]
    for check in checks:
        lines.append(
            f"- V_{check.direction}: {symbols['demand']} = {shown(check.V_demand, 'force', imperial, 2)} {unit}, "
            f"{symbols['capacity']} = {shown(check.V_capacity, 'force', imperial, 2)} {unit} → "
            f"{formatter.DCR(round(check.DCR, 3))}"
        )
    _show("\n".join(lines))


def column_shear_results_detailed(column: "ColumnSection", force: Optional["Forces"] = None) -> None:
    """The shear of one combination, both directions, as console tables."""
    force, checks = combination(column, force)
    language = get_language()
    print(translate("===== COLUMN SHEAR DETAILED RESULTS =====", language))
    TablePrinter("MATERIALS", language).print_table_data(materials_table(column), headers="keys")
    TablePrinter("GEOMETRY", language).print_table_data(geometry_table(column), headers="keys")
    TablePrinter("FORCES", language).print_table_data(forces_table(column, force), headers="keys")
    for check in checks:
        print(translate("----- Shear along {direction} -----", language, direction=check.direction))
        TablePrinter("MAX AND MIN LIMIT CHECKS", language).print_table_data(limits_table(column, check), headers="keys")
        TablePrinter("SHEAR STRENGTH", language).print_table_data(strength_table(column, check), headers="keys")


def column_shear_results_detailed_doc(column: "ColumnSection", force: Optional["Forces"] = None) -> None:
    """The same as a Word document, saved in the working directory."""
    force, checks = combination(column, force)
    doc = DocumentBuilder(title="Concrete column shear check", language=get_language())
    doc.add_heading("Column {label} shear check", level=1, label=column.label)
    doc.add_text(
        "Made with mento {version}. Design code: {design_code}",
        version=MENTO_VERSION,
        design_code=column.concrete.design_code,
    )
    doc.add_heading("Section Data", level=2)
    doc.add_table_data(pd.DataFrame(materials_table(column)))
    doc.add_table_data(pd.DataFrame(geometry_table(column)))
    doc.add_table_data(pd.DataFrame(forces_table(column, force)))
    for check in checks:
        doc.add_heading("Shear along {direction}", level=2, direction=check.direction)
        doc.add_table_min_max(pd.DataFrame(limits_table(column, check)))
        doc.add_table_dcr(pd.DataFrame(strength_table(column, check)))
    doc.save(f"Column {column.label} shear check {column.concrete.design_code}.docx")
