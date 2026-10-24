"""The shear report of a column: console tables, Word document and notebook lines.

Built from the column's results (mento/reports/columns.py); these tests hold
what it shows, which combination it picks and that every label it writes has
a Spanish entry.
"""

from pathlib import Path
from typing import Any, Iterator

import pytest

import mento
from mento import (
    CircularColumn,
    Concrete_ACI_318_19,
    Concrete_EN_1992_2004,
    DesignNotRunError,
    Forces,
    Node,
    RectangularColumn,
    SteelBar,
)
from mento.i18n import ES
from mento.reports import columns as reports
from mento.units import MPa, cm, kN, kNm, mm


@pytest.fixture(autouse=True)
def english() -> Iterator[None]:
    mento.set_language("en")
    yield
    mento.set_language("en")


FORCES = [
    Forces(label="1.2D+1.6L", N_x=300 * kN, V_y=60 * kN, V_z=150 * kN),
    Forces(label="0.9D+E", N_x=-200 * kN, V_z=180 * kN),
]


def _rectangle() -> RectangularColumn:
    column = RectangularColumn(
        label="C1",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=30 * cm,
        height=50 * cm,
        c_c=25 * mm,
    )
    column.set_longitudinal_rebar(d_b_corner=20 * mm, n_width=1, d_b_width=16 * mm, n_height=2, d_b_height=16 * mm)
    column.set_ties(d_b=10 * mm, s=20 * cm, crossties_height=1, d_b_crosstie=8 * mm)
    return column


def _circle() -> CircularColumn:
    column = CircularColumn(
        label="C2",
        concrete=Concrete_EN_1992_2004(name="C30", f_c=30 * MPa),
        steel_bar=SteelBar(name="B500S", f_y=500 * MPa),
        diameter=50 * cm,
        c_c=30 * mm,
    )
    column.set_longitudinal_rebar(n=8, d_b=20 * mm)
    column.set_hoops(d_b=10 * mm, s=15 * cm, spiral=True)
    column.set_lever_arm(d=40 * cm, z=30 * cm)
    return column


def test_the_report_shows_both_directions_of_the_governing_combination(capsys: pytest.CaptureFixture) -> None:
    """The tension combination governs (DCR 1.017 in z, against 0.641 of the other): its tables are printed."""
    column = _rectangle()
    column.check(FORCES)
    force, checks = reports.combination(column)
    assert force is FORCES[1]
    assert max(c.DCR for c in checks) == max(c.DCR for c in column.shear_checks)

    column.shear_results_detailed()
    out = capsys.readouterr().out
    assert "===== COLUMN SHEAR DETAILED RESULTS =====" in out
    assert "----- Shear along y -----" in out and "----- Shear along z -----" in out
    assert "0.5·ØVc (§10.6.2.1)" in out
    assert "-200" in out


def test_a_combination_can_be_asked_for(capsys: pytest.CaptureFixture) -> None:
    column = _rectangle()
    Node(column, FORCES).check()
    assert reports.combination(column, FORCES[0])[0] is FORCES[0]
    Node(column, FORCES).shear_results_detailed(FORCES[0])
    assert "300" in capsys.readouterr().out
    with pytest.raises(ValueError, match="No results"):
        column.shear_results_detailed(Forces(label="other", V_z=1 * kN))


def test_no_report_before_a_check() -> None:
    with pytest.raises(DesignNotRunError):
        _rectangle().shear_results_detailed()


def test_the_limit_table_holds_the_ties_to_what_the_shear_asks() -> None:
    """V_y of the first combination: the Ø10 every 20 cm passes A_v,req but not
    d/2 = 12.8 cm, nor the 25.6 cm across between the tie's side and the crosstie."""
    column = _rectangle()
    column.check(FORCES)
    table = reports.limits_table(column, column.shear_checks[0])
    verdicts = dict(zip(table["Check"], table["Ok?"]))
    assert verdicts == {
        "Shear reinforcement": "✅",
        "Tie spacing along the column": "❌",
        "Leg spacing across the section": "❌",
        "Maximum shear capacity": "✅",
    }
    strength = reports.strength_table(column, column.shear_checks[0])
    assert strength["Shear strength"][-1] == "Demand Capacity Ratio"
    assert strength["Value"][-1] == round(column.shear_checks[0].DCR, 3)


def test_an_en_circle_shows_its_symbols_and_its_lever_arm(capsys: pytest.CaptureFixture) -> None:
    column = _circle()
    column.check([Forces(label="ULS", N_x=1000 * kN, M_y=200 * kNm, V_z=250 * kN)])
    column.shear_results_detailed()
    out = capsys.readouterr().out
    assert "VRd,c (§6.2.1(3))" in out
    assert "Depth to the tension chord" in out and "Spiral" in out


def test_the_report_follows_the_language(capsys: pytest.CaptureFixture) -> None:
    column = _circle()
    column.check([Forces(label="ELU", N_x=1000 * kN, M_y=200 * kNm, V_z=250 * kN)])
    mento.set_language("es")
    column.shear_results_detailed()
    out = capsys.readouterr().out
    assert ES["===== COLUMN SHEAR DETAILED RESULTS ====="] in out
    assert "Corte según z" in out
    assert "Sí" in out


@pytest.mark.parametrize("language", ["en", "es"])
def test_the_word_report_is_saved(language: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    mento.set_language(language)
    column = _rectangle()
    column.check(FORCES)
    column.shear_results_detailed_doc(FORCES[0])
    saved = list(tmp_path.glob("*.docx"))
    assert [p.name for p in saved] == ["Column C1 shear check ACI 318-19.docx"]


def test_the_notebook_lines(monkeypatch: pytest.MonkeyPatch) -> None:
    shown: list[str] = []
    monkeypatch.setattr(reports, "_show", shown.append)
    column = _rectangle()
    column.data
    column.check(FORCES)
    column.results
    assert shown[0].startswith("Column C1: Section width 30 cm")
    assert "combination 0.9D+E" in shown[1] and "- V_z:" in shown[1]


def test_a_column_has_no_flexure_report() -> None:
    column = _rectangle()
    with pytest.raises(NotImplementedError, match="column-shear-roadmap"):
        column.flexure_results_detailed()


def _labels(table: dict) -> list[Any]:
    first = next(iter(table))
    return [*table.keys(), *table[first]]


def test_every_label_the_report_writes_is_in_the_spanish_catalog() -> None:
    rectangle, circle = _rectangle(), _circle()
    rectangle.set_ties(d_b=10 * mm, s=20 * cm, crossties_width=1, crossties_height=1, d_b_crosstie=8 * mm)
    rectangle.check(FORCES)
    circle.check([Forces(label="ULS", N_x=1000 * kN, M_y=200 * kNm, V_z=250 * kN)])
    labels = []
    for column in (rectangle, circle):
        force, checks = reports.combination(column)
        labels += _labels(reports.materials_table(column)) + _labels(reports.geometry_table(column))
        labels += _labels(reports.forces_table(column, force))
        for check in checks:
            labels += _labels(reports.limits_table(column, check)) + _labels(reports.strength_table(column, check))
    labels += ["Yes", "No", "===== COLUMN SHEAR DETAILED RESULTS =====", "----- Shear along {direction} -----"]
    labels += ["Concrete column shear check", "Column {label} shear check", "Shear along {direction}"]
    missing = sorted({label for label in labels if isinstance(label, str) and label and label not in ES})
    assert missing == []


def test_the_notebook_line_renders_and_bare_ties_show_as_a_dash() -> None:
    """The real Markdown display, outside a notebook, and a column with no ties."""
    column = _rectangle()
    column.set_ties(d_b=0 * mm, s=0 * cm)
    column.data
    geometry = reports.geometry_table(column)
    assert dict(zip(geometry["Geometry"], geometry["Value"]))["Ties"] == "-"


def test_a_column_has_no_flexure_document() -> None:
    with pytest.raises(NotImplementedError, match="column-shear-roadmap"):
        _rectangle().flexure_results_detailed_doc()
