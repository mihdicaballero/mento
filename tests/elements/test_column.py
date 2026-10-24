"""RectangularColumn and CircularColumn: the shear of a column section in each direction.

The plan is docs/architecture/column-shear-roadmap.md, Phase B. The expected
values are worked by hand in each docstring; the equivalence tests hold the
column to RectangularBeam on the same rectangle, which is what the roadmap
asks for the clauses the two share.
"""

import math
from typing import Any, Iterator

import pytest

import mento
from mento import (
    CircularColumn,
    ColumnShearCheck,
    Concrete_ACI_318_19,
    Concrete_CIRSOC_201_25,
    Concrete_EN_1992_2004,
    DesignNotRunError,
    Forces,
    Node,
    RectangularBeam,
    RectangularColumn,
    SteelBar,
)
from mento.i18n import ES
from mento.units import MPa, cm, inch, kip, kN, kNm, ksi, mm, psi

STEEL = SteelBar(name="ADN 420", f_y=420 * MPa)


@pytest.fixture(autouse=True)
def english() -> Iterator[None]:
    mento.set_language("en")
    yield
    mento.set_language("en")


def _column(concrete: Any = None, width: Any = 30 * cm, height: Any = 50 * cm) -> RectangularColumn:
    column = RectangularColumn(
        label="C1",
        concrete=concrete or Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=STEEL,
        width=width,
        height=height,
        c_c=25 * mm,
    )
    column.set_longitudinal_rebar(d_b_corner=20 * mm, n_width=1, d_b_width=16 * mm, n_height=2, d_b_height=16 * mm)
    return column


def _z(checks: tuple[ColumnShearCheck, ...]) -> ColumnShearCheck:
    return next(check for check in checks if check.direction == "z")


def _y(checks: tuple[ColumnShearCheck, ...]) -> ColumnShearCheck:
    return next(check for check in checks if check.direction == "y")


# ---------------------------------------------------------------------------
# The same numbers as the beam where the clauses are the same
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("concrete_class", [Concrete_ACI_318_19, Concrete_CIRSOC_201_25], ids=["ACI", "CIRSOC"])
@pytest.mark.parametrize("N_kN", [300.0, 0.0, -300.0], ids=["compression", "no-axial", "tension"])
@pytest.mark.parametrize("V_kN", [50.0, 150.0, 400.0])
@pytest.mark.parametrize("ties", [True, False], ids=["ties", "bare"])
def test_the_z_direction_is_the_beam_of_the_same_rectangle(
    concrete_class: type, N_kN: float, V_kN: float, ties: bool
) -> None:
    """V_z of a 30x50 column, the corners Ø20 and one Ø16 between them on the
    faces of length width, is the shear of a 30x50 beam with 2Ø20 + 1Ø16 on
    each face and the same tie: V_c with N_u/(6·A_g), V_s, the section limit
    and the spacing rows of Table 9.7.6.2.2 that Table 10.7.6.5.2 sends to.
    The area asked for is compared past both thresholds, where the column's
    §10.6.2.1 and the beam's §9.6.3.1 agree on requiring it."""
    concrete = concrete_class(name="H25", f_c=25 * MPa)
    column = _column(concrete)
    beam = RectangularBeam(label="B", concrete=concrete, steel_bar=STEEL, width=30 * cm, height=50 * cm, c_c=25 * mm)
    beam.set_longitudinal_rebar_bot(2, 20 * mm, 1, 16 * mm)
    beam.set_longitudinal_rebar_top(2, 20 * mm, 1, 16 * mm)
    if ties:
        column.set_ties(d_b=10 * mm, s=20 * cm)
        beam.set_transverse_rebar(d_b=10 * mm, s_l=20 * cm, legs=2)
    else:
        # The beam assumes a starting stirrup in its depth; the column with no
        # ties has none, so the beam is told so to stand for it.
        beam.set_transverse_rebar(0, 0 * mm, 0 * cm)
    force = Forces(label="U", N_x=N_kN * kN, V_z=V_kN * kN)

    (on_beam,) = beam.shear_check_results([force])
    on_column = _z(column.shear_check_results([force]))

    assert on_column.V_concrete.to("kN").magnitude == pytest.approx(on_beam.V_concrete.to("kN").magnitude, rel=1e-9)
    assert on_column.V_capacity.to("kN").magnitude == pytest.approx(on_beam.V_capacity.to("kN").magnitude, rel=1e-9)
    assert on_column.DCR == pytest.approx(on_beam.DCR, rel=1e-9)
    if V_kN >= 150 and ties:
        assert on_column.reinforcement_required
        assert on_column.A_v_req.to("cm**2/m").magnitude == pytest.approx(
            on_beam.A_v_req.to("cm**2/m").magnitude, rel=1e-9
        )
        assert on_column.s_max_l.to("cm").magnitude == pytest.approx(on_beam.s_max_l_table.to("cm").magnitude)
        assert on_column.s_max_w.to("cm").magnitude == pytest.approx(on_beam.s_max_w.to("cm").magnitude)


def test_the_y_direction_is_the_beam_turned() -> None:
    """V_y reads the column turned: a 50x30 beam with 2Ø20 + 2Ø16 on each face."""
    column = _column()
    column.set_ties(d_b=10 * mm, s=20 * cm)
    beam = RectangularBeam(
        label="B", concrete=column.concrete, steel_bar=STEEL, width=50 * cm, height=30 * cm, c_c=25 * mm
    )
    beam.set_longitudinal_rebar_bot(2, 20 * mm, 2, 16 * mm)
    beam.set_longitudinal_rebar_top(2, 20 * mm, 2, 16 * mm)
    beam.set_transverse_rebar(d_b=10 * mm, s_l=20 * cm, legs=2)

    (on_beam,) = beam.shear_check_results([Forces(label="U", N_x=200 * kN, V_z=180 * kN)])
    on_column = _y(column.shear_check_results([Forces(label="U", N_x=200 * kN, V_y=180 * kN)]))

    assert on_column.DCR == pytest.approx(on_beam.DCR, rel=1e-9)
    assert on_column.A_v_req.to("cm**2/m").magnitude == pytest.approx(on_beam.A_v_req.to("cm**2/m").magnitude)


# ---------------------------------------------------------------------------
# By hand
# ---------------------------------------------------------------------------


def test_a_compressed_rectangle_with_two_diameters_of_leg_by_hand() -> None:
    """30x50, H25, cover 25 mm, Ø20 corners and 1Ø16 between them on the faces
    of length 30; tie Ø10 every 20 cm with one Ø6 crosstie on those faces.
    N_u = 300 kN, V_z = 150 kN.

    V_z: b_w = 300, far face 2Ø20 + 1Ø16, A_s = 628.3 + 201.1 = 829.4 mm²,
    centroid (628.3·10 + 201.1·8)/829.4 = 9.52 mm, d = 500 - 25 - 10 - 9.52 =
    455.5 mm. The legs: the tie's two sides and the Ø6, A_v = (2·78.54 +
    28.27)/200 = 0.927 mm²/mm, past A_v,min = 0.35/420·300 = 0.25 (Table
    9.6.3.4), so Table 22.5.5.1(a): v_c = 0.17·5 + 300000/(6·150000) =
    1.183 MPa, ØVc = 0.75·1.183·300·455.5 = 121.3 kN. 150 kN > 0.5·ØVc: shear
    reinforcement required (§10.6.2.1). V_s,req = (150 - 121.3)/0.75 = 38.3 kN,
    A_v,req = max(38.3e3/(420·455.5), 0.25) = 0.25 mm²/mm = 2.50 cm²/m.
    ØVn = 0.75·(161.7 + 0.927·420·455.5) = 254.3 kN.
    """
    column = _column()
    column.set_ties(d_b=10 * mm, s=20 * cm, crossties_width=1, d_b_crosstie=6 * mm)
    check = _z(column.shear_check_results([Forces(label="U", N_x=300 * kN, V_z=150 * kN)]))

    A_s = 2 * math.pi * 20**2 / 4 + math.pi * 16**2 / 4
    d = 500 - 25 - 10 - (2 * math.pi * 20**2 / 4 * 10 + math.pi * 16**2 / 4 * 8) / A_s
    A_v = (2 * math.pi * 10**2 / 4 + math.pi * 6**2 / 4) / 200
    V_c = (0.17 * 5 + 300e3 / (6 * 300 * 500)) * 300 * d
    assert check.d.to("mm").magnitude == pytest.approx(d)
    assert check.A_v.to("cm**2/m").magnitude == pytest.approx(A_v * 10)
    assert check.V_concrete.to("kN").magnitude == pytest.approx(0.75 * V_c / 1e3)
    assert check.reinforcement_required
    assert check.A_v_req.to("cm**2/m").magnitude == pytest.approx(0.35 / 420 * 300 * 10)
    assert check.V_capacity.to("kN").magnitude == pytest.approx(0.75 * (V_c + A_v * 420 * d) / 1e3)
    assert check.DCR == pytest.approx(150e3 / (0.75 * (V_c + A_v * 420 * d)))


def test_a_compressed_column_needs_no_shear_reinforcement_below_half_phi_V_c() -> None:
    """§10.6.2.1 against §9.6.3.1. The column of the test above with a tie Ø8
    every 25 cm (A_v = 0.402 mm²/mm, past A_v,min = 0.25: Table 22.5.5.1(a)),
    d = 500 - 25 - 8 - 9.52 = 457.5 mm, under N_u = 300 kN, V_z = 55 kN:
    ØVc = 0.75·1.183·300·457.5 = 121.8 kN and 0.5·ØVc = 60.9 kN, so the column
    asks for no shear reinforcement. A beam of the same rectangle is held to
    φ·0.083·√f'c·b_w·d = 0.75·0.083·5·300·457.5 = 42.7 kN, which ignores the
    axial load, and asks it for A_v,min."""
    column = _column()
    column.set_ties(d_b=8 * mm, s=25 * cm)
    beam = RectangularBeam(
        label="B", concrete=column.concrete, steel_bar=STEEL, width=30 * cm, height=50 * cm, c_c=25 * mm
    )
    beam.set_longitudinal_rebar_bot(2, 20 * mm, 1, 16 * mm)
    beam.set_longitudinal_rebar_top(2, 20 * mm, 1, 16 * mm)
    beam.set_transverse_rebar(d_b=8 * mm, s_l=25 * cm, legs=2)
    force = Forces(label="U", N_x=300 * kN, V_z=55 * kN)

    check = _z(column.shear_check_results([force]))
    (on_beam,) = beam.shear_check_results([force])

    assert not check.reinforcement_required
    assert check.A_v_req.magnitude == 0 and check.A_v_min.magnitude == 0
    assert check.s_max_l is None and check.s_max_w is None
    assert on_beam.A_v_req.magnitude > 0
    assert column.warnings == ()


def test_a_circle_in_tension_by_hand() -> None:
    """Ø50, H25, 8Ø20, hoop Ø10 every 15 cm, N_u = -500 kN, V = 200 kN.

    §22.5.2.2: b_w = 500, d = 400 mm. A_g = π·500²/4 = 196 350 mm². The bars
    past the centre line, (8 - 1)//2 = 3Ø20, A_s = 942.5 mm². The hoop is two
    legs: A_v = 2·78.54/150 = 1.047 mm²/mm, past A_v,min = 0.35/420·500 = 0.417.
    Table 22.5.5.1(a): 0.17·5 - 500000/(6·196350) = 0.425 MPa; (b): 0.66 ·
    (942.5/200000)^(1/3)·5 - 0.424 = 0.129 MPa. ØVc = 0.75·0.425·500·400 =
    63.8 kN, the same in both directions.
    """
    column = CircularColumn(
        label="C2",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=STEEL,
        diameter=50 * cm,
        c_c=30 * mm,
    )
    column.set_longitudinal_rebar(n=8, d_b=20 * mm)
    column.set_hoops(d_b=10 * mm, s=15 * cm)
    checks = column.shear_check_results([Forces(label="T", N_x=-500 * kN, V_y=200 * kN, V_z=200 * kN)])

    v_c = 0.17 * 5 - 500e3 / (6 * math.pi * 500**2 / 4)
    for check in checks:
        assert check.b_w.to("mm").magnitude == pytest.approx(500)
        assert check.d.to("mm").magnitude == pytest.approx(400)
        assert check.A_v.to("cm**2/m").magnitude == pytest.approx(2 * math.pi * 10**2 / 4 / 150 * 10)
        assert check.V_concrete.to("kN").magnitude == pytest.approx(0.75 * v_c * 500 * 400 / 1e3)
        assert check.s_w is None and check.s_max_w is None
        assert math.isfinite(check.DCR)


def test_a_column_keeps_the_cap_on_the_root_of_f_c() -> None:
    """§22.5.3.1 caps √f'c at 8.3 MPa for V_c and §22.5.3.2 lifts it for beams
    and joists only. A 30x50 column of H80 with ties past A_v,min, N = 0:
    v_c = 0.17·8.3 (Table 22.5.5.1(a)), not 0.17·√80."""
    column = _column(Concrete_ACI_318_19(name="H80", f_c=80 * MPa))
    column.set_ties(d_b=10 * mm, s=10 * cm)
    check = _z(column.shear_check_results([Forces(label="U", V_z=100 * kN)]))
    d = check.d.to("mm").magnitude
    assert check.V_concrete.to("kN").magnitude == pytest.approx(0.75 * 0.17 * 8.3 * 300 * d / 1e3)


def test_the_cirsoc_column_takes_its_own_spacing_caps() -> None:
    """Table 10.7.6.5.2 is Table 9.7.6.2.2: d/2 along, capped at 600 mm (ACI
    318-19) or 400 mm (CIRSOC 201-25). A 100x100 column has d/2 near 475 mm."""
    force = Forces(label="U", V_z=1500 * kN)
    caps = []
    for concrete in (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa)):
        column = _column(concrete, width=100 * cm, height=100 * cm)
        column.set_ties(d_b=12 * mm, s=15 * cm)
        caps.append(_z(column.shear_check_results([force])).s_max_l.to("mm").magnitude)
    d = 1000 - 25 - 12 - (2 * 400 * 10 + 256 * 8) / (2 * 400 + 256)
    assert caps[0] == pytest.approx(d / 2)
    assert caps[1] == pytest.approx(400)


def test_the_legs_across_are_the_tie_and_the_crossties_at_their_bars() -> None:
    """V_y of the 30x50 reads the faces of length 50: two Ø16 between the Ø20
    corners, their centres at ±(250 - 25 - 10 - 10)/3 = ±68.3 mm, and one
    crosstie on the second of them. The tie's sides at ±(250 - 25 - 5) =
    ±220 mm, so the legs stand at -220, 68.3 and 220: the widest gap 288.3 mm.
    With two crossties it is 220 - 68.3 = 151.7 mm."""
    column = _column()
    column.set_ties(d_b=10 * mm, s=20 * cm, crossties_height=1)
    one = _y(column.shear_check_results([Forces(label="U", V_y=10 * kN)])).s_w
    column.set_ties(d_b=10 * mm, s=20 * cm, crossties_height=2)
    two = _y(column.shear_check_results([Forces(label="U", V_y=10 * kN)])).s_w
    bar = 205 / 3
    assert one.to("mm").magnitude == pytest.approx(220 + bar)
    assert two.to("mm").magnitude == pytest.approx(220 - bar)


def test_a_us_customary_column_reports_in_kip_and_inches() -> None:
    concrete = Concrete_ACI_318_19(name="C4000", f_c=4000 * psi)
    column = RectangularColumn(
        label="C",
        concrete=concrete,
        steel_bar=SteelBar(name="G60", f_y=60 * ksi),
        width=16 * inch,
        height=24 * inch,
        c_c=1.5 * inch,
    )
    column.set_longitudinal_rebar(d_b_corner=1 * inch, n_width=1, n_height=2)
    column.set_ties(d_b=0.5 * inch, s=6 * inch)
    check = _z(column.shear_check_results([Forces(label="U", N_x=100 * kip, V_z=60 * kip)]))
    assert check.V_capacity.units == (1 * kip).units
    assert check.d.to("inch").magnitude == pytest.approx(24 - 1.5 - 0.5 - 0.5)
    assert 0 < check.DCR < 1


# ---------------------------------------------------------------------------
# The element
# ---------------------------------------------------------------------------


def test_the_node_checks_both_directions_of_every_combination() -> None:
    column = _column()
    column.set_ties(d_b=10 * mm, s=15 * cm)
    forces = [Forces(label="A", N_x=100 * kN, V_y=40 * kN, V_z=80 * kN), Forces(label="B", N_x=50 * kN, V_z=20 * kN)]
    Node(column, forces).check()
    assert [(c.label, c.direction) for c in column.shear_checks] == [("A", "y"), ("A", "z"), ("B", "y"), ("B", "z")]
    assert column.shear_checks[0].V_demand.to("kN").magnitude == pytest.approx(40)


def test_results_are_dropped_with_the_ties_they_described() -> None:
    column = _column()
    column.set_ties(d_b=10 * mm, s=15 * cm)
    with pytest.raises(DesignNotRunError):
        column.shear_checks
    column.check([Forces(label="A", V_z=80 * kN)])
    assert column.shear_checks
    column.set_ties(d_b=8 * mm, s=15 * cm)
    with pytest.raises(DesignNotRunError):
        column.shear_checks
    assert column.warnings == ()


def test_what_the_column_does_not_do_says_so() -> None:
    column = _column()
    column.set_ties(d_b=10 * mm, s=15 * cm)
    force = [Forces(label="A", V_z=80 * kN)]
    for call in (column.design, column.design_shear, column.check_flexure, column.design_flexure):
        with pytest.raises(NotImplementedError, match="column-shear-roadmap"):
            call(force)
    en = _column(Concrete_EN_1992_2004(name="C25", f_c=25 * MPa))
    with pytest.raises(NotImplementedError, match="EN 1992-2004"):
        en.shear_check_results(force)


@pytest.mark.parametrize(
    "build, message",
    [
        (lambda c: c.set_ties(d_b=10 * mm, s=15 * cm, crossties_width=2), "crossties_width"),
        (lambda c: c.set_ties(d_b=0 * mm, s=0 * cm, crossties_height=1), "perimeter tie"),
        (lambda c: c.set_ties(d_b=10 * mm, s=0 * cm), "both"),
    ],
    ids=["more crossties than bars", "crossties without a tie", "a tie with no spacing"],
)
def test_impossible_ties_are_refused(build: Any, message: str) -> None:
    column = _column()
    with pytest.raises(ValueError, match=message):
        build(column)
        column.shear_check_results([Forces(label="A", V_z=80 * kN)])


def test_a_column_without_bars_is_refused() -> None:
    column = RectangularColumn(
        label="C",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=STEEL,
        width=30 * cm,
        height=30 * cm,
        c_c=25 * mm,
    )
    with pytest.raises(ValueError, match="longitudinal bars"):
        column.shear_check_results([Forces(label="A", V_z=80 * kN)])


# ---------------------------------------------------------------------------
# Warnings
# ---------------------------------------------------------------------------


def test_the_ties_a_shear_needs_are_warned_per_direction() -> None:
    """The 30x50 with a bare Ø6 tie every 30 cm under N = 0, V_y = 120 kN:
    reinforcement required, its area short, its spacing past d/2 and its legs
    past d across; V_z = 0 asks nothing. The moment is not checked."""
    column = _column()
    column.set_ties(d_b=6 * mm, s=30 * cm)
    column.check([Forces(label="U", V_y=120 * kN, M_y=40 * kNm)])
    codes = {(w.code, w.values.get("direction")) for w in column.warnings}
    assert codes == {
        ("column_ties_below_required", "y"),
        ("column_tie_spacing_exceeds_max", "y"),
        ("column_leg_spacing_exceeds_max", "y"),
        ("column_force_component_not_checked", None),
    }
    shortfall = next(w for w in column.warnings if w.code == "column_ties_below_required")
    assert shortfall.message.startswith("Ties for V_y")
    assert shortfall.combinations == ("U",)
    mento.set_language("es")
    assert next(w for w in column.warnings if w.code == "column_ties_below_required").message.startswith(
        "Estribos para V_y"
    )


def test_a_shear_past_the_section_is_warned() -> None:
    column = _column()
    column.set_ties(d_b=12 * mm, s=10 * cm)
    column.check([Forces(label="U", V_z=900 * kN)])
    limit = next(w for w in column.warnings if w.code == "column_shear_exceeds_section_limit")
    assert limit.values["direction"] == "z"
    assert limit.values["V"].to("kN").magnitude == pytest.approx(900)


def test_every_column_message_is_in_the_spanish_catalog() -> None:
    from mento.design_warnings import _MESSAGES

    assert [m for code, m in _MESSAGES.items() if code.startswith("column_") and m not in ES] == []
