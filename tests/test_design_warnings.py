"""Structured warnings: the detailing limits a section misses, as data."""

from types import MappingProxyType
from typing import Any, Iterator

import pytest

import mento
from mento import (
    Concrete_ACI_318_19,
    Concrete_EN_1992_2004,
    DesignWarning,
    Forces,
    Node,
    RectangularBeam,
    SteelBar,
)
from mento.codes.registry import design_code
from mento.units import MPa, Quantity, cm, kN, kNm, mm

FORCES = [
    Forces(label="1.2D+1.6L", V_z=250 * kN, M_y=150 * kNm),
    Forces(label="1.4D", V_z=120 * kN, M_y=-40 * kNm),
]


@pytest.fixture(autouse=True)
def english() -> Iterator[None]:
    mento.set_language("en")
    yield
    mento.set_language("en")


def _beam(width: Quantity = 20 * cm, height: Quantity = 60 * cm) -> RectangularBeam:
    return RectangularBeam(
        label="V101",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=width,
        height=height,
        c_c=25 * mm,
    )


def _poorly_detailed() -> tuple[RectangularBeam, Node]:
    beam = _beam()
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=10 * mm)
    beam.set_longitudinal_rebar_top(n1=6, d_b1=25 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=35 * cm)
    node = Node(section=beam, forces=FORCES)
    node.check()
    return beam, node


def _by_code(warnings: tuple[DesignWarning, ...]) -> dict[str, DesignWarning]:
    return {w.code: w for w in warnings}


def test_a_sound_design_has_no_warnings() -> None:
    beam = _beam()
    node = Node(section=beam, forces=FORCES)
    node.design()

    assert node.warnings == ()


def test_each_missed_limit_is_one_warning_with_a_stable_code() -> None:
    _, node = _poorly_detailed()
    found = _by_code(node.warnings)

    assert set(found) == {
        "As_below_min",
        "As_above_max",
        "bars_do_not_fit",
        "Av_below_min",
        "stirrup_spacing_exceeds_max",
    }
    assert found["As_below_min"].face == "bottom"
    assert found["As_above_max"].face == "top"
    assert found["bars_do_not_fit"].face == "top"
    # Only the positive moment asks the bottom face for steel.
    assert found["As_below_min"].combinations == ("1.2D+1.6L",)
    assert found["stirrup_spacing_exceeds_max"].combinations == ("1.2D+1.6L", "1.4D")


def test_values_are_quantities_in_report_units() -> None:
    beam, node = _poorly_detailed()
    found = _by_code(node.warnings)

    below = found["As_below_min"].values
    assert isinstance(below, MappingProxyType)
    assert below["A_s"].units == (1 * cm**2).units
    assert below["A_s"].magnitude == pytest.approx(beam._A_s_bot.to("cm**2").magnitude)
    assert below["A_s"] < below["A_s_min"]

    spacing = found["stirrup_spacing_exceeds_max"].values
    assert spacing["s"] == 35 * cm
    # The governing combination is the one furthest past the limit: the high
    # shear halves the spacing limit to d/4.
    assert spacing["s_max"].to("cm").magnitude == pytest.approx(beam._d_shear.to("cm").magnitude / 4)


def test_messages_follow_the_language() -> None:
    _, node = _poorly_detailed()
    english = _by_code(node.warnings)["As_below_min"].message
    mento.set_language("es")
    spanish = _by_code(node.warnings)["As_below_min"].message

    assert english.startswith("Steel on the bottom face: A_s = 1.57 cm²")
    assert spanish.startswith("Armadura en la cara inferior: A_s = 1.57 cm²")
    assert all("{" not in w.message for w in node.warnings)


def test_a_warning_does_not_change_the_dcr() -> None:
    """A strong but badly spaced cage: the check passes, the spacing does not."""
    beam = _beam()
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    beam.set_transverse_rebar(n_stirrups=2, d_b=12 * mm, s_l=40 * cm)
    node = Node(section=beam, forces=[Forces(label="V", V_z=120 * kN, M_y=100 * kNm)])
    node.check()

    assert "stirrup_spacing_exceeds_max" in _by_code(node.warnings)
    assert beam.shear_design.DCR < 1
    assert beam.shear_design.DCR == max(check.DCR for check in beam.shear_checks)


def test_no_stirrups_under_shear_asks_for_them() -> None:
    beam = _beam()
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam, forces=[Forces(label="V", V_z=120 * kN)])
    node.check_shear()

    found = _by_code(node.warnings)
    assert "stirrups_required" in found
    assert found["stirrups_required"].values["A_v_req"].magnitude > 0


def test_a_section_too_small_for_the_shear_says_so() -> None:
    beam = _beam(height=40 * cm)
    node = Node(section=beam, forces=[Forces(label="big", V_z=600 * kN, M_y=100 * kNm)])
    node.design()

    found = _by_code(node.warnings)
    assert "shear_exceeds_section_limit" in found
    values = found["shear_exceeds_section_limit"].values
    assert values["V"] > values["V_max"]
    # The design placed at least A_v,min, so the limit is the ØVmax the
    # report prints: rows (a)/(b) of Table 22.5.5.1 are already in it.
    assert values["V_max"].to("kN").magnitude == pytest.approx(beam._phi_V_max.to("kN").magnitude)


def _bare_aci_beam() -> RectangularBeam:
    """20x60, f'c 25, 2Ø12 below, 2Ø10 above, no stirrups: the review's e14."""
    beam = _beam()
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=12 * mm)
    beam.set_longitudinal_rebar_top(n1=2, d_b1=10 * mm)
    beam.set_transverse_rebar(0, 0 * mm, 0 * cm)
    return beam


def test_aci_section_limit_is_read_with_the_stirrups_the_code_requires() -> None:
    """20x60, f'c 25, 2Ø12 below, 2Ø10 above, c_c 25 mm, no stirrups.

    d = 600 - 25 - 6 = 569 mm, A_cv = 113 800 mm², ρw = 226/113 800 =
    0.00199, λs = sqrt(2/(1 + 0.004*569)) = 0.781. As the section is, row (c)
    of ACI 318-19 Table 22.5.5.1 gives V_c = 0.66*0.781*0.00199^(1/3)*sqrt(25)*
    A_cv = 36.9 kN, and Eq. (22.5.1.2) reads phi*(V_c + 0.66*sqrt(f'c)*bw*d) =
    0.75*(36.9 + 375.5) = 309 kN: under 320 kN bd94d2f said "enlarge the
    section", while 2eØ10/10 on the same section passes at DCR 0.92. Once
    the section carries the A_v,min that §9.6.3.1 requires of it anyway, rows
    (a)/(b) apply: V_c = max(0.17, 0.66*0.00199^(1/3))*sqrt(25)*A_cv = 0.85 MPa*
    A_cv = 96.7 kN and the limit is 0.75*(96.7 + 375.5) = 354 kN. 320 kN is
    short of it -- stirrups are what is missing -- and 360 kN is past it
    whatever the stirrups.
    """
    beam = _bare_aci_beam()
    node = Node(section=beam, forces=[Forces(label="U", V_z=320 * kN, M_y=20 * kNm)])
    node.check_shear()
    found = _by_code(node.warnings)
    assert "stirrups_required" in found
    assert "shear_exceeds_section_limit" not in found

    beam = _bare_aci_beam()
    node = Node(section=beam, forces=[Forces(label="U", V_z=360 * kN, M_y=20 * kNm)])
    node.check_shear()
    limit = _by_code(node.warnings)["shear_exceeds_section_limit"]
    assert limit.values["V"].to("kN").magnitude == pytest.approx(360)
    assert limit.values["V_max"].to("kN").magnitude == pytest.approx(354.2, abs=0.5)
    assert limit.combinations == ("U",)


def _bare_en_beam(stirrups: bool = False) -> RectangularBeam:
    """EN C25/30, 30x50, 3Ø16 below, 2Ø12 above, c_c 25 mm: the review's e15."""
    beam = RectangularBeam(
        label="V",
        concrete=Concrete_EN_1992_2004(name="C25/30", f_c=25 * MPa),
        steel_bar=SteelBar(name="B500S", f_y=500 * MPa),
        width=30 * cm,
        height=50 * cm,
        c_c=25 * mm,
    )
    beam.set_longitudinal_rebar_bot(n1=3, d_b1=16 * mm)
    beam.set_longitudinal_rebar_top(n1=2, d_b1=12 * mm)
    if stirrups:
        beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=15 * cm)
    else:
        beam.set_transverse_rebar(0, 0 * mm, 0 * cm)
    return beam


def test_en_beam_without_stirrups_is_asked_for_what_the_demand_needs() -> None:
    """EN C25/30, 30x50, 3Ø16 below, 2Ø12 above, c_c 25 mm, no stirrups.

    d = 500 - 25 - 8 = 467 mm, z = 0.9*d = 420.3 mm, f_ywd = 500/1.15 =
    434.8 MPa. V_Rd,c of Eq. (6.2.a) = 0.12*1.654*(100*0.00431*25)^(1/3)*
    300*467 = 61.4 kN. Under 150 kN EN 1992-1-1 §6.2.1(5) asks for enough
    stirrups that V_Ed ≤ V_Rd: with cot θ = 2.5 (150 kN is far under
    V_Rd,max at 21.8°, 391 kN) that is A_sw/s = V_Ed/(z*f_ywd*cot θ) =
    150 000/(420.3*434.8*2.5) = 0.328 mm²/mm = 3.28 cm²/m, not the
    2.40 cm²/m of Eq. (9.5N) (0.08*sqrt(25)/500*300) bd94d2f quoted whatever
    the shear. Under 30 kN ≤ V_Rd,c §6.2.1(3) needs none calculated and (4)
    the minimum, which is what a beam is asked for.
    """
    beam = _bare_en_beam()
    node = Node(section=beam, forces=[Forces(label="ULS", V_z=150 * kN, M_y=50 * kNm)])
    node.check_shear()
    required = _by_code(node.warnings)["stirrups_required"]
    assert required.values["A_v_req"].to("cm**2/m").magnitude == pytest.approx(3.28, abs=0.01)
    assert required.values["A_v_min"].to("cm**2/m").magnitude == pytest.approx(2.40, abs=0.01)
    assert beam.shear_checks[0].A_v_req.to("cm**2/m").magnitude == pytest.approx(3.28, abs=0.01)
    # The capacity of the section as it is stays the concrete's alone.
    assert beam.shear_checks[0].V_capacity.to("kN").magnitude == pytest.approx(61.4, abs=0.1)

    beam = _bare_en_beam()
    node = Node(section=beam, forces=[Forces(label="ULS", V_z=30 * kN, M_y=50 * kNm)])
    node.check_shear()
    assert _by_code(node.warnings)["stirrups_required"].values["A_v_req"].to("cm**2/m").magnitude == pytest.approx(
        2.40, abs=0.01
    )


def test_en_section_limit_is_the_strut_at_45_degrees_with_or_without_stirrups() -> None:
    """Same section. V_Rd,max of EN 1992-1-1 §6.2.3(3), Eq. (6.9), is largest
    at θ = 45°, where cot θ + tan θ = 2: αcw*bw*z*ν1*fcd/2 = 1*300*420.3*
    0.54*16.67/2 = 567 kN without stirrups (ν1 = 0.6*(1 - 25/250), fcd =
    25/1.5) and 300*413.1*0.54*16.67/2 = 558 kN with 1eØ8/15, whose legs
    lower d to 459 mm. 150 kN is under both -- 1eØ8/15 carries it at DCR
    0.50 -- so the bare section is not told to grow, as bd94d2f did by
    quoting V_Rd,c = 61.4 kN as the most it can carry. 600 kN is past both,
    and no stirrup changes that.
    """
    beam = _bare_en_beam()
    node = Node(section=beam, forces=[Forces(label="ULS", V_z=150 * kN, M_y=50 * kNm)])
    node.check_shear()
    assert "shear_exceeds_section_limit" not in _by_code(node.warnings)

    beam = _bare_en_beam()
    node = Node(section=beam, forces=[Forces(label="ULS", V_z=600 * kN, M_y=50 * kNm)])
    node.check_shear()
    limit = _by_code(node.warnings)["shear_exceeds_section_limit"]
    assert limit.values["V_max"].to("kN").magnitude == pytest.approx(567.4, abs=0.5)

    beam = _bare_en_beam(stirrups=True)
    node = Node(section=beam, forces=[Forces(label="ULS", V_z=600 * kN, M_y=50 * kNm)])
    node.check_shear()
    limit = _by_code(node.warnings)["shear_exceeds_section_limit"]
    assert limit.values["V_max"].to("kN").magnitude == pytest.approx(557.7, abs=0.5)
    assert limit.values["V_max"].to("kN").magnitude == pytest.approx(beam._V_Rd_max.to("kN").magnitude)


def test_values_only_checks_record_warnings_too() -> None:
    beam = _beam()
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=10 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=10 * mm, s_l=35 * cm)
    beam.flexure_check_results(FORCES)
    beam.shear_check_results(FORCES)

    assert {"As_below_min", "stirrup_spacing_exceeds_max"} <= set(_by_code(beam.warnings))


def test_en_beam_reports_its_own_limits() -> None:
    beam = RectangularBeam(
        label="V",
        concrete=Concrete_EN_1992_2004(name="C25", f_c=25 * MPa),
        steel_bar=SteelBar(name="B500S", f_y=500 * MPa),
        width=20 * cm,
        height=60 * cm,
        c_c=25 * mm,
    )
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=60 * cm)
    node = Node(section=beam, forces=[Forces(label="V", V_z=100 * kN, M_y=80 * kNm)])
    node.check()

    codes = set(_by_code(node.warnings))
    assert "stirrup_spacing_exceeds_max" in codes
    assert "stirrup_diameter_below_min" not in codes


@pytest.mark.parametrize("concrete", [Concrete_ACI_318_19, mento.Concrete_CIRSOC_201_25], ids=["ACI", "CIRSOC"])
def test_a_thin_stirrup_placed_for_shear_alone_is_not_below_any_code_minimum(concrete: type) -> None:
    """20x50, 3Ø16 below, 2Ø10 above, 1eØ8/15, M = 60 kNm, V = 80 kN: singly
    reinforced (DCR 0.61 in flexure, 0.52 in shear), so no bar is compression
    steel that a stirrup has to support. ACI 318-19 §9.7.6.4.2 / CIRSOC 201-25
    §9.7.6.4.2, Tabla 9.7.6.4.2 size only the stirrups of §9.7.6.4.1, those
    laterally supporting compression reinforcement; neither code states a
    minimum for a stirrup placed for shear. The 10 mm bd94d2f quoted as "the
    minimum" is the bottom of mento's ACI catalogue -- a design preference,
    which stays one -- and the 6 mm of the CIRSOC hook the bottom of its
    catalogue. The Ø6 of ``_poorly_detailed`` is not a code limit either."""
    beam = RectangularBeam(
        label="V",
        concrete=concrete(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=20 * cm,
        height=50 * cm,
        c_c=25 * mm,
    )
    beam.set_longitudinal_rebar_bot(n1=3, d_b1=16 * mm)
    beam.set_longitudinal_rebar_top(n1=2, d_b1=10 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=15 * cm)
    node = Node(section=beam, forces=[Forces(label="U", V_z=80 * kN, M_y=60 * kNm)])
    node.check()

    assert not beam._doubly_reinforced
    assert beam.flexure_design.bottom.DCR < 1 and beam.shear_design.DCR < 1
    assert "stirrup_diameter_below_min" not in _by_code(node.warnings)
    assert "stirrup_diameter_below_min" not in mento.design_warnings._MESSAGES


def test_warnings_are_empty_before_any_check() -> None:
    assert _beam().warnings == ()


def test_slab_bar_spacing_limits_are_warned() -> None:
    from mento import OneWaySlab
    from mento.units import m

    slab = OneWaySlab(
        label="L1",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=1 * m,
        height=20 * cm,
        c_c=25 * mm,
    )
    slab.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=60 * cm)
    slab.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=2 * cm)
    Node(section=slab, forces=[Forces(label="M", M_y=10 * kNm)]).check_flexure()

    found = {(w.code, w.face) for w in slab.warnings}
    assert ("bar_spacing_exceeds_max", "bottom") in found
    assert ("bar_spacing_below_min", "top") in found


def test_wall_mesh_spacing_names_the_limit_as_mentos_own() -> None:
    """ACI 318-19 wall 25x150, hw = 3 m, Ø16/40 both ways, Vu = 100 kN (DCR 0.14).

    §11.7.3.1 caps the horizontal spacing at the lesser of 3h = 750 mm and
    450 mm, and adds lw/5 = 300 mm only "if shear reinforcement is required
    for in-plane strength", which at DCR 0.14 it is not; §11.7.2.1 does the
    same with lw/3 = 500 mm. mento takes lw/5 and lw/3 always, so 400 mm is
    past the 300 mm it applies horizontally and within the 450 mm vertically.
    The message says whose limit it is, instead of "exceeds the maximum", which
    read as the clause's.
    """
    from mento import ShearWall
    from mento.units import m

    wall = ShearWall(
        label="W",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        thickness=25 * cm,
        length=1.5 * m,
        height=3.0 * m,
        c_c=20 * mm,
    )
    wall.set_horizontal_rebar(d_b=16 * mm, s=40 * cm)
    wall.set_vertical_rebar(d_b=16 * mm, s=40 * cm)
    wall.shear_check_results([Forces(label="E", V_z=100 * kN)])

    spacing = [w for w in wall.warnings if w.code == "mesh_spacing_exceeds_max"]
    assert len(spacing) == 1
    assert spacing[0].values["s"].to("mm").magnitude == pytest.approx(400)
    assert spacing[0].values["s_max"].to("mm").magnitude == pytest.approx(300)
    assert spacing[0].message.startswith("Horizontal wall mesh spacing: 40 cm exceeds the limit mento applies, 30 cm")
    assert "lw/5" in spacing[0].message
    mento.set_language("es")
    assert "el límite que aplica mento" in wall.warnings[0].message


def test_stirrup_legs_too_far_apart_across_the_width_are_worded_as_such() -> None:
    """ACI 318-19 80x40 H25, 4Ø16 at the bottom, 1eØ8/15, Vu = 250 kN.

    One closed stirrup has two legs, and across 80 cm they sit
    800 - 2*25 - 8 = 742 mm apart, centre to centre. With d = 400 - 25 - 8 -
    8 = 359 mm, phi*Vc = 0.75*0.17*sqrt(25)*800*359 = 183.1 kN and
    Vs,req = (250 - 183.1)/0.75 = 89.2 kN, under the 0.33*sqrt(25)*800*359 =
    473.9 kN of Table 9.7.6.2.2: the legs may be up to d = 359 mm apart
    across the width. The along-the-member spacing of 15 cm is within d/2.
    The width direction carries a template of its own, in both languages.
    """
    beam = _beam(width=80 * cm, height=40 * cm)
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=15 * cm)
    node = Node(section=beam, forces=[Forces(label="V", V_z=250 * kN)])
    node.check()

    spacing = [w for w in node.warnings if w.code == "stirrup_spacing_exceeds_max"]
    assert len(spacing) == 1
    assert spacing[0].values["s"].to("mm").magnitude == pytest.approx(742)
    assert spacing[0].values["s_max"].to("mm").magnitude == pytest.approx(359)
    assert spacing[0].message == "Stirrup leg spacing across the width: 74.2 cm exceeds the maximum 35.9 cm."
    mento.set_language("es")
    assert _by_code(node.warnings)["stirrup_spacing_exceeds_max"].message == (
        "Separación de las ramas de estribo en el ancho: 74.2 cm supera la máxima 35.9 cm."
    )


def test_the_other_direction_of_each_wall_mesh_limit_is_worded_as_such() -> None:
    """ACI 318-19 wall 25x150, hw = 3 m, Ø10/20 horizontal, Ø16/60 vertical, Vu = 900 kN.

    Horizontal: 2 curtains of Ø10 at 200 mm give rho_t = 2*78.54/(250*200) =
    0.00314. hw/lw = 2, so alpha_c = 0.17 and phi*Vn >= Vu asks for
    rho_t = (900/0.75 - 0.17*sqrt(25)*250*1500)/(250*1500*420) = (1200 -
    318.75)/157.5 = 0.0056 (§11.5.4.3). Vertical: §11.7.2.1 with lw/3 taken
    always caps the spacing at min(3*250, 450, 1500/3) = 450 mm, and the bars
    sit at 600. Each is the direction the other wall test does not word.
    """
    from mento import ShearWall
    from mento.units import m

    wall = ShearWall(
        label="W",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        thickness=25 * cm,
        length=1.5 * m,
        height=3.0 * m,
        c_c=20 * mm,
    )
    wall.set_horizontal_rebar(d_b=10 * mm, s=20 * cm)
    wall.set_vertical_rebar(d_b=16 * mm, s=60 * cm)
    wall.shear_check_results([Forces(label="E", V_z=900 * kN)])

    found = _by_code(wall.warnings)
    assert (
        found["mesh_ratio_below_min"].message == "Horizontal wall mesh: ρt = 0.00314 is below the required ρt = 0.0056."
    )
    assert found["mesh_spacing_exceeds_max"].message == (
        "Vertical wall mesh spacing: 60 cm exceeds the limit mento applies, 45 cm "
        "(§11.7.2.1 with lw/3 taken always: conservative)."
    )
    mento.set_language("es")
    found = _by_code(wall.warnings)
    assert found["mesh_ratio_below_min"].message == (
        "Malla horizontal del muro: ρt = 0.00314 es menor que la requerida ρt = 0.0056."
    )
    assert found["mesh_spacing_exceeds_max"].message == (
        "Separación de la malla vertical del muro: 60 cm supera el límite que aplica mento, 45 cm "
        "(§11.7.2.1 con lw/3 siempre: conservador)."
    )


def test_an_unlabelled_combination_is_named_by_its_position() -> None:
    """The same poorly detailed beam under the same two forces without labels.

    ``Forces.label`` defaults to ``None``, and bd94d2f dropped it, so every
    per-combination warning came out with ``combinations == ()`` -- the value
    the docstring reserves for a limit of the section alone. Now the first
    force is ``#1`` and the second ``#2``, and only the spacing warning, which
    is the section's, stays empty.
    """
    beam = _beam()
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=10 * mm)
    beam.set_longitudinal_rebar_top(n1=6, d_b1=25 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=35 * cm)
    forces = [Forces(V_z=250 * kN, M_y=150 * kNm), Forces(V_z=120 * kN, M_y=-40 * kNm)]
    node = Node(section=beam, forces=forces)
    node.check()

    found = _by_code(node.warnings)
    assert found["As_below_min"].combinations == ("#1",)
    assert found["stirrup_spacing_exceeds_max"].combinations == ("#1", "#2")
    assert found["bars_do_not_fit"].combinations == ()

    # The values-only entry points name them the same way.
    beam.flexure_check_results(forces)
    beam.shear_check_results(forces)
    found = _by_code(beam.warnings)
    assert found["As_below_min"].combinations == ("#1",)
    assert found["stirrup_spacing_exceeds_max"].combinations == ("#1", "#2")

    # A label given is kept, and a missing one still counts its position.
    node = Node(
        section=beam, forces=[Forces(label="D", V_z=250 * kN, M_y=150 * kNm), Forces(V_z=120 * kN, M_y=-40 * kNm)]
    )
    node.check()
    assert _by_code(node.warnings)["stirrup_spacing_exceeds_max"].combinations == ("D", "#2")


def test_an_unlabelled_wall_combination_is_named_by_its_position() -> None:
    from mento import ShearWall
    from mento.units import m

    wall = ShearWall(
        label="W",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        thickness=25 * cm,
        length=1.5 * m,
        height=3.0 * m,
        c_c=20 * mm,
    )
    wall.set_horizontal_rebar(d_b=6 * mm, s=40 * cm)
    wall.set_vertical_rebar(d_b=6 * mm, s=40 * cm)
    wall.shear_check_results([Forces(V_z=100 * kN), Forces(V_z=3000 * kN)])

    found = _by_code(wall.warnings)
    assert found["shear_exceeds_section_limit"].combinations == ("#2",)
    assert found["mesh_ratio_below_min"].combinations[0] in ("#1", "#2")


def test_a_value_just_past_its_limit_prints_with_enough_digits_to_show_it() -> None:
    """Three significant figures print 13.00 cm and 12.95 cm both as "13 cm",
    and bd94d2f read "Stirrup spacing along the member: 13 cm exceeds the
    maximum 13 cm." A message quotes as many digits as it takes for two
    values that differ to print differently, and no more than three where
    they already do."""
    from mento.design_warnings import _Raw, collect

    tight = collect([_Raw("stirrup_spacing_exceeds_max", {"s": 13.0 * cm, "s_max": 12.96 * cm, "direction": "l"})])
    assert tight[0].message == "Stirrup spacing along the member: 13 cm exceeds the maximum 12.96 cm."
    closer = collect([_Raw("stirrup_spacing_exceeds_max", {"s": 13.001 * cm, "s_max": 13.0 * cm, "direction": "l"})])
    assert closer[0].message == "Stirrup spacing along the member: 13.001 cm exceeds the maximum 13 cm."
    loose = collect([_Raw("stirrup_spacing_exceeds_max", {"s": 35.0 * cm, "s_max": 13.912 * cm, "direction": "l"})])
    assert loose[0].message == "Stirrup spacing along the member: 35 cm exceeds the maximum 13.9 cm."
    ratio = collect([_Raw("mesh_ratio_below_min", {"direction": "v", "rho": 0.0025, "rho_min": 0.00251})])
    assert ratio[0].message.endswith("ρl = 0.0025 is below the minimum ρl,min = 0.00251.")


def test_the_stirrup_spacing_of_the_review_prints_its_limit() -> None:
    """10x30 designed for 80 kNm gets 1eØ10/13 with d = 300 - 25 - 10 - 4 =
    261 mm (2Ø8 left below), d/2 = 13.05 cm; 1Ø12 set below by hand lowers
    d to 259 mm and d/2 to 12.95 cm, which 13 cm now exceeds -- by 0.05 cm,
    which the message has to show."""
    beam = _beam(width=10 * cm, height=30 * cm)
    node = Node(section=beam, forces=[Forces(label="U", V_z=20 * kN, M_y=80 * kNm)])
    node.design()
    beam.set_longitudinal_rebar_bot(n1=1, d_b1=12 * mm)
    node.check()

    spacing = _by_code(node.warnings)["stirrup_spacing_exceeds_max"]
    assert spacing.values["s"].to("cm").magnitude == pytest.approx(13.0)
    assert spacing.values["s_max"].to("cm").magnitude == pytest.approx(12.95)
    assert spacing.message == "Stirrup spacing along the member: 13 cm exceeds the maximum 12.95 cm."


def test_a_warning_prints_as_its_message() -> None:
    _, node = _poorly_detailed()
    warning = node.warnings[0]

    assert str(warning) == warning.message


def test_bars_that_do_not_fit_are_warned_after_a_design() -> None:
    """A 10x20 web cannot take what 200 kNm asks for: the search finds no layout at all."""
    beam = _beam(width=10 * cm, height=20 * cm)
    node = Node(section=beam, forces=[Forces(label="M", M_y=200 * kNm, V_z=20 * kN)])
    node.design()

    # The moment makes it doubly reinforced, so neither the tension steel below
    # nor the compression steel above has a layout that fits.
    not_fitting = [w for w in node.warnings if w.code == "bars_do_not_fit"]
    assert {w.face for w in not_fitting} == {"bottom", "top"}
    assert all(w.values == {} for w in not_fitting)


def test_a_design_short_of_the_moment_is_warned_until_the_bars_reach_it() -> None:
    """A 12x30 web takes 2Ø16 + 2Ø16 = 8.04 cm² at most beside the 1eØ6 it ends with; 60 kNm asks for 8.14 below.

    The moment needs compression steel too, 7.27 cm² above, where 2Ø12 + 2Ø12
    = 4.52 fit: both faces are short, and each is warned for as long as it
    carries what the design left. (At 40 kNm this web used to be declared
    short because Ø16 did not fit beside the 8 mm starter stirrup; a full
    design now redoes the flexure with the stirrup the shear design chose,
    and 2Ø16 + 2Ø12 carry it.)
    """
    beam = RectangularBeam(
        label="101",
        concrete=mento.Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=12 * cm,
        height=30 * cm,
        c_c=25 * mm,
    )
    node = Node(section=beam, forces=[Forces(label="1.4D", V_z=50 * kN, M_y=60 * kNm)])
    node.design()

    short = {w.face: w for w in node.warnings if w.code == "As_below_required"}
    assert set(short) == {"bottom", "top"}
    assert short["bottom"].values["A_s"].to("cm**2").magnitude == pytest.approx(8.04, rel=1e-3)
    assert short["bottom"].values["A_s_req"].to("cm**2").magnitude == pytest.approx(8.14, rel=1e-2)
    mento.set_language("es")
    assert "agrandar la sección" in _by_code(node.warnings)["As_below_required"].message

    # Bars set by hand that reach the area clear it, face by face.
    beam.set_longitudinal_rebar_bot(2, 20 * mm, 0, None, 2, 20 * mm)
    assert {w.face for w in node.warnings if w.code == "As_below_required"} == {"top"}
    beam.set_longitudinal_rebar_top(2, 20 * mm, 0, None, 2, 20 * mm)
    assert "As_below_required" not in {w.code for w in node.warnings}


def test_the_maximum_is_only_read_on_the_face_in_tension() -> None:
    """29.45 cm² on the bottom is past A_s,max, but only a positive moment pulls it."""

    def heavy_bottom() -> RectangularBeam:
        beam = _beam()
        beam.set_transverse_rebar(n_stirrups=1, d_b=10 * mm, s_l=20 * cm)
        beam.set_longitudinal_rebar_bot(2, 25 * mm, 1, 25 * mm, 2, 25 * mm, 1, 25 * mm)
        beam.set_longitudinal_rebar_top(2, 12 * mm)
        return beam

    # A negative moment makes those bars compression steel, and no moment
    # pulls neither face: no warning, and the report row passes.
    beam = heavy_bottom()
    node = Node(section=beam, forces=[Forces(label="neg", M_y=-50 * kNm), Forces(label="zero", M_y=0 * kNm)])
    node.check_flexure()
    assert "As_above_max" not in _by_code(node.warnings)
    assert beam._data_min_max_flexure["Ok?"][2] == "✅"

    beam = heavy_bottom()
    node = Node(section=beam, forces=[Forces(label="pos", M_y=150 * kNm)])
    node.check_flexure()
    over = _by_code(node.warnings)["As_above_max"]
    assert (over.face, over.combinations) == ("bottom", ("pos",))
    assert beam._data_min_max_flexure["Ok?"][2] == "❌"


def test_en_holds_both_faces_to_its_maximum() -> None:
    """EN 1992-1-1 §9.2.1.1(3) caps tension OR compression steel, so the face a
    moment compresses is read against its 4 % too -- unlike ACI 318-19, whose
    maximum is the ductility limit of the face in tension."""
    beam = RectangularBeam(
        label="E",
        concrete=Concrete_EN_1992_2004(name="C25", f_c=25 * MPa),
        steel_bar=SteelBar(name="B500S", f_y=500 * MPa),
        width=20 * cm,
        height=40 * cm,
        c_c=25 * mm,
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=20 * cm)
    beam.set_longitudinal_rebar_bot(2, 32 * mm, 1, 32 * mm, 2, 32 * mm, 1, 32 * mm)  # 48.3 cm² > 4 %
    beam.set_longitudinal_rebar_top(2, 16 * mm)
    node = Node(section=beam, forces=[Forces(label="neg", M_y=-40 * kNm)])
    node.check_flexure()
    over = [w for w in node.warnings if w.code == "As_above_max"]
    assert [w.face for w in over] == ["bottom"]
    assert over[0].combinations == ("neg",)


def test_en_reads_its_maximum_on_the_gross_section() -> None:
    """EN 1992-1-1 §9.2.1.1(3): A_s,max = 0.04 Ac, and Ac is the concrete section, b*h (§1.6).

    20x60 C30/37 B500S, 1eØ8/15, M_Ed = 200 kNm. 0.04*200*600 = 4800 mm² =
    48.0 cm², the same on both faces. 3Ø32 + 2Ø32 + 1Ø20 = 43.35 cm² below,
    with d = 526.2 mm, were read against 0.04*200*526.2 = 42.10 cm² -- b*d,
    10 % tighter than the clause -- and warned ``As_above_max``; the
    design's ``flexure_admissible`` gate, which drops alternatives and
    decides when no layout is left, said no as well. 3Ø32 + 3Ø32 = 48.25
    cm² is past the 48.0 and still warned.
    """

    def beam_with(*bottom: Any) -> RectangularBeam:
        beam = RectangularBeam(
            label="E",
            concrete=Concrete_EN_1992_2004(name="C30", f_c=30 * MPa),
            steel_bar=SteelBar(name="B500S", f_y=500 * MPa),
            width=20 * cm,
            height=60 * cm,
            c_c=25 * mm,
        )
        beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=15 * cm)
        beam.set_longitudinal_rebar_bot(*bottom)
        beam.set_longitudinal_rebar_top(2, 12 * mm)
        return beam

    admissible = design_code(Concrete_EN_1992_2004(name="C30", f_c=30 * MPa)).flexure_admissible
    assert admissible is not None
    forces = [Forces(label="ELU", M_y=200 * kNm)]

    beam = beam_with(3, 32 * mm, 0, None, 2, 32 * mm, 1, 20 * mm)
    node = Node(section=beam, forces=forces)
    node.check_flexure()
    assert beam.flexure_design.bottom.A_s.to("cm**2").magnitude == pytest.approx(43.35, abs=0.005)
    assert beam._d_bot.to("mm").magnitude == pytest.approx(526.2, abs=0.05)
    assert beam.flexure_design.bottom.A_s_max == beam.flexure_design.top.A_s_max
    assert beam.flexure_design.bottom.A_s_max.to("cm**2").magnitude == pytest.approx(48.0)
    assert "As_above_max" not in _by_code(node.warnings)
    assert admissible(beam, "bot")

    beam = beam_with(3, 32 * mm, 0, None, 3, 32 * mm)
    node = Node(section=beam, forces=forces)
    node.check_flexure()
    over = _by_code(node.warnings)["As_above_max"]
    assert over.values["A_s"].to("cm**2").magnitude == pytest.approx(48.25, abs=0.005)
    assert over.values["A_s_max"].to("cm**2").magnitude == pytest.approx(48.0)
    assert not admissible(beam, "bot")


# ---------------------------------------------------------------------------
# The minimum a face has to meet: ACI 318-19 / CIRSOC 201-25 §9.6.1.3
# ---------------------------------------------------------------------------


def _cirsoc_section(kind: str) -> RectangularBeam:
    """The two sections the relief was first seen short on, CIRSOC 201-25."""
    concrete = mento.Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    if kind == "slab":
        return mento.OneWaySlab(
            label="L1", concrete=concrete, steel_bar=steel, width=100 * cm, height=20 * cm, c_c=25 * mm
        )
    return RectangularBeam(label="V1", concrete=concrete, steel_bar=steel, width=20 * cm, height=60 * cm, c_c=25 * mm)


@pytest.mark.parametrize("kind, M_u", [("beam", 45), ("slab", 22)])
@pytest.mark.parametrize("sign", [1, -1], ids=["bottom", "top"])
def test_designed_face_never_below_the_minimum_it_has_to_meet(kind: str, M_u: float, sign: int) -> None:
    section = _cirsoc_section(kind)
    node = Node(section=section, forces=[Forces(label="U", M_y=sign * M_u * kNm)])
    node.design()

    face = section.flexure_design.bottom if sign > 0 else section.flexure_design.top
    if kind == "beam":
        # Light enough that the 4/3 relief governs: the face sits below the
        # A_s,min of §9.6.1.2 and complies all the same.
        assert face.A_s < face.A_s_min
        assert face.A_s_min_eff.to("cm**2").magnitude == pytest.approx((4 * face.A_s_calc / 3).to("cm**2").magnitude)
    else:
        # A slab takes 0.0018*b*h of §7.6.1.1, which §9.6.1.3 does not relieve.
        assert face.A_s_min.to("cm**2").magnitude == pytest.approx(0.0018 * 100 * 20)
        assert face.A_s_min_eff == face.A_s_min
        assert face.A_s >= face.A_s_min
    for f in (section.flexure_design.bottom, section.flexure_design.top):
        assert f.A_s_min_eff <= f.A_s_min
        assert f.A_s.to("cm**2").magnitude >= f.A_s_min_eff.to("cm**2").magnitude - 1e-9
    assert "As_below_min" not in [w.code for w in node.warnings]


def test_relieved_minimum_reads_the_detailed_report_as_passing() -> None:
    beam = _cirsoc_section("beam")
    node = Node(section=beam, forces=[Forces(label="U", M_y=-45 * kNm)])
    node.design()
    node.flexure_results_detailed()
    assert beam._data_min_max_flexure["Ok?"][0] == "✅ 9.6.1.3"


def test_layout_short_of_four_thirds_still_warns() -> None:
    # The requirement adopts 4/3 A_s_calc, so the old flag was set; the bars
    # checked by hand carry more than A_s_calc but less than 4/3 of it, which
    # the relief does not cover.
    beam = _cirsoc_section("beam")
    node = Node(section=beam, forces=[Forces(label="U", M_y=-45 * kNm)])
    beam.set_longitudinal_rebar_top(n1=2, d_b1=12 * mm)  # 2.26 cm²: A_s_calc < 2.26 < 4/3 A_s_calc
    node.check_flexure()
    top = beam.flexure_design.top
    assert top.A_s_calc < top.A_s < top.A_s_min_eff
    warning = next(w for w in node.warnings if w.code == "As_below_min")
    assert warning.values["A_s_min_eff"].to("cm**2").magnitude == pytest.approx(top.A_s_min_eff.to("cm**2").magnitude)


def test_the_warning_names_the_clause_minimum_and_the_one_the_face_has_to_meet() -> None:
    """CIRSOC 201-25 20x60, Mu = -45 kNm, top 2Ø12 = 2.26 cm², 1eØ8.

    d = 600 - 25 - 8 - 6 = 561 mm. §9.6.1.2: A_s,min = max(0.25*sqrt(25)/420,
    1.4/420)*200*561 = 0.003333*112 200 = 3.74 cm²; A_s,calc = 2.16 cm² and
    §9.6.1.3 relieves the face to 4/3 of it, A_s,min,eff = 2.88 cm². The
    same name cannot carry both numbers: ``flexure_design.top.A_s_min`` and
    the report's As,min column carry 3.74, so ``values["A_s_min"]`` does
    too, and the one the face is short of -- 2.88, which bd94d2f filed under
    ``A_s_min`` -- is ``A_s_min_eff``, which the message quotes.
    """
    beam = _cirsoc_section("beam")
    beam.set_longitudinal_rebar_top(n1=2, d_b1=12 * mm)
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=12 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=20 * cm)
    node = Node(section=beam, forces=[Forces(label="U", V_z=40 * kN, M_y=-45 * kNm)])
    node.check()

    warning = _by_code(node.warnings)["As_below_min"]
    assert warning.face == "top"
    assert warning.values["A_s"].to("cm**2").magnitude == pytest.approx(2.26, abs=0.01)
    assert warning.values["A_s_min"].to("cm**2").magnitude == pytest.approx(3.74, abs=0.01)
    assert warning.values["A_s_min_eff"].to("cm**2").magnitude == pytest.approx(2.88, abs=0.01)
    assert warning.values["A_s_min"] == beam.flexure_design.top.A_s_min
    assert "A_s,min,eff = 2.88 cm²" in warning.message
    assert "3.74" not in warning.message
    mento.set_language("es")
    assert "A_s,mín,ef = 2.88 cm²" in _by_code(node.warnings)["As_below_min"].message


def test_en_minimum_is_not_relieved() -> None:
    beam = RectangularBeam(
        label="V1",
        concrete=Concrete_EN_1992_2004(name="C25", f_c=25 * MPa),
        steel_bar=SteelBar(name="B500S", f_y=500 * MPa),
        width=20 * cm,
        height=60 * cm,
        c_c=25 * mm,
    )
    node = Node(section=beam, forces=[Forces(label="U", M_y=-45 * kNm)])
    node.design()
    top = beam.flexure_design.top
    assert top.A_s_min_eff == top.A_s_min


# ---------------------------------------------------------------------------
# The spacing is read off the section as it is now
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bars_first", [True, False], ids=["bars then stirrups", "stirrups then bars"])
def test_clear_spacing_follows_the_stirrup_whatever_the_call_order(bars_first: bool) -> None:
    """20x50, 4Ø12 below, 1eØ16, c_c 25 mm.

    The legs sit between the cover and the bars, so the clear space left for
    the four bars is (200 - 2*(25 + 16) - 4*12)/3 = 23.3 mm, under the 25 mm of
    ACI 318-19 §25.2.1 the settings ask for. With the Ø8 the settings assume
    it is (200 - 2*(25 + 8) - 48)/3 = 28.7 mm and passes. The warning is
    read off the section, so it cannot depend on whether the stirrups were
    set before or after the bars, nor wait for a reporting check to refresh
    it: on bd94d2f the "bars then stirrups" order kept the 28.7 mm and stayed
    silent until ``check_flexure`` ran.
    """
    beam = _beam(height=50 * cm)
    if bars_first:
        beam.set_longitudinal_rebar_bot(n1=4, d_b1=12 * mm)
        beam.set_transverse_rebar(n_stirrups=1, d_b=16 * mm, s_l=20 * cm)
    else:
        beam.set_transverse_rebar(n_stirrups=1, d_b=16 * mm, s_l=20 * cm)
        beam.set_longitudinal_rebar_bot(n1=4, d_b1=12 * mm)

    spacing = _by_code(beam.warnings)["clear_spacing_below_min"]
    assert spacing.face == "bottom"
    assert spacing.values["s"].to("mm").magnitude == pytest.approx(23.33, abs=0.01)
    assert spacing.values["s_min"].to("mm").magnitude == pytest.approx(25.0)

    # A lighter stirrup set afterwards widens the space again, and the
    # warning goes with it -- with no check in between.
    beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=20 * cm)
    assert "clear_spacing_below_min" not in _by_code(beam.warnings)
    assert beam._available_s_bot.to("mm").magnitude == pytest.approx(28.67, abs=0.01)


def test_bars_set_by_hand_that_fit_clear_the_flag_of_their_face() -> None:
    """A 10x30 web finds no layout for 80 kNm on either face (it is doubly
    reinforced and neither the tension nor the compression steel fits), so
    the design flags both. The flag is the search's verdict on the width,
    not a property of the bars: one Ø12 set by hand has 100 - 2*(25 + 10) -
    12 = 18 mm beside it and fits. On bd94d2f the flag outlived the bars and
    the face was reported as not fitting with a single bar on it; the other
    face keeps its flag until it is set too.
    """
    beam = _beam(width=10 * cm, height=30 * cm)
    node = Node(section=beam, forces=[Forces(label="U", V_z=20 * kN, M_y=80 * kNm)])
    node.design()

    def not_fitting() -> set[str | None]:
        return {w.face for w in node.warnings if w.code == "bars_do_not_fit"}

    assert not_fitting() == {"bottom", "top"}
    beam.set_longitudinal_rebar_bot(n1=1, d_b1=12 * mm)
    assert not_fitting() == {"top"}
    beam.set_longitudinal_rebar_top(n1=1, d_b1=12 * mm)
    assert not_fitting() == set()
    # A design run again reads the width again, and flags it again.
    node.design()
    assert not_fitting() == {"bottom", "top"}


def test_slab_bars_set_by_hand_clear_the_flag_of_their_face() -> None:
    """The slab setters detail a spacing, not a count, but they are the same
    hand: a face given bars by hand is no longer the face the search gave
    up on."""
    from mento import OneWaySlab
    from mento.units import m

    slab = OneWaySlab(
        label="L1",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=1 * m,
        height=20 * cm,
        c_c=25 * mm,
    )
    slab._infeasible_faces = {"bot", "top"}
    slab.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=15 * cm)
    assert slab._infeasible_faces == {"top"}
    slab.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=15 * cm)
    assert slab._infeasible_faces == set()


def test_a_face_with_one_bar_has_no_clear_spacing_to_miss() -> None:
    """10x30, 1Ø12 top and bottom, 1eØ10/10. A layer with one bar leaves
    100 - 2*(25 + 10) - 12 = 18 mm beside it, which is not a distance between
    bars: there is no pair to hold to the 25 mm of ACI 318-19 §25.2.1 (30 mm
    on top, the vibrator). bd94d2f reported both faces as "clear spacing
    between the bars below the minimum". Two Ø12 in the same web are 100 -
    70 - 24 = 6 mm apart and are warned, as before."""
    beam = _beam(width=10 * cm, height=30 * cm)
    beam.set_longitudinal_rebar_bot(n1=1, d_b1=12 * mm)
    beam.set_longitudinal_rebar_top(n1=1, d_b1=12 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=10 * mm, s_l=10 * cm)
    node = Node(section=beam, forces=[Forces(label="U", V_z=10 * kN, M_y=5 * kNm)])
    node.check()
    assert "clear_spacing_below_min" not in _by_code(node.warnings)

    beam.set_longitudinal_rebar_bot(n1=2, d_b1=12 * mm)
    spacing = _by_code(beam.warnings)["clear_spacing_below_min"]
    assert spacing.face == "bottom"
    assert spacing.values["s"].to("mm").magnitude == pytest.approx(6.0, abs=0.01)


def test_one_bar_per_layer_is_still_no_pair() -> None:
    """1Ø12 in the first layer and 1Ø12 in the second sit one above the other:
    neither layer has two bars side by side, so there is no clear spacing to
    report on the width."""
    beam = _beam(width=10 * cm, height=30 * cm)
    beam.set_longitudinal_rebar_bot(1, 12 * mm, 0, None, 1, 12 * mm)
    beam.set_longitudinal_rebar_top(n1=1, d_b1=12 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=10 * mm, s_l=10 * cm)
    assert "clear_spacing_below_min" not in _by_code(beam.warnings)


def test_clearing_the_stirrups_widens_the_space_for_the_bars() -> None:
    """Without stirrups the bars sit against the cover: (200 - 50 - 48)/3 = 34 mm."""
    beam = _beam(height=50 * cm)
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=12 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=16 * mm, s_l=20 * cm)
    beam.set_transverse_rebar(0, 0 * mm, 0 * cm)

    assert beam._available_s_bot.to("mm").magnitude == pytest.approx(34.0, abs=0.01)
    assert "clear_spacing_below_min" not in _by_code(beam.warnings)


# ---------------------------------------------------------------------------
# Stirrups that brace compression reinforcement, ACI 318-19 / CIRSOC 201-25 §9.7.6.4
# ---------------------------------------------------------------------------


def test_stirrups_of_a_doubly_reinforced_beam_are_held_to_its_compression_bars() -> None:
    """ACI 20x50 H25 ADN 420, Mu 260 kNm, Vu 60 kN: 2Ø20 + 1Ø20 in two layers below, 2Ø16 + 1Ø16 above.

    The top bars are compression steel. §9.7.6.4.3 caps the stirrups at
    min(16*16 = 256, 48*10 = 480, 200) = 200 mm, tighter than the d/2 =
    216 mm of Table 9.7.6.2.2, so the design details 1eØ10/20 where it used
    to detail 1eØ10/21. Set by hand at 21 cm the check says so, against the
    Ø16 bar; set at Ø8, below the No. 10 (9.5 mm) that §9.7.6.4.2(a) asks
    for a bar up to No. 32, it says that. Neither belongs to a combination:
    the bars are the section's.
    """
    beam = _beam(height=50 * cm)
    node = Node(section=beam, forces=[Forces(label="ELU", V_z=60 * kN, M_y=260 * kNm)])
    node.design()

    assert str(beam.reinforcement.top) == "2Ø16 mm + 1Ø16 mm"
    assert (beam.shear_design.d_b, beam.shear_design.s_l) == (10 * mm, 20 * cm)
    assert node.warnings == ()

    beam.set_transverse_rebar(n_stirrups=1, d_b=10 * mm, s_l=21 * cm)
    node.check()
    found = _by_code(node.warnings)
    spacing = found["stirrup_spacing_exceeds_compression_support"]
    assert (spacing.values["s"], spacing.values["s_max"], spacing.values["d_b_comp"]) == (21 * cm, 20 * cm, 16 * mm)
    assert spacing.combinations == ()
    assert "stirrup_spacing_exceeds_max" not in found  # Table 9.7.6.2.2 alone allows 21.6 cm
    assert "stirrup_diameter_below_compression_support" not in found

    beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=20 * cm)
    node.check()
    found = _by_code(node.warnings)
    diameter = found["stirrup_diameter_below_compression_support"]
    assert (diameter.values["d_b"], diameter.values["d_b_min"], diameter.values["d_b_comp"]) == (
        8 * mm,
        9.5 * mm,
        16 * mm,
    )
    assert "stirrup_spacing_exceeds_compression_support" not in found


def test_a_negative_moment_puts_the_braced_bars_at_the_bottom() -> None:
    """The same beam under Mu = -260 kNm: 2Ø16 + 1Ø16 at the bottom are the compression steel."""
    beam = _beam(height=50 * cm)
    node = Node(section=beam, forces=[Forces(label="ELU", V_z=60 * kN, M_y=-260 * kNm)])
    node.design()

    assert str(beam.reinforcement.bottom) == "2Ø16 mm + 1Ø16 mm"
    assert beam.shear_design.s_l == 20 * cm

    beam.set_transverse_rebar(n_stirrups=1, d_b=10 * mm, s_l=21 * cm)
    node.check()
    spacing = _by_code(node.warnings)["stirrup_spacing_exceeds_compression_support"]
    assert (spacing.values["s_max"], spacing.values["d_b_comp"]) == (20 * cm, 16 * mm)


def test_a_stirrup_that_makes_the_section_doubly_reinforced_is_spaced_for_it() -> None:
    """ACI 15x50 H40 ADN 420, c_c 25 mm, Mu = 211.5 kNm, Vu = 0: 2Ø25 + 2Ø20 in two layers, 2Ø10 on top.

    The flexure is designed at the depth of the Ø8 starter stirrup, where
    the section is singly reinforced: A_s,req = 14.94 cm² against the
    tension-controlled A_s,max = 14.98 cm² there. The shear design picks
    Ø10, which sinks the bars 2 mm: A_s,max falls to 14.92 cm², the 16.10
    cm² placed are past it, and the section now relies on the 2Ø10 on top
    as compression steel. §9.7.6.4.3 then caps the stirrups at min(16*10,
    48*10, 150) = 150 mm. Spaced with the compression faces of the Ø8 depth,
    the stirrups came out at 1eØ10/21: in silence before the cap existed,
    and with ``stirrup_spacing_exceeds_compression_support`` on the design
    itself once it did. Each diameter is now read with the compression steel
    its own depth relies on.
    """
    beam = RectangularBeam(
        label="V",
        concrete=Concrete_ACI_318_19(name="H40", f_c=40 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=15 * cm,
        height=50 * cm,
        c_c=25 * mm,
    )
    node = Node(section=beam, forces=[Forces(label="ELU", M_y=211.5 * kNm)])
    node.design()

    assert str(beam.reinforcement.bottom) == "2Ø25 mm + 2Ø20 mm"
    assert str(beam.reinforcement.top) == "2Ø10 mm"
    assert beam._compression_faces == {"top"}
    assert beam.shear_design.d_b == 10 * mm
    assert beam.shear_design.s_l <= 15 * cm
    assert node.warnings == ()


def test_cirsoc_grades_the_bracing_stirrup_with_the_compression_bar() -> None:
    """CIRSOC 20x40 H25 ADN 420, Mu 200 kNm, Vu 60 kN: 2Ø20 + 1Ø20 on top as compression steel.

    Tabla 9.7.6.4.2 asks 8 mm for a bar over 16 and up to 25 mm, so the
    design passes over the Ø6 its catalogue starts at -- it used to detail
    1eØ6/16 -- and details 1eØ8/16; §9.7.6.4.3 allows min(320, 384, 200) =
    200 mm, looser than the 16.7 cm of Tabla 9.7.6.2.2 at that depth. A Ø6
    set by hand is warned, against the Ø20 bar.
    """
    from mento import Concrete_CIRSOC_201_25

    beam = RectangularBeam(
        label="V",
        concrete=Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=20 * cm,
        height=40 * cm,
        c_c=25 * mm,
    )
    node = Node(section=beam, forces=[Forces(label="ELU", V_z=60 * kN, M_y=200 * kNm)])
    node.design()

    assert str(beam.reinforcement.top) == "2Ø20 mm + 1Ø20 mm"
    assert (beam.shear_design.d_b, beam.shear_design.s_l) == (8 * mm, 16 * cm)
    assert node.warnings == ()

    beam.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=16 * cm)
    node.check()
    found = _by_code(node.warnings)
    diameter = found["stirrup_diameter_below_compression_support"]
    assert (diameter.values["d_b"], diameter.values["d_b_min"], diameter.values["d_b_comp"]) == (
        6 * mm,
        8 * mm,
        20 * mm,
    )
    assert "stirrup_spacing_exceeds_compression_support" not in found
    mento.set_language("es")
    assert _by_code(node.warnings)["stirrup_diameter_below_compression_support"].message == (
        "Diámetro de estribo 6 mm menor que el mínimo 8 mm que exige el arriostramiento de barras comprimidas Ø20 mm."
    )


def test_a_singly_reinforced_beam_owes_its_stirrups_nothing_for_compression() -> None:
    """20x60 under +150 / -40 kNm needs no compression steel: the bracing limits do not apply, whatever the stirrups."""
    beam = _beam()
    node = Node(section=beam, forces=FORCES)
    node.design()
    beam.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=25 * cm)
    node.check()
    codes = set(_by_code(node.warnings))

    assert "stirrup_spacing_exceeds_max" in codes
    assert "stirrup_spacing_exceeds_compression_support" not in codes
    assert "stirrup_diameter_below_compression_support" not in codes


@pytest.mark.parametrize(
    "concrete, top, d_b_min",
    [
        (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), 16 * mm, 9.5 * mm),
        (mento.Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa), 20 * mm, 8 * mm),
    ],
    ids=["ACI", "CIRSOC"],
)
def test_a_doubly_reinforced_beam_with_no_stirrups_is_told_it_needs_them(
    concrete: Concrete_ACI_318_19, top: Quantity, d_b_min: Quantity
) -> None:
    """20x50 H25 ADN 420, no stirrups, Mu = 260 kNm, Vu = 0: 2Ø25 + 2Ø25 in two layers below, 3 bars above.

    With no stirrup the layers sit at 25 + 12.5 = 37.5 and 37.5 + 25 + 25 =
    87.5 mm, d = 500 - 62.5 = 437.5 mm, and the tension-controlled limit is
    A_s,max = 0.85*25*200*0.85*(0.003/0.0081)*437.5/420 = 13.94 cm²: the
    19.63 cm² below lean on the top bars as compression steel (flexure DCR
    0.92 under ACI, 0.89 under CIRSOC). §9.7.6.4.1 of both codes asks for
    transverse reinforcement wherever compression reinforcement is
    required, whatever the shear, and the check said nothing: with no
    stirrups it returned before reading the compression steel, and with
    Vu = 0 there was no ``stirrups_required`` either. Now it quotes the
    smallest stirrup §9.7.6.4.2 allows (No. 10 = 9.5 mm under ACI for a bar
    up to No. 32; 8 mm under CIRSOC Tabla 9.7.6.4.2 for a Ø20) and the
    spacing §9.7.6.4.3 gives it: min(16*16 = 256 or 16*20 = 320, 48*d_b,
    200) = 200 mm. Under 150 kNm the section is singly reinforced and owes
    nothing. A one-way slab is not held to it: §9.7.6.4 is the beams'.
    """
    beam = RectangularBeam(
        label="V",
        concrete=concrete,
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=20 * cm,
        height=50 * cm,
        c_c=25 * mm,
    )
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=25 * mm, n3=2, d_b3=25 * mm)
    beam.set_longitudinal_rebar_top(n1=3, d_b1=top)
    beam.set_transverse_rebar(0, 0 * mm, 0 * cm)
    node = Node(section=beam, forces=[Forces(label="ELU", M_y=260 * kNm)])
    node.check()

    assert beam.flexure_design.bottom.A_s_max.to("cm**2").magnitude == pytest.approx(13.94, abs=0.005)
    assert beam.flexure_design.DCR < 1
    found = _by_code(node.warnings)
    assert set(found) == {"stirrups_required_for_compression_support"}
    bracing = found["stirrups_required_for_compression_support"]
    assert (bracing.values["d_b_comp"], bracing.values["d_b_min"], bracing.values["s_max"]) == (top, d_b_min, 20 * cm)
    assert bracing.combinations == ()
    mento.set_language("es")
    assert _by_code(node.warnings)["stirrups_required_for_compression_support"].message == (
        f"La sección depende de barras comprimidas Ø{top.magnitude:g} mm y no tiene estribos que las arriostren: "
        f"hacen falta estribos cerrados de al menos {d_b_min.magnitude:g} mm separados a lo sumo 20 cm."
    )

    node = Node(section=beam, forces=[Forces(label="ELU", M_y=150 * kNm)])
    node.check()
    assert node.warnings == ()


def test_a_doubly_reinforced_slab_strip_owes_no_stirrups_for_its_compression_bars() -> None:
    """ACI 100x15 slab, Ø20/8 below and Ø12/15 above, Mu = 80 kNm: the top bars act in compression.

    ACI 318-19 §9.7.6.4 is a beam provision; a one-way slab sends its shear
    reinforcement to §9.7.6.2 alone (§7.7.5.1) and may be built with none.
    The strip carries Ø20/8 = 39.27 cm²/m against the tension-controlled
    A_s,max = 0.85*25*1000*0.85*(0.003/0.0081)*115/420 = 18.32 cm² at d =
    150 - 25 - 10 = 115 mm, so it relies on the top bars; no bracing
    warning is raised for it.
    """
    from mento import OneWaySlab
    from mento.units import m

    slab = OneWaySlab(
        label="L1",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=1 * m,
        height=15 * cm,
        c_c=25 * mm,
    )
    slab.set_slab_longitudinal_rebar_bot(d_b1=20 * mm, s_b1=8 * cm)
    slab.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=15 * cm)
    node = Node(section=slab, forces=[Forces(label="ELU", M_y=80 * kNm)])
    node.check()

    assert slab._compression_faces == {"top"}
    assert "stirrups_required_for_compression_support" not in _by_code(node.warnings)
