"""Structured warnings: the detailing limits a section misses, as data."""

from types import MappingProxyType
from typing import Iterator

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
    """A 12x30 web takes 4Ø12 at most, and 40 kNm asks for 5.18 cm² below."""
    beam = RectangularBeam(
        label="101",
        concrete=mento.Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=12 * cm,
        height=30 * cm,
        c_c=25 * mm,
    )
    node = Node(section=beam, forces=[Forces(label="1.4D", V_z=50 * kN, M_y=40 * kNm)])
    node.design()

    short = _by_code(node.warnings)["As_below_required"]
    assert short.face == "bottom"
    assert short.values["A_s"].to("cm**2").magnitude == pytest.approx(4.52, rel=1e-3)
    assert short.values["A_s_req"].to("cm**2").magnitude == pytest.approx(5.18, rel=1e-2)
    mento.set_language("es")
    assert "agrandar la sección" in _by_code(node.warnings)["As_below_required"].message

    # Bars set by hand that reach the area clear it.
    beam.set_longitudinal_rebar_bot(2, 16 * mm, 0, None, 2, 16 * mm)
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


def test_clearing_the_stirrups_widens_the_space_for_the_bars() -> None:
    """Without stirrups the bars sit against the cover: (200 - 50 - 48)/3 = 34 mm."""
    beam = _beam(height=50 * cm)
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=12 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=16 * mm, s_l=20 * cm)
    beam.set_transverse_rebar(0, 0 * mm, 0 * cm)

    assert beam._available_s_bot.to("mm").magnitude == pytest.approx(34.0, abs=0.01)
    assert "clear_spacing_below_min" not in _by_code(beam.warnings)
