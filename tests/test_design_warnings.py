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
        "stirrup_diameter_below_min",
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

    diameter = found["stirrup_diameter_below_min"].values
    assert (diameter["d_b"], diameter["d_b_min"]) == (6 * mm, 10 * mm)


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
    # EN 1992-1-1 states no minimum stirrup diameter, so none is reported.
    assert "stirrup_diameter_below_min" not in codes


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
    assert warning.values["A_s_min"].to("cm**2").magnitude == pytest.approx(top.A_s_min_eff.to("cm**2").magnitude)


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
