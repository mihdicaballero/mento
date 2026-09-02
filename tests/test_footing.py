"""Footing: the minimum longitudinal reinforcement of a member on the ground.

A :class:`~mento.slab.Footing` is a :class:`~mento.slab.OneWaySlab` in every
respect but one, so these tests are about that one respect. They pin the two
codes' minimum-reinforcement clauses to numbers, and they pin the property that
motivated the class: the design returns the largest applicable minimum already
applied, so a consumer never has to correct the engine's answer.
"""

import warnings

import pytest

from mento.codes.aci_318_19.equations import flexure as aci_flexure
from mento.codes.en_1992_2004.equations import flexure as en_flexure
from mento.forces import Forces
from mento.material import (
    Concrete_ACI_318_19,
    Concrete_CIRSOC_201_25,
    Concrete_EN_1992_2004,
    SteelBar,
)
from mento.node import Node
from mento.codes.EN_1992_2004_beam import _minimum_flexural_reinforcement_area_EN_1992_2004
from mento.slab import Footing, OneWaySlab
from mento.units import cm, inch, kNm, kip, m, mm, MPa, psi

# A metre-wide strip 600 mm deep: the geometry of the reference case the
# clauses below were verified against.
_WIDTH = 1 * m
_HEIGHT = 60 * cm
_COVER = 50 * mm


@pytest.fixture
def steel_b500s() -> SteelBar:
    return SteelBar(name="B500S", f_y=500 * MPa)


def _strip(cls, concrete, steel, label, height=_HEIGHT):  # type: ignore[no-untyped-def]
    return cls(
        label=label,
        concrete=concrete,
        steel_bar=steel,
        width=_WIDTH,
        height=height,
        c_c=_COVER,
    )


def _design(section, M_y=50 * kNm):  # type: ignore[no-untyped-def]
    Node(section=section, forces=[Forces(label="C1", M_y=M_y)]).design()
    return section


# ---------------------------------------------------------------------------
# The equations themselves
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "f_y, expected",
    [
        (420.0, 0.0018),  # the reference grade the table is written at
        (500.0, 0.0018 * 420 / 500),  # scales inversely above it
        (420.0 / 0.7, 0.0014),  # would scale below the floor, so the floor holds
        (700.0, 0.0014),
    ],
)
def test_aci_shrinkage_and_temperature_ratio(f_y: float, expected: float) -> None:
    """ACI 318-19 Table 24.4.3.2, metric."""
    assert aci_flexure.shrinkage_and_temperature_ratio(f_y) == pytest.approx(expected)


def test_aci_shrinkage_and_temperature_ratio_imperial() -> None:
    """The same table, anchored at 60 ksi instead of 420 MPa."""
    assert aci_flexure.shrinkage_and_temperature_ratio(60000.0, is_imperial=True) == pytest.approx(0.0018)
    assert aci_flexure.shrinkage_and_temperature_ratio(75000.0, is_imperial=True) == pytest.approx(0.00144)


@pytest.mark.parametrize(
    "h, expected",
    [
        (200.0, 1.00),  # saturated below 300 mm
        (300.0, 1.00),
        (600.0, 0.79),  # the reference case of the office
        (800.0, 0.65),
        (1200.0, 0.65),  # saturated above 800 mm
    ],
)
def test_en_crack_control_coefficient_k(h: float, expected: float) -> None:
    """EN 1992-1-1 §7.3.2(2): k interpolated with the depth."""
    assert en_flexure.crack_control_coefficient_k(h) == pytest.approx(expected)


@pytest.mark.parametrize(
    "f_yk, expected",
    [(300.0, 0.0010), (400.0, 0.0010), (450.0, 0.00095), (500.0, 0.0009), (600.0, 0.0009)],
)
def test_en_foundation_min_reinforcement_ratio(f_yk: float, expected: float) -> None:
    assert en_flexure.foundation_min_reinforcement_ratio(f_yk) == pytest.approx(expected)


def test_en_crack_control_min_reinforcement_is_equation_7_1() -> None:
    """A_s,min * sigma_s = k_c * k * f_ct,eff * A_ct, rearranged."""
    A_s_min = en_flexure.crack_control_min_reinforcement(0.4, 0.79, 2.56, 300_000.0, 500 / 1.15)
    assert A_s_min == pytest.approx(0.4 * 0.79 * 2.56 * 300_000.0 / (500 / 1.15))


# ---------------------------------------------------------------------------
# ACI 318-19 and CIRSOC 201-25
# ---------------------------------------------------------------------------


def test_aci_footing_minimum_is_the_gross_section_rule() -> None:
    """0.0018 * b * h at f_y = 420 MPa, and nothing about the effective depth."""
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    footing = _design(_strip(Footing, concrete, steel, "Z1"))

    expected = 0.0018 * (1000 * mm) * (600 * mm)
    assert footing._A_s_min_bot.to("cm**2").magnitude == pytest.approx(expected.to("cm**2").magnitude)


def test_aci_footing_minimum_is_lower_than_the_slab_minimum() -> None:
    """The whole point of the exemption: §9.6.1.1(b) relieves a member on the
    ground of the flexural minimum a slab spanning between supports carries."""
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    footing = _design(_strip(Footing, concrete, steel, "Z1"))
    slab = _design(_strip(OneWaySlab, concrete, steel, "L1"))

    assert footing._A_s_min_bot < slab._A_s_min_bot


def test_aci_footing_minimum_scales_with_the_steel_grade() -> None:
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    footing = _design(_strip(Footing, concrete, SteelBar(name="B500S", f_y=500 * MPa), "Z1"))

    expected = 0.0018 * 420 / 500 * (1000 * mm) * (600 * mm)
    assert footing._A_s_min_bot.to("cm**2").magnitude == pytest.approx(expected.to("cm**2").magnitude)


def test_cirsoc_shares_the_aci_clause() -> None:
    """CIRSOC 201-25 publishes ACI's formulas, this one included."""
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    aci = _design(_strip(Footing, Concrete_ACI_318_19(name="H25", f_c=25 * MPa), steel, "Z1"))
    cirsoc = _design(_strip(Footing, Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa), steel, "Z2"))

    assert cirsoc._A_s_min_bot.to("cm**2").magnitude == pytest.approx(aci._A_s_min_bot.to("cm**2").magnitude)


def test_aci_footing_minimum_in_imperial_units() -> None:
    """Same clause, anchored at 60 ksi, and still written on b*h."""
    concrete = Concrete_ACI_318_19(name="C4000", f_c=4000 * psi)
    steel = SteelBar(name="G60", f_y=60000 * psi)
    footing = Footing(
        label="Z1",
        concrete=concrete,
        steel_bar=steel,
        width=36 * inch,
        height=24 * inch,
        c_c=2 * inch,
    )
    _design(footing, M_y=40 * kip * inch)

    expected = 0.0018 * 36 * 24  # in²
    assert footing._A_s_min_bot.to("in**2").magnitude == pytest.approx(expected)


def test_aci_footing_with_no_moment_still_gets_the_minimum() -> None:
    """Shrinkage and temperature steel does not depend on a moment being there."""
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    footing = _design(_strip(Footing, concrete, steel, "Z1"), M_y=0 * kNm)

    expected = 0.0018 * (1000 * mm) * (600 * mm)
    assert footing._A_s_bot.to("cm**2").magnitude >= expected.to("cm**2").magnitude


# ---------------------------------------------------------------------------
# EN 1992-2004
# ---------------------------------------------------------------------------


def test_en_footing_minimum_is_the_crack_control_rule(steel_b500s: SteelBar) -> None:
    """The reference case: h = 600 mm, C25 (f_ctm = 2.56 MPa), B500S.

    k = 0.79 at that depth, k_c = 0.40 in bending, A_ct is the half of the
    section in tension and sigma_s is f_yd, so Eq. (7.1) gives the minimum --
    and it is above the geometric one, so it governs.
    """
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _design(_strip(Footing, concrete, steel_b500s, "Z1"))

    f_ctm = concrete.f_ctm.to(MPa).magnitude
    expected_mm2 = 0.4 * 0.79 * f_ctm * (1000 * 600 / 2) / (500 / 1.15)
    geometric_mm2 = 0.0009 * 1000 * 600
    assert expected_mm2 > geometric_mm2  # the crack-control rule is the governing one here

    assert footing._A_s_min_bot.to("mm**2").magnitude == pytest.approx(expected_mm2)


def test_en_footing_minimum_is_the_larger_of_the_two_rules(steel_b500s: SteelBar) -> None:
    """Which of the two rules governs depends on the depth.

    Both are proportional to h, so only k separates them: it decays past
    300 mm, and by 900 mm it has saturated at 0.65 and the crack-control
    minimum has fallen below the geometric one, which then governs.
    """
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _design(_strip(Footing, concrete, steel_b500s, "Z1", height=90 * cm), M_y=100 * kNm)

    f_ctm = concrete.f_ctm.to(MPa).magnitude
    crack_mm2 = 0.4 * 0.65 * f_ctm * (1000 * 900 / 2) / (500 / 1.15)
    geometric_mm2 = 0.0009 * 1000 * 900
    assert geometric_mm2 > crack_mm2  # the geometric rule is the governing one here

    assert footing._A_s_min_bot.to("mm**2").magnitude == pytest.approx(max(crack_mm2, geometric_mm2))


def test_en_footing_minimum_is_lower_than_the_slab_minimum(steel_b500s: SteelBar) -> None:
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _design(_strip(Footing, concrete, steel_b500s, "Z1"))
    slab = _design(_strip(OneWaySlab, concrete, steel_b500s, "L1"))

    assert footing._A_s_min_bot < slab._A_s_min_bot


# ---------------------------------------------------------------------------
# What the class is, and what a design returns
# ---------------------------------------------------------------------------


def test_footing_is_a_one_way_slab(steel_b500s: SteelBar) -> None:
    """Same reinforcement parameterisation, same checks -- only the minimum differs."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1")

    assert isinstance(footing, OneWaySlab)
    assert footing.mode == "slab"
    assert hasattr(footing, "set_slab_longitudinal_rebar_bot")


def test_support_says_which_clause_applies(steel_b500s: SteelBar) -> None:
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    assert _strip(Footing, concrete, steel_b500s, "Z1").support == "soil"
    assert _strip(OneWaySlab, concrete, steel_b500s, "L1").support == "free"


@pytest.mark.parametrize(
    "concrete, steel",
    [
        (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), SteelBar(name="ADN 420", f_y=420 * MPa)),
        (Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), SteelBar(name="B500S", f_y=500 * MPa)),
    ],
    ids=["aci", "en"],
)
def test_design_already_covers_the_minimum(concrete, steel) -> None:  # type: ignore[no-untyped-def]
    """The reason this lives in mento: what the design returns needs no
    correction by the consumer. A small moment is used so the minimum, and not
    the demand, is what sizes the steel."""
    footing = _design(_strip(Footing, concrete, steel, "Z1"), M_y=10 * kNm)

    assert footing._A_s_bot >= footing._A_s_min_bot
    assert footing._A_s_bot >= footing._A_s_req_bot


def test_a_footing_still_carries_its_moment(steel_b500s: SteelBar) -> None:
    """The relaxed minimum only ever relaxes the minimum: a demand above it is
    still what governs, and the designed section resists it."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1")
    combo = [Forces(label="C1", M_y=300 * kNm)]
    Node(section=footing, forces=combo).design()

    assert footing._A_s_bot > footing._A_s_min_bot
    assert footing.flexure_check_results(combo)[0].bottom.DCR <= 1.0


# ---------------------------------------------------------------------------
# Detailing: how far apart the bars sit
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "concrete, steel",
    [
        (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), SteelBar(name="ADN 420", f_y=420 * MPa)),
        (Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), SteelBar(name="B500S", f_y=500 * MPa)),
    ],
    ids=["aci", "en"],
)
def test_footing_bars_are_capped_at_300_mm(concrete, steel) -> None:  # type: ignore[no-untyped-def]
    """3h stops binding on a section this thick, so the footing cap is what holds."""
    footing = _strip(Footing, concrete, steel, "Z1")

    assert footing._max_bar_spacing().to("mm").magnitude == pytest.approx(300.0)


@pytest.mark.parametrize(
    "concrete, steel",
    [
        (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), SteelBar(name="ADN 420", f_y=420 * MPa)),
        (Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), SteelBar(name="B500S", f_y=500 * MPa)),
    ],
    ids=["aci", "en"],
)
def test_footing_bars_are_floored_at_100_mm(concrete, steel) -> None:  # type: ignore[no-untyped-def]
    assert _strip(Footing, concrete, steel, "Z1")._min_bar_spacing().to("mm").magnitude == pytest.approx(100.0)


def test_a_slab_keeps_the_wider_slab_limits(steel_b500s: SteelBar) -> None:
    """The tighter pair is a footing's, not every slab's."""
    slab = _strip(OneWaySlab, Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), steel_b500s, "L1")

    assert slab._max_bar_spacing().to("mm").magnitude == pytest.approx(400.0)
    assert slab._min_bar_spacing() is None


@pytest.mark.parametrize("M_y", [50, 150, 300, 450, 600, 750, 900, 1100])
def test_a_designed_footing_stays_inside_the_spacing_range(steel_b500s: SteelBar, M_y: int) -> None:
    """Across the range a metre strip can carry, the design answers with a
    bigger bar rather than with bars closer than 100 mm."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _design(_strip(Footing, concrete, steel_b500s, "Z1"), M_y=M_y * kNm)

    spacing = footing._s_b1_b.to("mm").magnitude
    assert 100.0 <= spacing <= 300.0
    assert footing.flexure_check_results([Forces(label="C1", M_y=M_y * kNm)])[0].bottom.DCR <= 1.0


def test_the_floor_is_what_holds_the_footing_apart(steel_b500s: SteelBar) -> None:
    """The same demand on a slab, which has no floor, is detailed tighter."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _design(_strip(Footing, concrete, steel_b500s, "Z1"), M_y=900 * kNm)
    slab = _design(_strip(OneWaySlab, concrete, steel_b500s, "L1"), M_y=900 * kNm)

    assert slab._s_b1_b.to("mm").magnitude < 100.0
    assert footing._s_b1_b.to("mm").magnitude >= 100.0


def test_the_spacing_row_of_a_footing_reports_both_bounds(steel_b500s: SteelBar) -> None:
    """The check table carries the range, so a footing detailed outside it fails."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _design(_strip(Footing, concrete, steel_b500s, "Z1"))

    row = footing._data_min_max_flexure
    index = row["Check"].index("Bar spacing bottom")
    assert row["Min."][index] == pytest.approx(100.0)
    assert row["Max."][index] == pytest.approx(300.0)


# ---------------------------------------------------------------------------
# Detailing: how thin the section may be
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "concrete, steel, height, limit",
    [
        (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), SteelBar(name="ADN 420", f_y=420 * MPa), 18 * cm, "200"),
        (Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), SteelBar(name="B500S", f_y=500 * MPa), 22 * cm, "250"),
    ],
    ids=["aci", "en"],
)
def test_a_thin_footing_warns(concrete, steel, height, limit) -> None:  # type: ignore[no-untyped-def]
    with pytest.warns(UserWarning, match=f"below the {limit} mm"):
        _strip(Footing, concrete, steel, "Z1", height=height)


@pytest.mark.parametrize(
    "concrete, steel",
    [
        (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), SteelBar(name="ADN 420", f_y=420 * MPa)),
        (Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), SteelBar(name="B500S", f_y=500 * MPa)),
    ],
    ids=["aci", "en"],
)
def test_a_footing_of_the_usual_depth_does_not_warn(concrete, steel) -> None:  # type: ignore[no-untyped-def]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _strip(Footing, concrete, steel, "Z1")


def test_a_thin_slab_does_not_warn(steel_b500s: SteelBar) -> None:
    """The thickness rule is a footing's; a suspended slab is as thin as it is."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _strip(OneWaySlab, Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), steel_b500s, "L1", height=15 * cm)


def test_a_thin_footing_is_still_designed(steel_b500s: SteelBar) -> None:
    """The warning is advice. The thickness is the engineer's to choose, and the
    design that follows is still the right design for the section given."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    with pytest.warns(UserWarning):
        footing = _strip(Footing, concrete, steel_b500s, "Z1", height=20 * cm)
    _design(footing, M_y=20 * kNm)

    assert footing._A_s_bot >= footing._A_s_min_bot
    assert footing._A_s_bot.magnitude > 0


def test_the_warning_names_the_footing_and_the_code() -> None:
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    with pytest.warns(UserWarning) as caught:
        _strip(Footing, concrete, steel, "Z7", height=15 * cm)

    message = str(caught[0].message)
    assert "Z7" in message
    assert "ACI 318-19" in message


# ---------------------------------------------------------------------------
# A section with nothing on the tension face
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("face", ["bot", "top"], ids=["bottom", "top"])
@pytest.mark.parametrize(
    "concrete, steel",
    [
        (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), SteelBar(name="ADN 420", f_y=420 * MPa)),
        (Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), SteelBar(name="B500S", f_y=500 * MPa)),
    ],
    ids=["aci", "en"],
)
def test_an_unreinforced_face_reports_a_failing_dcr(concrete, steel, face) -> None:  # type: ignore[no-untyped-def]
    """Both codes answer a moment on a bare face with a DCR far above 1.

    EN used to divide by zero here, on either face. It is what a footing too
    thin to be reinforced within its spacing range now reports, so it has to be
    a number -- and a slab reinforced on one face only is checked for the
    combination that puts the other one in tension, so both branches are real.
    """
    section = _strip(OneWaySlab, concrete, steel, "L1")
    # A positive moment puts the bottom in tension, a negative one the top.
    combo = [Forces(label="C1", M_y=(50 if face == "bot" else -50) * kNm)]

    result = section.flexure_check_results(combo)[0]
    assert getattr(result, "bottom" if face == "bot" else "top").DCR > 1


# ---------------------------------------------------------------------------
# A design code that states none of these rules
# ---------------------------------------------------------------------------


def test_a_code_without_the_footing_rules_imposes_none() -> None:
    """A rule a code does not have is not a rule an element can fail.

    ``min_bar_spacing_slab`` and ``min_thickness_on_soil`` are optional hooks:
    both registered codes fill them in, so this registers one that does not and
    drives a footing through it. Nothing is capped and nothing is warned about
    -- silence, rather than an error about a limit that was never stated.
    """
    import dataclasses

    from mento.codes.registry import _REGISTRY, DesignCode, design_code, register

    reference = design_code(Concrete_ACI_318_19(name="H25", f_c=25 * MPa))
    invented = DesignCode(
        **{
            **{f.name: getattr(reference, f.name) for f in dataclasses.fields(reference)},
            "title": "NBR 6118-2023",
            "year": 2023,
            "max_bar_spacing_slab": None,
            "min_bar_spacing_slab": None,
            "min_thickness_on_soil": None,
        }
    )
    register(invented)
    try:
        concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
        concrete.design_code = invented.title
        steel = SteelBar(name="ADN 420", f_y=420 * MPa)

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            # Thin enough that either registered code would have warned.
            footing = _strip(Footing, concrete, steel, "Z1", height=15 * cm)

        assert footing._min_bar_spacing() is None
        assert footing._max_bar_spacing() is None
        # And the bar count is left at the slab default, since nothing bounds it.
        assert footing.settings.max_bars_per_layer >= 200
        # With no bounds there is no space to search a mat in, so none is offered
        # and the design keeps whatever the per-face search arrived at.
        assert list(footing._mat_candidates(10 * cm**2, 5 * cm**2)) == []
    finally:
        _REGISTRY.pop(invented.title, None)


# ---------------------------------------------------------------------------
# One mat: both faces at one spacing
# ---------------------------------------------------------------------------


def _design_envelope(section, M_bot, M_top):  # type: ignore[no-untyped-def]
    """Design for an envelope with a moment on each face, and return it."""
    combo = [Forces(label="C1", M_y=M_bot), Forces(label="C2", M_y=-M_top)]
    Node(section=section, forces=combo).design()
    return section, combo


@pytest.mark.parametrize(
    "concrete, steel",
    [
        (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), SteelBar(name="ADN 420", f_y=420 * MPa)),
        (Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), SteelBar(name="B500S", f_y=500 * MPa)),
    ],
    ids=["aci", "en"],
)
def test_a_designed_footing_is_one_mat(concrete, steel) -> None:  # type: ignore[no-untyped-def]
    """Both faces come out on one module: the top at the bottom's spacing, or at
    exactly twice it so one top bar lands on every second bottom bar.

    The bottom carries five times the top's moment here, so the two faces would
    otherwise be detailed at spacings nobody would draw on one mat.
    """
    footing, _ = _design_envelope(_strip(Footing, concrete, steel, "Z1"), 600 * kNm, 120 * kNm)

    s_bot = footing._s_b1_b.to("mm").magnitude
    s_top = footing._s_b1_t.to("mm").magnitude
    assert s_top / s_bot in (1.0, 2.0)


def test_the_fallback_takes_the_smaller_of_the_two_spacings(steel_b500s: SteelBar) -> None:
    """What the design falls back to when no searched mat verifies.

    Closing the wider face up only ever adds bars, which is why it needs no
    verification of its own. Driven from a layout set by hand, so that the two
    spacings going in are known and the rule can be read off the result.
    """
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1")
    footing.set_slab_longitudinal_rebar_bot(d_b1=16 * mm, s_b1=120 * mm)
    footing.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=250 * mm)
    A_s_top_before = footing._A_s_top

    footing._match_spacings_downward()

    assert footing._s_b1_b.to("mm").magnitude == pytest.approx(120.0)
    assert footing._s_b1_t.to("mm").magnitude == pytest.approx(120.0)
    # Closing the top up added bars, and left its bar alone.
    assert footing._A_s_top > A_s_top_before
    assert footing._d_b1_t.to("mm").magnitude == pytest.approx(12.0)


def test_the_mat_is_lighter_than_reconciling_the_two_spacings(steel_b500s: SteelBar) -> None:
    """Why the search exists: taking the smaller of the two spacings buys steel
    the section never asked for, and choosing the module instead does not."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    searched, _ = _design_envelope(_strip(Footing, concrete, steel_b500s, "Z1"), 600 * kNm, 120 * kNm)

    reconciled, _ = _design_envelope(_strip(Footing, concrete, steel_b500s, "Z2"), 600 * kNm, 120 * kNm)
    reconciled._match_spacings_downward()

    assert searched._A_s_bot + searched._A_s_top <= reconciled._A_s_bot + reconciled._A_s_top


def test_matching_the_mat_never_undoes_the_design(steel_b500s: SteelBar) -> None:
    """Closing the wider face up only adds steel, so both faces still pass."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing, combo = _design_envelope(_strip(Footing, concrete, steel_b500s, "Z1"), 600 * kNm, 120 * kNm)

    results = footing.flexure_check_results(combo)
    assert max(r.bottom.DCR for r in results) <= 1.0
    assert max(r.top.DCR for r in results) <= 1.0
    assert footing._A_s_bot >= footing._A_s_min_bot
    assert footing._A_s_top >= footing._A_s_min_top


def test_the_mat_keeps_each_face_within_the_spacing_range(steel_b500s: SteelBar) -> None:
    """Matching down cannot push a face below the 100 mm floor: the module is
    one of the two spacings the design already chose."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing, _ = _design_envelope(_strip(Footing, concrete, steel_b500s, "Z1"), 900 * kNm, 100 * kNm)

    for spacing in (footing._s_b1_b, footing._s_b1_t):
        assert 100.0 <= spacing.to("mm").magnitude <= 300.0


def test_each_face_keeps_its_own_diameter(steel_b500s: SteelBar) -> None:
    """Only the module is shared. A face carrying a fifth of the moment is
    still allowed its own, thinner bar."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing, _ = _design_envelope(_strip(Footing, concrete, steel_b500s, "Z1"), 600 * kNm, 120 * kNm)

    assert footing._d_b1_t <= footing._d_b1_b


def test_no_mat_is_detailed_below_the_practical_bar(steel_b500s: SteelBar) -> None:
    """A footing mesh is not drawn with the thinnest bar in the catalogue, even
    where the arithmetic would allow it on the lightly loaded face."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing, _ = _design_envelope(_strip(Footing, concrete, steel_b500s, "Z1"), 600 * kNm, 20 * kNm)

    assert footing._d_b1_t.to("mm").magnitude >= 10.0
    assert footing._d_b1_b.to("mm").magnitude >= 10.0


def test_a_footing_reinforced_on_one_face_gets_no_second_grid(steel_b500s: SteelBar) -> None:
    """Nothing to reconcile, and a grid is not invented to reconcile with."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _design(_strip(Footing, concrete, steel_b500s, "Z1"), M_y=300 * kNm)

    assert footing._s_b1_b.magnitude > 0
    assert footing._s_b1_t.magnitude == 0
    assert footing._d_b1_t.magnitude == 0


def test_a_slab_still_details_its_faces_independently(steel_b500s: SteelBar) -> None:
    """The mat is a footing's; a suspended slab keeps the spacing each face
    actually needed."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    slab, _ = _design_envelope(_strip(OneWaySlab, concrete, steel_b500s, "L1"), 600 * kNm, 120 * kNm)

    assert slab._s_b1_b != slab._s_b1_t


def test_reconciling_an_existing_mat_is_a_no_op(steel_b500s: SteelBar) -> None:
    """The fallback run on a mat that is already one leaves it alone."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1")
    footing.set_slab_longitudinal_rebar_bot(d_b1=16 * mm, s_b1=150 * mm)
    footing.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=150 * mm)
    before = (footing._s_b1_b, footing._s_b1_t, footing._A_s_bot, footing._A_s_top)

    footing._match_spacings_downward()

    assert (footing._s_b1_b, footing._s_b1_t, footing._A_s_bot, footing._A_s_top) == before


def test_the_top_may_be_set_out_at_twice_the_bottom(steel_b500s: SteelBar) -> None:
    """The modular option is real and gets used.

    A top face carrying a tenth of the bottom's moment is drawn at double the
    bottom's module -- one top bar on every second bottom bar -- rather than
    matched to it and given steel it has no use for.
    """
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    doubled = 0
    for h, M_bot, M_top in [(50, 300, 30), (60, 400, 40), (40, 200, 20), (70, 600, 60)]:
        footing, _ = _design_envelope(
            _strip(Footing, concrete, steel_b500s, "Z1", height=h * cm), M_bot * kNm, M_top * kNm
        )
        ratio = footing._s_b1_t.to("mm").magnitude / footing._s_b1_b.to("mm").magnitude
        assert ratio in (1.0, 2.0)
        doubled += ratio == 2.0

    assert doubled, "the double module was never taken, so the option is not doing anything"


# ---------------------------------------------------------------------------
# The edges of the mat search
# ---------------------------------------------------------------------------


def test_a_section_needing_nothing_is_offered_no_mat(steel_b500s: SteelBar) -> None:
    """No steel required on either face is not a mat problem."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1")

    assert list(footing._mat_candidates(0 * cm**2, 0 * cm**2)) == []


def test_a_face_no_bar_can_cover_drops_that_module(steel_b500s: SteelBar) -> None:
    """A top face wanting more steel than the widest module can carry does not
    veto the search -- that module is simply not one of the candidates."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1")

    # 60 cm² over the four bars a 300 mm module fits in the strip is beyond any
    # bar in the catalogue, so the wide modules drop out and the close ones stay.
    candidates = list(footing._mat_candidates(20 * cm**2, 60 * cm**2))
    assert candidates
    assert all(s_top <= 150 for (_, s_top), _ in candidates)


def test_the_fallback_leaves_a_single_grid_alone(steel_b500s: SteelBar) -> None:
    """Reconciling needs two grids; with one there is nothing to reconcile."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1")
    footing.set_slab_longitudinal_rebar_bot(d_b1=16 * mm, s_b1=150 * mm)
    before = (footing._s_b1_b, footing._d_b1_b, footing._A_s_bot)

    footing._match_spacings_downward()

    assert (footing._s_b1_b, footing._d_b1_b, footing._A_s_bot) == before


def test_a_mat_that_does_not_verify_is_passed_over(steel_b500s: SteelBar) -> None:
    """The search proposes, the capacity check disposes.

    Driven with a verifier that refuses everything, so the design falls back to
    reconciling the two spacings it had already proved -- which is the point of
    keeping that rule around.
    """
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1")
    footing.set_slab_longitudinal_rebar_bot(d_b1=16 * mm, s_b1=120 * mm)
    footing.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=250 * mm)

    footing._finalize_longitudinal_design(20 * cm**2, 5 * cm**2, lambda: False)

    # Not a searched mat: the layout it started from, with the spacings matched.
    assert footing._d_b1_b.to("mm").magnitude == pytest.approx(16.0)
    assert footing._d_b1_t.to("mm").magnitude == pytest.approx(12.0)
    assert footing._s_b1_b == footing._s_b1_t == 120 * mm


@pytest.mark.parametrize(
    "height, M_bot, M_top",
    [(25 * cm, 200 * kNm, -20 * kNm), (30 * cm, 100 * kNm, -600 * kNm)],
    ids=["bottom-short", "top-short"],
)
def test_a_section_that_cannot_carry_the_moment_falls_back(steel_b500s, height, M_bot, M_top) -> None:  # type: ignore[no-untyped-def]
    """No mat verifies when the section is too small for the demand on a face.

    Each candidate is applied and rejected in turn -- the first case runs the
    bottom out of capacity, the second the top -- and the design still returns
    a layout, leaving the check to report DCR > 1 rather than raising.
    """
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1", height=height)
    combo = [Forces(label="C1", M_y=M_bot), Forces(label="C2", M_y=M_top)]
    Node(section=footing, forces=combo).design()

    results = footing.flexure_check_results(combo)
    assert max(max(r.bottom.DCR, r.top.DCR) for r in results) > 1


# ---------------------------------------------------------------------------
# Which of the two EN rules governs, and at what depth
# ---------------------------------------------------------------------------


def _en_minimum_per_mille(footing):  # type: ignore[no-untyped-def]
    """A_s,min on the tension face of a metre strip, per mille of the gross section."""
    d = (footing.height - footing.c_c - 10 * mm).to("mm").magnitude
    area = _minimum_flexural_reinforcement_area_EN_1992_2004(footing, d)
    gross = (footing.width * footing.height).to("mm**2").magnitude
    return area / gross * 1000


@pytest.mark.parametrize(
    "height, expected",
    [(40 * cm, 1.097), (60 * cm, 0.932), (90 * cm, 0.900)],
    ids=["h040", "h060", "h090"],
)
def test_the_crack_rule_governs_the_thin_footings_not_the_thick(steel_b500s, height, expected) -> None:  # type: ignore[no-untyped-def]
    """Crack control is written on A_ct = b*h/2, so its ratio goes with k/2, and
    k decays from 1.00 to 0.65 between 300 and 800 mm.

    Both rules are proportional to h, so only k separates them: crack control
    governs while the section is thin, and the 0.9 per mille geometric minimum
    takes over once k has decayed — by 0.90 m it already has.
    """
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1", height=height)

    assert _en_minimum_per_mille(footing) == pytest.approx(expected, abs=0.001)


def test_the_excess_over_the_geometric_minimum_shrinks_with_depth(steel_b500s: SteelBar) -> None:
    """The same thing read as a trend rather than three numbers."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    excess = [
        _en_minimum_per_mille(_strip(Footing, concrete, steel_b500s, "Z1", height=h * cm)) - 0.900 for h in (40, 60, 90)
    ]

    assert excess[0] > excess[1] > 0
    assert excess[2] == pytest.approx(0.0, abs=1e-9)


def test_a_footing_with_no_top_moment_is_left_without_top_steel(steel_b500s: SteelBar) -> None:
    """A face that is not in tension is given no minimum, not a smaller one.

    Whether a footing carries top steel is the consumer's call — it is the only
    one that sees both orthogonal sections of the element — so a design given no
    negative moment leaves the top empty rather than reasoning about what
    minimum an unloaded face would be owed.
    """
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1", height=40 * cm)
    combo = [Forces(label="C1", M_y=200 * kNm)]
    Node(section=footing, forces=combo).design()

    assert footing._A_s_top.magnitude == 0
    assert footing._d_b1_t.magnitude == 0
    assert footing._s_b1_t.magnitude == 0
    assert max(r.top.A_s_min for r in footing.flexure_check_results(combo)) == 0


def test_both_faces_in_tension_both_get_the_crack_minimum(steel_b500s: SteelBar) -> None:
    """A footing bending both ways has a tension zone on each face, so each is
    owed the crack-control minimum."""
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    footing = _strip(Footing, concrete, steel_b500s, "Z1", height=40 * cm)
    combo = [Forces(label="M+", M_y=150 * kNm), Forces(label="M-", M_y=-60 * kNm)]
    Node(section=footing, forces=combo).design()

    results = footing.flexure_check_results(combo)
    gross = (footing.width * footing.height).to("mm**2").magnitude
    for face in ("bottom", "top"):
        governing = max(getattr(r, face).A_s_min for r in results)
        assert governing / gross * 1000 == pytest.approx(1.097, abs=0.001)
