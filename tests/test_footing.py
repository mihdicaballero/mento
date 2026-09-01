"""Footing: the minimum longitudinal reinforcement of a member on the ground.

A :class:`~mento.slab.Footing` is a :class:`~mento.slab.OneWaySlab` in every
respect but one, so these tests are about that one respect. They pin the two
codes' minimum-reinforcement clauses to numbers, and they pin the property that
motivated the class: the design returns the largest applicable minimum already
applied, so a consumer never has to correct the engine's answer.
"""

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
