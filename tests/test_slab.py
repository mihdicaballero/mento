import math

import pytest
from pint import Quantity

from mento.beam import RectangularBeam
from mento.node import Node
from mento.slab import OneWaySlab, _bars_at_spacing
from mento.material import Concrete_ACI_318_19, SteelBar, Concrete_EN_1992_2004
from mento.units import kip, inch, mm, cm, kN, MPa, kNm, ksi
from mento.forces import Forces


@pytest.fixture()
def slab_example_EN_1992_2004() -> OneWaySlab:
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    steelBar = SteelBar(name="B500S", f_y=420 * MPa)
    slab = OneWaySlab(
        label="Slab 01",
        concrete=concrete,
        steel_bar=steelBar,
        width=100 * cm,
        height=20 * cm,
        c_c=2.5 * cm,
    )
    return slab


@pytest.fixture()
def slab_example_ACI_318_19() -> OneWaySlab:
    concrete = Concrete_ACI_318_19(name="H25", f_c=4 * ksi)
    steelBar = SteelBar(name="ADN 420", f_y=60 * ksi)
    slab = OneWaySlab(
        label="Slab 02",
        concrete=concrete,
        steel_bar=steelBar,
        width=12 * inch,
        height=7 * inch,
        c_c=0.75 * inch,
    )
    return slab


@pytest.fixture()
def slab_example_ACI_318_19_metric() -> OneWaySlab:
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    slab = OneWaySlab(
        label="Slab 03",
        concrete=concrete,
        steel_bar=steelBar,
        width=100 * cm,
        height=20 * cm,
        c_c=2.5 * cm,
    )
    return slab


def build_slab() -> OneWaySlab:
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    return OneWaySlab(
        label="Slab extra",
        concrete=concrete,
        steel_bar=steel,
        width=120 * cm,
        height=20 * cm,
        c_c=20 * mm,
    )


def test_slab_initialization_sets_slab_mode_and_defaults() -> None:
    slab = build_slab()

    assert slab.mode == "slab"
    assert slab.reinforcement.transverse.d_b.magnitude == 0
    assert slab.settings.max_bars_per_layer >= 200


def test_set_slab_transverse_rebar_computes_shear_area() -> None:
    slab = build_slab()

    slab.set_slab_transverse_rebar(d_b=10 * mm, s_long=20 * cm, s_trans=25 * cm)

    transverse = slab.reinforcement.transverse
    assert transverse.A_v.magnitude > 0
    assert transverse.s_l == 20 * cm


def test_set_slab_transverse_rebar_zero_spacing_clears_shear_area() -> None:
    slab = build_slab()
    slab.set_slab_transverse_rebar(d_b=10 * mm, s_long=20 * cm, s_trans=25 * cm)

    # Clearing the stirrups must not divide by the zero spacing
    slab.set_slab_transverse_rebar()

    transverse = slab.reinforcement.transverse
    assert transverse.s_l.magnitude == 0
    assert transverse.A_v.to("cm**2/m").magnitude == 0


def test_longitudinal_rebar_spacing_updates_counts() -> None:
    slab = build_slab()

    slab.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=30 * cm, d_b3=10 * mm, s_b3=25 * cm)
    slab.set_slab_longitudinal_rebar_top(d_b1=8 * mm, s_b1=40 * cm)

    # A slab only ever fills positions 1 and 3, so the layers of a face are
    # position 1 first and position 3 second, with the empty ones omitted.
    bottom = slab.reinforcement.bottom
    top = slab.reinforcement.top
    assert bottom.layers[0].n == 4  # ceil(120/30)
    assert bottom.layers[1].n == 5  # ceil(120/25)
    assert top.layers[0].n == 3  # ceil(120/40)
    assert len(top.layers) == 1  # position 3 spacing left as default zero

    # ensure defaults are preserved when zero values are provided
    # The per-position spacing has no equivalent on the public reinforcement view.
    previous_spacing = slab._s_b1_t
    slab.set_slab_longitudinal_rebar_top(s_b3=0 * mm)
    assert slab._s_b1_t == previous_spacing


def test_shear_check_ACI_318_19_1(slab_example_ACI_318_19: OneWaySlab) -> None:
    # Example from Two-Way Flat Plate Concrete Floor System Analysis and Design (ACI 318-14) adjusted to ACI 318-19.
    # With guidance from CRSI Design Guide on ACI 318-19
    # See calcpad: ACI 318-19 Slab Shear 01 - Imperial.cpd
    f = Forces(V_z=1.52 * kip, N_x=0 * kip)
    node = Node(section=slab_example_ACI_318_19, forces=f)
    slab_example_ACI_318_19.set_slab_longitudinal_rebar_bot(d_b1=0.5 * inch, s_b1=10 * inch)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    # d_shear = 6.00 in = 152.4 mm, so the size effect factor is capped:
    # lambda_s = min(sqrt(2/(1 + d/10in)), 1.0) = min(1.118, 1.0) = 1.0
    # (ACI 318-19 Eq. 22.5.5.1.3). The factor only ever reduces V_c; the
    # reference sheet omitted the cap, which gave lambda_s = 1.118 and
    # phi*V_c = 23.92 kN instead of 21.39 kN.
    assert results.iloc[1]["Av,min"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["ØVc"] == pytest.approx(21.39, rel=1e-3)
    assert results.iloc[1]["ØVs"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["ØVn"] == pytest.approx(21.39, rel=1e-3)
    # phi*V_max = phi_v*(V_c + 0.66*lambda*sqrt(f_c)*A_cv) carries V_c, so it
    # drops by the same 2.53 kN: 145.45 -> 142.93.
    assert results.iloc[1]["ØVmax"] == pytest.approx(142.93, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.317, rel=1e-2)

    # Assert non-numeric values directly
    assert results.iloc[1]["Vu≤ØVmax"] is True
    assert results.iloc[1]["Vu≤ØVn"] is True


def test_check_flexure_ACI_318_19_1(slab_example_ACI_318_19_metric: OneWaySlab) -> None:
    # Testing the check of the reinforced slab with simple reinforcement
    # See calcpad: ACI 318-19 Slab Flexure 01 - Metric.cpd
    f = Forces(label="C1", M_y=20 * kNm)
    slab_example_ACI_318_19_metric.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=20 * cm)
    node = Node(section=slab_example_ACI_318_19_metric, forces=f)
    results = node.check_flexure()

    print(results)

    assert results.iloc[1]["Label"] == "Slab 03"
    assert results.iloc[1]["Comb."] == "C1"
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["As,min"] == pytest.approx(5.63, rel=1e-2)
    assert results.iloc[1]["As,req bot"] == pytest.approx(4.25, rel=1e-3)
    assert results.iloc[1]["As,req top"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["As"] == pytest.approx(5.65, rel=1e-2)
    assert results.iloc[1]["Mu"] == pytest.approx(20, rel=1e-5)
    assert results.iloc[1]["DCR"] == pytest.approx(0.573, rel=1e-5)


def test_check_flexure_ACI_318_19_2(slab_example_ACI_318_19_metric: OneWaySlab) -> None:
    # Testing the check of the reinforced slab with simple reinforcement
    # See calcpad: ACI 318-19 Slab Flexure 01 - Metric.cpd
    f = Forces(label="C1", M_y=50 * kNm)
    slab_example_ACI_318_19_metric.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=20 * cm)
    node = Node(section=slab_example_ACI_318_19_metric, forces=f)
    results = node.check_flexure()

    print(results)

    assert results.iloc[1]["Label"] == "Slab 03"
    assert results.iloc[1]["Comb."] == "C1"
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["As,min"] == pytest.approx(5.63, rel=1e-2)
    assert results.iloc[1]["As,req bot"] == pytest.approx(8.22, rel=1e-3)
    assert results.iloc[1]["As,req top"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["As"] == pytest.approx(5.65, rel=1e-2)
    assert results.iloc[1]["Mu"] == pytest.approx(50, rel=1e-5)
    assert results.iloc[1]["DCR"] == pytest.approx(1.431, rel=1e-5)


@pytest.mark.parametrize(
    ("height", "cover", "expected_capped"),
    [
        (7 * inch, 0.75 * inch, True),  # d = 6.00 in  -> uncapped 1.118
        (24 * inch, 1.5 * inch, False),  # d = 22.1 in  -> uncapped 0.672
    ],
)
def test_lambda_s_is_capped_at_one_imperial(height, cover, expected_capped) -> None:
    """
    ACI 318-19 Eq. 22.5.5.1.3 caps the size effect factor at 1.0.

    The factor exists to reduce V_c in deep members. Without the cap a shallow
    member (d < 10 in / 250 mm) gets lambda_s > 1, which inflates V_c instead —
    the opposite of the provision's purpose, and on the unsafe side. Slabs are
    the practical case: they are shallow AND have no minimum stirrups, which is
    the only branch of Table 22.5.5.1 where lambda_s applies.
    """
    slab = OneWaySlab(
        label="lambda_s",
        concrete=Concrete_ACI_318_19(name="C4", f_c=4 * ksi),
        steel_bar=SteelBar(name="fy60", f_y=60 * ksi),
        width=12 * inch,
        height=height,
        c_c=cover,
    )
    slab.set_slab_longitudinal_rebar_bot(d_b1=0.5 * inch, s_b1=10 * inch)
    Node(section=slab, forces=Forces(V_z=1.52 * kip)).check_shear()

    # d_shear and lambda_s have no equivalent on the public shear result yet.
    d_in = slab._d_shear.to("inch").magnitude
    uncapped = math.sqrt(2 / (1 + d_in / 10))
    assert (uncapped > 1.0) is expected_capped
    assert slab._lambda_s <= 1.0
    assert slab._lambda_s == pytest.approx(min(uncapped, 1.0), rel=1e-9)


def test_lambda_s_is_capped_at_one_metric() -> None:
    """Metric branch of the same cap: a 15 cm slab has d well under 250 mm."""
    slab = OneWaySlab(
        label="lambda_s_metric",
        concrete=Concrete_ACI_318_19(name="C25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=100 * cm,
        height=15 * cm,
        c_c=25 * mm,
    )
    slab.set_slab_longitudinal_rebar_bot(d_b1=10 * mm, s_b1=20 * cm)
    Node(section=slab, forces=Forces(V_z=30 * kN)).check_shear()

    # d_shear and lambda_s have no equivalent on the public shear result yet.
    d_mm = slab._d_shear.to("mm").magnitude
    assert d_mm < 250  # the regime where the cap bites
    assert math.sqrt(2 / (1 + 0.004 * d_mm)) > 1.0
    assert slab._lambda_s == pytest.approx(1.0, rel=1e-9)


def test_shear_check_with_no_tension_reinforcement_warns_and_does_not_raise() -> None:
    """
    A slab with no longitudinal steel has no shear capacity under ACI, and the
    check must say so instead of crashing.

    Slabs start with zero reinforcement on both faces (unlike beams, which fall
    back to starter bars). With no stirrups either, Table 22.5.5.1 gives
    V_c proportional to rho_w**(1/3) = 0, so phi*Vn is exactly zero. Dividing
    the demand by that used to raise ZeroDivisionError, breaking the contract
    that a check reports an insufficient section rather than blowing up.

    The UserWarning is the intended signal; the infinite DCR is the graceful
    way to report zero capacity.
    """
    slab = OneWaySlab(
        label="no-rebar",
        concrete=Concrete_ACI_318_19(name="C25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=100 * cm,
        height=20 * cm,
        c_c=25 * mm,
    )
    assert slab.reinforcement.bottom.A_s.magnitude == 0

    with pytest.warns(UserWarning, match="Longitudinal rebar As cannot be zero"):
        Node(section=slab, forces=Forces(V_z=5 * kN)).check_shear()

    assert slab.V_c.magnitude == 0
    assert slab.shear_design.DCR == float("inf")


def test_flexure_check_with_no_bottom_reinforcement_floors_phi_Mn() -> None:
    """
    Companion of the shear case on the flexure side.

    With no bottom steel and a non-negative moment, phi*Mn on the bottom face is
    zero. The check floors it at 0.01 kNm so the DCR stays a finite number
    instead of raising, and reports As = 0 so the cause is visible.
    """
    slab = OneWaySlab(
        label="no-rebar-flexure",
        concrete=Concrete_ACI_318_19(name="C25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=100 * cm,
        height=20 * cm,
        c_c=25 * mm,
    )
    results = Node(section=slab, forces=Forces(M_y=0 * kNm)).check_flexure()

    # phi*Mn has no equivalent on the public flexure result yet.
    assert slab._phi_M_n_bot.to("kN*m").magnitude == pytest.approx(0.01)
    assert slab.flexure_design.bottom.DCR == 0
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["As"] == 0


##########################################################
# SHEAR DESIGN: A SLAB MAY BE BUILT WITHOUT STIRRUPS
##########################################################


def _designed_slab(concrete: Concrete_ACI_318_19 | Concrete_EN_1992_2004, V_z: Quantity) -> OneWaySlab:
    """A 100x20 slab strip designed for one combination, bending plus shear."""
    slab = OneWaySlab(
        label="Slab shear",
        concrete=concrete,
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=100 * cm,
        height=20 * cm,
        c_c=20 * mm,
    )
    Node(section=slab, forces=Forces(label="C1", M_y=30 * kNm, V_z=V_z)).design()
    return slab


@pytest.mark.parametrize(
    "concrete",
    [
        Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        Concrete_EN_1992_2004(name="C25", f_c=25 * MPa),
    ],
    ids=["ACI_318_19", "EN_1992_2004"],
)
def test_slab_shear_design_places_no_stirrups_when_concrete_carries_the_shear(
    concrete: Concrete_ACI_318_19 | Concrete_EN_1992_2004,
) -> None:
    """ACI 318-19 7.6.3.1 and EN 1992-1-1 6.2.1(4) both let a slab go without.

    Designing for a moment and no shear used to hand the slab a full stirrup
    cage anyway, because the beam minimum was applied to it: A_v_min came out
    positive, the stirrup designer was always run, and the result was a
    reinforcement the demand never asked for and the code does not require.
    """
    slab = _designed_slab(concrete, V_z=0 * kN)

    transverse = slab.reinforcement.transverse
    assert transverse.n_stirrups == 0
    assert transverse.A_v.to("cm**2/m").magnitude == 0
    assert slab.shear_design.A_v_req.to("cm**2/m").magnitude == 0
    assert slab.shear_design.A_v_min.to("cm**2/m").magnitude == 0
    # The concrete alone still resists, so the section checks out.
    assert slab.shear_design.DCR == 0


@pytest.mark.parametrize(
    "concrete",
    [
        Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        Concrete_EN_1992_2004(name="C25", f_c=25 * MPa),
    ],
    ids=["ACI_318_19", "EN_1992_2004"],
)
def test_slab_shear_design_places_stirrups_once_the_demand_exceeds_the_concrete(
    concrete: Concrete_ACI_318_19 | Concrete_EN_1992_2004,
) -> None:
    """Waiving the minimum is not waiving the design: above V_c stirrups appear."""
    slab = _designed_slab(concrete, V_z=200 * kN)

    transverse = slab.reinforcement.transverse
    assert transverse.n_stirrups > 0
    assert transverse.A_v >= slab.shear_design.A_v_req
    assert slab.shear_design.DCR <= 1


# ---------------------------------------------------------------------------
# Shear design, where the demand passes what the concrete alone carries
# ---------------------------------------------------------------------------


def _slab_needing_stirrups(concrete_cls) -> OneWaySlab:
    """A 100 x 20 cm strip, which is too shallow for 100 kN of shear on its own."""
    return OneWaySlab(
        label="Slab shear",
        concrete=concrete_cls(name="C25", f_c=25 * MPa),
        steel_bar=SteelBar(name="B500S", f_y=420 * MPa),
        width=100 * cm,
        height=20 * cm,
        c_c=25 * mm,
    )


@pytest.mark.parametrize(
    ("concrete_cls", "d_b", "s_l", "s_w", "A_v"),
    [
        pytest.param(Concrete_ACI_318_19, 10 * mm, 7 * cm, 15 * cm, 74.80, id="ACI-318-19"),
        pytest.param(Concrete_EN_1992_2004, 6 * mm, 12 * cm, 12 * cm, 19.63, id="EN-1992-2004"),
    ],
)
def test_shear_design_of_a_slab_that_needs_stirrups(concrete_cls, d_b, s_l, s_w, A_v) -> None:
    """A strip carrying more shear than its concrete does gets a buildable grid of legs.

    The forces are the same in both codes; what differs is the detailing rules,
    so the expected bar and spacings are per code.

    The number this pins is the one the beam-shaped search got wrong. That
    search spreads an even number of legs between the outermost bar centres,
    which on a metre-wide strip starts at eight legs and never comes back down:
    it returned 78.54 cm2/m under ACI and 23.56 cm2/m under EN, and its ACI
    answer sat at 8 cm along the span against a 7.95 cm limit. Neither figure
    is anywhere near A_v_req, and neither can be -- see the second half of the
    test, which is the part worth reading.
    """
    # The premise: with no shear reinforcement the strip does not pass.
    bare = _slab_needing_stirrups(concrete_cls)
    bare.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=20 * cm)
    forces = Forces(label="C1", M_y=30 * kNm, V_z=100 * kN)
    Node(section=bare, forces=forces).check_shear()
    assert bare.shear_design.DCR > 1

    slab = _slab_needing_stirrups(concrete_cls)
    Node(section=slab, forces=Forces(label="C1", M_y=30 * kNm, V_z=100 * kN)).design()

    shear = slab.shear_design
    assert shear.d_b == d_b
    assert shear.s_l == s_l
    # The transverse spacing has no equivalent on the public shear result yet.
    assert slab._stirrup_s_w == s_w
    assert shear.A_v.to("cm**2/m").magnitude == pytest.approx(A_v, rel=1e-3)

    # The design is inside the limits it will be checked against. It was not
    # before: the limits are written on d, and d moves by the whole diameter of
    # a stirrup a slab did not have until the design assigned one.
    assert slab._stirrup_s_l <= slab._stirrup_s_max_l
    assert slab._stirrup_s_w <= slab._stirrup_s_max_w

    # And it covers the demand, with room to spare that no search can remove:
    # s_max_l is d/2 on a 20 cm slab, so the least steel the detailing rules
    # allow is already several times what the shear asks for.
    assert shear.A_v >= shear.A_v_req
    assert shear.DCR <= 1


def test_designed_slab_is_reinforced_in_both_directions() -> None:
    """The design speaks the slab's own parameterisation, not the beam's.

    ``design_shear`` used to hand the row to ``set_transverse_rebar``, which
    reads it as a number of closed stirrups. A slab has none: it is detailed by
    a bar diameter and a spacing each way, which is what
    ``set_slab_transverse_rebar`` takes.
    """
    slab = _slab_needing_stirrups(Concrete_ACI_318_19)
    Node(section=slab, forces=Forces(label="C1", M_y=30 * kNm, V_z=100 * kN)).design()

    A_db = math.pi * slab._stirrup_d_b**2 / 4
    n_legs = slab.width / slab._stirrup_s_w
    assert slab._A_v.to("cm**2/m").magnitude == pytest.approx(
        (A_db * n_legs / slab._stirrup_s_l).to("cm**2/m").magnitude, rel=1e-9
    )


def test_shear_design_of_an_imperial_slab_strip() -> None:
    """The imperial branch of the slab search, on its own whole-inch grid."""
    slab = OneWaySlab(
        label="Slab shear imperial",
        concrete=Concrete_ACI_318_19(name="C4", f_c=4 * ksi),
        steel_bar=SteelBar(name="G60", f_y=60 * ksi),
        width=12 * inch,
        height=10 * inch,
        c_c=0.75 * inch,
    )
    bare = OneWaySlab(
        label="Slab shear imperial",
        concrete=Concrete_ACI_318_19(name="C4", f_c=4 * ksi),
        steel_bar=SteelBar(name="G60", f_y=60 * ksi),
        width=12 * inch,
        height=10 * inch,
        c_c=0.75 * inch,
    )
    bare.set_slab_longitudinal_rebar_bot(d_b1=0.5 * inch, s_b1=8 * inch)
    Node(section=bare, forces=Forces(label="C1", V_z=12 * kip)).check_shear()
    assert bare.shear_design.DCR > 1

    Node(section=slab, forces=Forces(label="C1", M_y=10 * kNm, V_z=12 * kip)).design()

    shear = slab.shear_design
    assert shear.d_b == 0.375 * inch
    assert shear.s_l == 4 * inch
    assert slab._stirrup_s_w == 8 * inch
    assert slab._stirrup_s_l <= slab._stirrup_s_max_l
    assert slab._stirrup_s_w <= slab._stirrup_s_max_w
    assert shear.A_v >= shear.A_v_req
    assert shear.DCR <= 1


def test_transverse_rebar_set_on_a_slab_is_credited_by_the_shear_check() -> None:
    """Both codes read ``_stirrup_n`` to decide whether a section carries stirrups.

    ``set_slab_transverse_rebar`` left it at zero, so a slab given a full grid
    of legs was checked as a section with none: it reported phi*V_s = 0 and a
    zero stirrup diameter, and the report's minimum-diameter row had nothing to
    compare against.
    """
    slab = _slab_needing_stirrups(Concrete_ACI_318_19)
    slab.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=20 * cm)
    slab.set_slab_transverse_rebar(d_b=10 * mm, s_long=7 * cm, s_trans=15 * cm)

    transverse = slab.reinforcement.transverse
    assert transverse.d_b == 10 * mm
    assert transverse.n_stirrups > 0

    Node(section=slab, forces=Forces(label="C1", M_y=30 * kNm, V_z=100 * kN)).check_shear()
    assert slab._phi_V_s > 0 * kN


def test_clearing_slab_transverse_rebar_leaves_no_stirrup_behind() -> None:
    """The zero state has to clear the diameter and the count as well as A_v."""
    slab = _slab_needing_stirrups(Concrete_ACI_318_19)
    slab.set_slab_transverse_rebar(d_b=10 * mm, s_long=7 * cm, s_trans=15 * cm)

    slab.set_slab_transverse_rebar()

    transverse = slab.reinforcement.transverse
    assert transverse.n_stirrups == 0
    assert transverse.d_b.magnitude == 0
    assert transverse.A_v.to("cm**2/m").magnitude == 0
    assert slab._stirrup_s_w.magnitude == 0


def test_slab_shear_design_is_the_least_steel_the_spacing_limits_allow() -> None:
    """No bar, at any spacing the code allows, covers A_v_req with less steel.

    Written against the whole grid rather than against a number, because the
    point of the slab branch is that it finds the optimum, not that it happens
    to return 74.8 cm2/m. The grid starts at 1 cm, below the floor the search
    stops at: a tighter spacing only ever adds steel, so including it makes the
    claim stronger rather than unfair.
    """
    slab = _slab_needing_stirrups(Concrete_ACI_318_19)
    Node(section=slab, forces=Forces(label="C1", M_y=30 * kNm, V_z=100 * kN)).design()

    A_v_req = slab.shear_design.A_v_req
    chosen = slab.shear_design.A_v.to("cm**2/m").magnitude

    for _, row in slab.shear_design_results.iterrows():
        A_db = math.pi * row["d_b"] ** 2 / 4
        s_l_max = int(row["s_max_l"].to("cm").magnitude)
        s_w_max = int(min(row["s_max_w"], slab.width).to("cm").magnitude)
        for s_l in range(1, s_l_max + 1):
            for s_w in range(1, s_w_max + 1):
                A_v = A_db * (slab.width / (s_w * cm)) / (s_l * cm)
                if A_v >= A_v_req:
                    assert A_v.to("cm**2/m").magnitude >= chosen * (1 - 1e-9), (
                        f"Ø{row['d_b']} at {s_l} cm x {s_w} cm covers A_v_req with less steel"
                    )


##########################################################
# HOW A DESIGN IS WRITTEN ONTO A SLAB
##########################################################


def test_a_designed_slab_is_detailed_by_a_spacing_not_by_a_bar_count() -> None:
    """The search is the beam's; what it produces has to be read as a slab.

    Its answer is groups of bars -- ``n_1`` and ``n_2`` are one and the same
    layer -- and the beam applies it as such. A slab that went through it came
    out split into "2Ø12 + 4Ø12", two layers that are really one, with the
    spacing that actually details a slab left empty: there was no way to read
    the result other than by counting bars.
    """
    slab = _designed_slab(Concrete_ACI_318_19(name="H25", f_c=25 * MPa), V_z=0 * kN)

    bottom = slab.reinforcement.bottom
    assert len(bottom.layers) == 1
    layer = bottom.layers[0]
    assert layer.s is not None and layer.s.to("cm").magnitude > 0
    assert str(bottom) == f"Ø{layer.d_b:.4g~P}/{layer.s:.4g~P}"
    # Positions 2 and 4 stay empty: a slab layer is one bar at one spacing.
    assert slab._n2_b == 0 and slab._n4_b == 0
    # The spacing is the state, and the count is what it produces.
    assert _bars_at_spacing(slab._s_b1_b, slab.width) == layer.n


def test_the_spacing_a_design_leaves_is_one_that_can_be_drawn() -> None:
    """Rounded to the centimetre, and up, so the bars the search chose stay."""
    slab = _designed_slab(Concrete_ACI_318_19(name="H25", f_c=25 * MPa), V_z=0 * kN)

    spacing = slab._s_b1_b.to("cm").magnitude
    assert spacing == pytest.approx(round(spacing))
    n = slab.reinforcement.bottom.layers[0].n
    # Rounding up cannot cost a bar: n bars at width/n rounded up still ask for n.
    assert spacing >= slab.width.to("cm").magnitude / n
    assert slab.reinforcement.bottom.A_s >= slab.flexure_design.bottom.A_s_req


def test_a_spacing_is_never_rounded_into_fewer_bars_than_the_design_chose() -> None:
    """Rounding up is the rule, but it does not always hold.

    ``ceil(width / n)`` normally asks for the same ``n`` bars back. Once the
    bars are close enough that a whole centimetre spans more than one of them
    -- from about eleven in a metre -- rounding up drops one, so the spacing is
    checked against the count it produces and rounded down instead.
    """
    slab = OneWaySlab(
        label="Slab spacing",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=100 * cm,
        height=20 * cm,
        c_c=20 * mm,
    )

    # 100/6 = 16.7 -> 17 cm, and 17 cm still asks for 6 bars.
    assert slab._spacing_for_bars(6).to("cm").magnitude == 17
    # 100/11 = 9.1 -> 10 cm would only ask for 10, so it rounds the other way.
    assert slab._spacing_for_bars(11).to("cm").magnitude == 9
    for n in range(1, 30):
        assert _bars_at_spacing(slab._spacing_for_bars(n), slab.width) >= n
    assert slab._spacing_for_bars(0).magnitude == 0


def test_an_imperial_slab_is_designed_to_a_whole_inch_spacing() -> None:
    slab = OneWaySlab(
        label="Slab imperial flexure",
        concrete=Concrete_ACI_318_19(name="C4", f_c=4 * ksi),
        steel_bar=SteelBar(name="G60", f_y=60 * ksi),
        width=36 * inch,
        height=10 * inch,
        c_c=0.75 * inch,
    )
    Node(section=slab, forces=Forces(label="C1", M_y=20 * kNm)).design()

    spacing = slab._s_b1_b.to("inch").magnitude
    assert spacing == pytest.approx(round(spacing))
    assert _bars_at_spacing(slab._s_b1_b, slab.width) == slab._n1_b


def test_designing_a_slab_that_needs_no_top_steel_clears_the_top_spacing() -> None:
    """A cleared face has to lose its spacing, or the bars come back.

    The beam clears a face by zeroing the bar counts. On a slab those counts are
    derived from the spacing, so a spacing left behind would be turned back into
    bars by the next recalculation -- the following design of the bottom face.
    """
    slab = OneWaySlab(
        label="Slab top cleared",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=100 * cm,
        height=20 * cm,
        c_c=20 * mm,
    )
    slab.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=15 * cm)

    # Sagging only: nothing is required on the top face.
    Node(section=slab, forces=Forces(label="C1", M_y=30 * kNm)).design()

    assert slab._s_b1_t.magnitude == 0
    assert slab._n1_t == 0
    assert slab.reinforcement.top.layers == ()
    assert slab.reinforcement.top.A_s.magnitude == 0

    # The same, through the path the design takes when the face has no layout
    # at all to apply rather than an empty one.
    slab.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=15 * cm)
    assert slab.reinforcement.top.n_bars > 0

    slab._clear_top_longitudinal()

    assert slab._s_b1_t.magnitude == 0
    assert slab._s_b3_t.magnitude == 0
    assert slab.reinforcement.top.layers == ()
    # And the bars stay gone once another face is recalculated.
    slab.set_slab_longitudinal_rebar_bot(d_b1=10 * mm, s_b1=20 * cm)
    assert slab.reinforcement.top.layers == ()


def test_a_hogging_slab_is_designed_on_its_top_face_by_a_spacing() -> None:
    slab = OneWaySlab(
        label="Slab hogging",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=100 * cm,
        height=20 * cm,
        c_c=20 * mm,
    )
    Node(section=slab, forces=Forces(label="C1", M_y=-30 * kNm)).design()

    top = slab.reinforcement.top
    assert len(top.layers) == 1
    assert top.layers[0].s is not None
    assert _bars_at_spacing(slab._s_b1_t, slab.width) == top.layers[0].n
    assert top.A_s >= slab.flexure_design.top.A_s_req


@pytest.mark.parametrize(
    ("concrete", "height", "expected_cm"),
    [
        (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), 12 * cm, 36),  # 3h governs
        (Concrete_ACI_318_19(name="H25", f_c=25 * MPa), 25 * cm, 45),  # 450 mm governs
        (Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), 12 * cm, 36),  # 3h governs
        (Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), 25 * cm, 40),  # 400 mm governs
    ],
    ids=["ACI_3h", "ACI_450mm", "EN_3h", "EN_400mm"],
)
def test_the_code_caps_how_far_apart_the_bars_of_a_slab_may_sit(
    concrete: Concrete_ACI_318_19 | Concrete_EN_1992_2004, height: Quantity, expected_cm: float
) -> None:
    """ACI 318-19 7.7.2.3 is 3h or 450 mm; EN 1992-1-1 9.3.1.1(3) is 3h or 400 mm."""
    slab = OneWaySlab(
        label="Slab s_max",
        concrete=concrete,
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=100 * cm,
        height=height,
        c_c=20 * mm,
    )

    assert slab._max_bar_spacing().to("cm").magnitude == pytest.approx(expected_cm)


def test_a_design_is_never_spaced_beyond_the_code_maximum() -> None:
    """Area alone is not a layout.

    A light strip needs so little steel that the search covers it with the two
    bars it starts from, which on a metre of slab is a bar every half metre:
    the area is there and most of the slab is not reinforced. The spacing is
    capped at what the code allows, which only ever adds bars.
    """
    slab = OneWaySlab(
        label="Slab lightly loaded",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=100 * cm,
        height=12 * cm,
        c_c=20 * mm,
    )
    Node(section=slab, forces=Forces(label="C1", M_y=2 * kNm)).design()

    s_max = slab._max_bar_spacing()
    assert s_max.to("cm").magnitude == 36  # 3h on a 12 cm slab
    assert slab._s_b1_b <= s_max
    assert slab.reinforcement.bottom.layers[0].n == _bars_at_spacing(s_max, slab.width)
    assert slab.reinforcement.bottom.A_s >= slab.flexure_design.bottom.A_s_req


def test_a_slab_spaced_beyond_the_code_maximum_fails_the_check() -> None:
    """The bars add up to the area and the slab between them is still bare."""
    slab = OneWaySlab(
        label="Slab too spread out",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=100 * cm,
        height=25 * cm,
        c_c=25 * mm,
    )
    slab.set_slab_longitudinal_rebar_bot(d_b1=20 * mm, s_b1=60 * cm)

    Node(section=slab, forces=Forces(label="C1", M_y=20 * kNm)).check_flexure()

    rows = slab._data_min_max_flexure
    assert rows["Check"][3] == "Bar spacing bottom"
    assert rows["Value"][3] == pytest.approx(600)
    assert rows["Max."][3] == pytest.approx(450)
    assert rows["Ok?"][3] == "❌"
    assert slab._all_flexure_checks_passed is False


def test_a_beam_still_reports_the_clear_distance_between_its_bars() -> None:
    """The max-spacing row is the slab's; a beam has no such limit to report."""
    beam = RectangularBeam(
        label="B1",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=20 * cm,
        height=50 * cm,
        c_c=25 * mm,
    )
    Node(section=beam, forces=Forces(label="C1", M_y=100 * kNm)).design()

    rows = beam._data_min_max_flexure
    assert rows["Check"][3] == "Minimum spacing bottom"
    assert rows["Value"][3] == pytest.approx(beam._available_s_bot.to("mm").magnitude, rel=1e-9)
    assert rows["Max."][3] == ""


def test_the_flexure_report_of_a_slab_names_the_spacing() -> None:
    """The count row gives way to the notation the strip is drawn in."""
    slab = _designed_slab(Concrete_ACI_318_19(name="H25", f_c=25 * MPa), V_z=0 * kN)

    rows = slab._flexure_capacity_bot
    assert rows["Variable"][:2] == ["Ø1/s1", "Ø3/s3"]
    assert rows["Value"][0] == str(slab.reinforcement.bottom.layers[0])
    assert rows["Value"][1] == "-"
    # Every column of the table still has one entry per row.
    assert len({len(column) for column in rows.values()}) == 1


##########################################################
# HOW THE TRANSVERSE REINFORCEMENT IS LABELLED
##########################################################


def test_a_designed_slab_is_labelled_by_its_two_spacings() -> None:
    """A slab strip has no cage, so the beam's stirrup count says nothing.

    The design chooses a spacing each way and the label used to report only the
    longitudinal one, prefixed by a stirrup count that does not exist -- a
    100 cm strip came out as "7eO10 mm/3 cm" with the 7 cm across the width
    nowhere to be seen. The diameter is written once: both directions are the
    same bar.
    """
    slab = _designed_slab(Concrete_ACI_318_19(name="H25", f_c=25 * MPa), V_z=150 * kN)

    transverse = slab.reinforcement.transverse
    assert transverse.layout == "grid"
    assert transverse.s_w.to("cm").magnitude > 0
    label = str(transverse)
    assert label.startswith("Ø")
    assert "e Ø" not in label and "eØ" not in label
    assert label.count("Ø") == 1
    assert f"{transverse.s_l:.4g~P}" in label
    assert f"{transverse.s_w:.4g~P}" in label
    # The shear result carries the same notation as the section does.
    assert str(slab.shear_design) == label


def test_a_slab_without_stirrups_says_so() -> None:
    slab = _designed_slab(Concrete_ACI_318_19(name="H25", f_c=25 * MPa), V_z=0 * kN)

    assert str(slab.reinforcement.transverse) == "no stirrups"


def test_the_shear_report_of_a_slab_names_both_spacings() -> None:
    """The count row gives way to the spacing that actually detailed the strip."""
    slab = _designed_slab(Concrete_ACI_318_19(name="H25", f_c=25 * MPa), V_z=150 * kN)

    rows = slab._shear_reinforcement
    assert rows["Variable"][:3] == ["db", "sl", "sw"]
    assert rows["Shear reinforcement strength"][:3] == [
        "Stirrup diameter",
        "Stirrup spacing along length",
        "Stirrup spacing along width",
    ]
    assert rows["Unit"][:3] == ["mm", "cm", "cm"]
    assert rows["Value"][2] == pytest.approx(slab._stirrup_s_w.to("cm").magnitude, rel=1e-9)
    # Every column of the table still has one entry per row.
    assert len({len(column) for column in rows.values()}) == 1
