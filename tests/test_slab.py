import math

import pytest

from mento.node import Node
from mento.slab import OneWaySlab
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
