import pytest

from mento.node import Node
from mento.slab import OneWaySlab
from mento.material import Concrete_ACI_318_19, SteelBar, Concrete_EN_1992_2004
from mento.units import kip, inch, mm, cm, MPa, kNm, ksi
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


def test_shear_check_ACI_318_19_1(slab_example_ACI_318_19: OneWaySlab) -> None:
    # Example from Two-Way Flat Plate Concrete Floor System Analysis and Design (ACI 318-14) adjusted to ACI 318-19.
    # With guidance from CRSI Design Guide on ACI 318-19
    # See calcpad: ACI 318-19 Slab Shear 01 - Imperial.cpd
    f = Forces(V_z=1.52 * kip, N_x=0 * kip)
    node = Node(section=slab_example_ACI_318_19, forces=f)
    slab_example_ACI_318_19.set_slab_longitudinal_rebar_bot(d_b1=0.5 * inch, s_b1=10 * inch)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["ØVc"] == pytest.approx(23.92, rel=1e-3)
    assert results.iloc[1]["ØVs"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["ØVn"] == pytest.approx(23.92, rel=1e-3)
    assert results.iloc[1]["ØVmax"] == pytest.approx(145.45, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.283, rel=1e-3)

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
