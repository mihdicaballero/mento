import pytest
import numpy as np

from mento.node import Node
from mento.slab import OneWaySlab
from mento.material import Concrete_ACI_318_19, SteelBar, Concrete_EN_1992_2004
from mento.units import kip, inch, mm, cm, MPa, ft, kNm, ksi
from mento.forces import Forces


@pytest.fixture()
def slab_example_EN_1992_2004() -> OneWaySlab:
    concrete = Concrete_EN_1992_2004(name="C25", f_ck=25 * MPa)
    steelBar = SteelBar(name="B500S", f_y=500 * MPa)
    custom_settings = {"clear_cover": 2.5 * cm}
    slab = OneWaySlab(
        label="Slab 01",
        concrete=concrete,
        steel_bar=steelBar,
        thickness=20 * cm,
        width=100 * cm,
        settings=custom_settings,
    )
    return slab


@pytest.fixture()
def slab_example_ACI_318_19() -> OneWaySlab:
    concrete = Concrete_ACI_318_19(name="H25", f_c=4 * ksi)
    steelBar = SteelBar(name="ADN 420", f_y=60 * ksi)
    custom_settings = {"clear_cover": 0.75 * inch}
    slab = OneWaySlab(
        label="Slab 01",
        concrete=concrete,
        steel_bar=steelBar,
        thickness=7 * inch,
        width=12 * inch,
        settings=custom_settings,
    )
    return slab


def test_shear_check_ACI_318_19_1(beam_example_imperial: OneWaySlab) -> None:
    # Example from Two-Way Flat Plate Concrete Floor System Analysis and Design (ACI 318-14) adjusted to ACI 318-19.
    # With guidance from CRSI Design Guide on ACI 318-19
    # See calcpad: ACI 318-19 Slab Shear 01 - Imperial.cpd
    f = Forces(V_z=1.52 * kip, N_x=0 * kip)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["ØVc"].magnitude == pytest.approx(58.288, rel=1e-3)
    assert results.iloc[0]["ØVs"].magnitude == pytest.approx(180.956, rel=1e-3)
    assert results.iloc[0]["ØVn"].magnitude == pytest.approx(239.247, rel=1e-3)
    assert results.iloc[0]["ØVmax"].magnitude == pytest.approx(291.44, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(0.70144, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["Vu<ØVmax"] is np.True_
    assert results.iloc[0]["Vu<ØVn"] is np.True_


def test_check_flexure_ACI_318_19_1(slab_example_ACI_318_19: OneWaySlab) -> None:
    # Testing the check of the reinforced slab with simple reinforcement
    # # See calcpad: ACI 318-19 Slab Flexure 01 - Metric.cpd
    f = Forces(label="C1", M_y=20 * kNm)
    slab_example_ACI_318_19.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=20 * cm)
    node = Node(section=slab_example_ACI_318_19, forces=f)
    results = node.check_flexure()

    assert results.iloc[0]["Section Label"] == "Slab 01"
    assert results.iloc[0]["Load Combo"] == "C1"
    assert results.iloc[0]["Position"] == "Bottom"
    assert results.iloc[0]["As,min"].to(cm**2).magnitude == pytest.approx(
        5.26, rel=1e-2
    )
    assert results.iloc[0]["As,req bot"].to(cm**2).magnitude == pytest.approx(
        33.1712108, rel=1e-3
    )
    assert results.iloc[0]["As,req top"].to(cm**2).magnitude == pytest.approx(
        4.9250577, rel=1e-3
    )
    assert results.iloc[0]["As"].to(inch**2).magnitude == pytest.approx(
        5.66, rel=1e-2
    )
    assert results.iloc[0]["Mu"].to(kip * ft).magnitude == pytest.approx(400, rel=1e-5)
