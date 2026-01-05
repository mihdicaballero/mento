import pytest
import numpy as np
import matplotlib
import pandas as pd
from unittest.mock import patch

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Rectangle, FancyBboxPatch

from mento.node import Node
from mento.beam import RectangularBeam
from mento.material import (
    Concrete_ACI_318_19,
    SteelBar,
    Concrete_EN_1992_2004,
    Concrete_CIRSOC_201_25,
)
from mento.units import psi, kip, inch, ksi, mm, kN, cm, MPa, ft, kNm
from mento.forces import Forces
from mento.codes.ACI_318_19_beam import (
    _determine_nominal_moment_simple_reinf_ACI_318_19,
    _determine_nominal_moment_double_reinf_ACI_318_19,
)
from mento.results import CUSTOM_COLORS
from mento.settings import BeamSettings
from mento.rebar import Rebar


@pytest.fixture()
def beam_example_imperial() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steelBar = SteelBar(name="ADN 420", f_y=60 * ksi)
    section = RectangularBeam(
        label="V101",
        concrete=concrete,
        steel_bar=steelBar,
        width=10 * inch,
        height=16 * inch,
        c_c=1.5 * inch,
    )
    return section


@pytest.fixture()
def beam_example_EN_1992_2004_01() -> RectangularBeam:
    # Example from Calcpad EN 1992-1-1_2004 Beam Flexure 01 - Metric v2
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    steelBar = SteelBar(name="B500S", f_y=500 * MPa)
    section = RectangularBeam(
        label="B_Example_EN_01",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=60 * cm,
        c_c=2.6 * cm,
    )
    return section


@pytest.fixture()
def beam_example_EN_1992_2004_02() -> RectangularBeam:
    concrete = Concrete_EN_1992_2004(name="C35", f_c=35 * MPa)
    steelBar = SteelBar(name="B500S", f_y=500 * MPa)
    section = RectangularBeam(
        label="B_Example_EN_02",
        concrete=concrete,
        steel_bar=steelBar,
        width=30 * cm,
        height=40 * cm,
        c_c=3.8 * cm,
    )
    return section


@pytest.fixture()
def beam_example_EN_1992_2004_03() -> RectangularBeam:
    concrete = Concrete_EN_1992_2004(name="C60", f_c=60 * MPa)
    steelBar = SteelBar(name="B400S", f_y=400 * MPa)
    section = RectangularBeam(
        label="B_Example_EN_03",
        concrete=concrete,
        steel_bar=steelBar,
        width=30 * cm,
        height=50 * cm,
        c_c=3 * cm,
    )
    return section


@pytest.fixture()
def beam_example_CIRSOC_201_2025() -> RectangularBeam:
    concrete = Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa)
    steel_bar = SteelBar(name="ADN 420", f_y=420 * MPa)
    section = RectangularBeam(
        label="Test",
        concrete=concrete,
        steel_bar=steel_bar,
        width=20 * cm,
        height=60 * cm,
        c_c=2.5 * cm,
    )
    return section


def build_metric_beam(settings: BeamSettings | None = None) -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    return RectangularBeam(
        label="additional",
        concrete=concrete,
        steel_bar=steel,
        width=30 * cm,
        height=50 * cm,
        c_c=40 * mm,
        settings=settings,
    )


def test_beam_respects_provided_settings_unit_system() -> None:
    settings = BeamSettings(unit_system="imperial", max_bars_per_layer=3)
    beam = build_metric_beam(settings=settings)

    assert beam.settings.unit_system == "metric"
    assert beam.settings.max_bars_per_layer == 3
    assert beam.mode == "beam"


def test_rebar_designer_factory_returns_rebar() -> None:
    beam = build_metric_beam()

    assert isinstance(beam._create_rebar_designer(), Rebar)


def test_initialize_longitudinal_rebar_attributes_metric_defaults() -> None:
    beam = build_metric_beam()

    assert beam._n1_b == 2
    assert beam._n1_t == 2
    assert beam._d_b1_b.magnitude == pytest.approx(8)
    assert beam._d_b1_t.magnitude == pytest.approx(8)


def test_shear_check_EN_1992_2004_rebar_1(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    f = Forces(V_z=100 * kN)
    beam_example_EN_1992_2004_01.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=25 * cm)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(1.83, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(2.262, rel=1e-3)
    assert results.iloc[1]["VEd,1"] == pytest.approx(100, rel=1e-3)
    assert results.iloc[1]["VEd,2"] == pytest.approx(100, rel=1e-3)
    assert results.iloc[1]["VRd,c"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VRd,s"] == pytest.approx(123.924, rel=1e-3)
    assert results.iloc[1]["VRd"] == pytest.approx(123.924, rel=1e-3)
    assert results.iloc[1]["VRd,max"] == pytest.approx(312.811, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.8069, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["VEd,1≤VRd,max"] is True
    assert results.iloc[1]["VEd,2≤VRd"] is True


def test_shear_check_EN_1992_2004_rebar_2(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    f = Forces(V_z=350 * kN)
    beam_example_EN_1992_2004_01.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=25 * cm)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(7.533, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(2.262, rel=1e-3)
    assert results.iloc[1]["VEd,1"] == pytest.approx(350, rel=1e-3)
    assert results.iloc[1]["VEd,2"] == pytest.approx(350, rel=1e-3)
    assert results.iloc[1]["VRd,c"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VRd,s"] == pytest.approx(105.099, rel=1e-3)
    assert results.iloc[1]["VRd"] == pytest.approx(105.099, rel=1e-3)
    assert results.iloc[1]["VRd,max"] == pytest.approx(350, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(3.33, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["VEd,1≤VRd,max"] is True
    assert results.iloc[1]["VEd,2≤VRd"] is False


def test_shear_check_EN_1992_2004_rebar_3(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    f = Forces(V_z=500 * kN)
    beam_example_EN_1992_2004_01.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=25 * cm)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(22.817, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(2.26, rel=1e-3)
    assert results.iloc[1]["VEd,1"] == pytest.approx(500, rel=1e-3)
    assert results.iloc[1]["VEd,2"] == pytest.approx(500, rel=1e-3)
    assert results.iloc[1]["VRd,c"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VRd,s"] == pytest.approx(49.566, rel=1e-3)
    assert results.iloc[1]["VRd"] == pytest.approx(49.566, rel=1e-3)
    assert results.iloc[1]["VRd,max"] == pytest.approx(453.6, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(10.088, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["VEd,1≤VRd,max"] is False
    assert results.iloc[1]["VEd,2≤VRd"] is False


def test_shear_check_EN_1992_2004_no_rebar_1(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    # Example from "EN 1992-1-1_2004 Beam Shear 01 - Metric.cpd"
    f = Forces(V_z=30 * kN)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VEd,1"] == pytest.approx(30, rel=1e-3)
    assert results.iloc[1]["VEd,2"] == pytest.approx(30, rel=1e-3)
    assert beam_example_EN_1992_2004_01._stirrup_d_b.to("mm").magnitude == 0
    assert beam_example_EN_1992_2004_01._d_shear.to("cm").magnitude == pytest.approx(56.6, rel=1e-3)
    assert results.iloc[1]["VRd,c"] == pytest.approx(56.51, rel=1e-3)
    assert results.iloc[1]["VRd,s"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VRd"] == pytest.approx(56.51, rel=1e-3)
    assert results.iloc[1]["VRd,max"] == pytest.approx(56.51, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.531, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["VEd,1≤VRd,max"] is True
    assert results.iloc[1]["VEd,2≤VRd"] is True


def test_shear_check_EN_1992_2004_no_rebar_2(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    f = Forces(V_z=30 * kN)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VEd,1"] == pytest.approx(30, rel=1e-3)
    assert results.iloc[1]["VEd,2"] == pytest.approx(30, rel=1e-3)
    assert beam_example_EN_1992_2004_01._stirrup_d_b.to("mm").magnitude == 0
    assert beam_example_EN_1992_2004_01._d_shear.to("cm").magnitude == pytest.approx(57, rel=1e-3)
    assert results.iloc[1]["VRd,c"] == pytest.approx(40.09, rel=1e-3)
    assert results.iloc[1]["VRd,s"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VRd"] == pytest.approx(40.09, rel=1e-3)
    assert results.iloc[1]["VRd,max"] == pytest.approx(40.09, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.748, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["VEd,1≤VRd,max"] is True
    assert results.iloc[1]["VEd,2≤VRd"] is True


def test_shear_check_EN_1992_2004_no_rebar_3(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    f = Forces(N_x=50 * kN, V_z=30 * kN)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VEd,1"] == pytest.approx(30, rel=1e-3)
    assert results.iloc[1]["VEd,2"] == pytest.approx(30, rel=1e-3)
    assert beam_example_EN_1992_2004_01._stirrup_d_b.to("mm").magnitude == 0
    assert beam_example_EN_1992_2004_01._d_shear.to("cm").magnitude == pytest.approx(56.6, rel=1e-3)
    assert results.iloc[1]["VRd,c"] == pytest.approx(63.59, rel=1e-3)
    assert results.iloc[1]["VRd,s"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VRd"] == pytest.approx(63.59, rel=1e-3)
    assert results.iloc[1]["VRd,max"] == pytest.approx(63.59, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.472, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["VEd,1≤VRd,max"] is True
    assert results.iloc[1]["VEd,2≤VRd"] is True


def test_shear_design_EN_1992_2004_1(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    f = Forces(N_x=0 * kN, V_z=30 * kN)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.design_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(1.62, rel=1e-3)
    assert results.iloc[1]["VEd,1"] == pytest.approx(30, rel=1e-3)
    assert results.iloc[1]["VEd,2"] == pytest.approx(30, rel=1e-3)
    assert results.iloc[1]["VRd,c"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VRd,s"] == pytest.approx(88.52, rel=1e-3)
    assert results.iloc[1]["VRd"] == pytest.approx(88.52, rel=1e-3)
    assert results.iloc[1]["VRd,max"] == pytest.approx(312.81, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.339, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["VEd,1≤VRd,max"] is True
    assert results.iloc[1]["VEd,2≤VRd"] is True


def test_shear_check_ACI_318_19_1(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=37.727 * kip, N_x=0 * kip)
    beam_example_imperial.set_transverse_rebar(n_stirrups=1, d_b=0.5 * inch, s_l=6 * inch)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(2.12, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(10.0623, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(16.624, rel=1e-3)
    assert results.iloc[1]["ØVc"] == pytest.approx(58.288, rel=1e-3)
    assert results.iloc[1]["ØVs"] == pytest.approx(180.956, rel=1e-3)
    assert results.iloc[1]["ØVn"] == pytest.approx(239.247, rel=1e-3)
    assert results.iloc[1]["ØVmax"] == pytest.approx(291.44, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.70144, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["Vu≤ØVmax"] is True
    assert results.iloc[1]["Vu≤ØVn"] is True


def test_shear_check_ACI_318_19_2(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=37.727 * kip, N_x=20 * kip)
    beam_example_imperial.set_transverse_rebar(n_stirrups=1, d_b=0.5 * inch, s_l=6 * inch)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(2.12, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(9.1803, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(16.624, rel=1e-3)
    assert results.iloc[1]["ØVc"] == pytest.approx(67.888, rel=1e-3)
    assert results.iloc[1]["ØVs"] == pytest.approx(180.959, rel=1e-3)
    assert results.iloc[1]["ØVn"] == pytest.approx(248.847, rel=1e-3)
    assert results.iloc[1]["ØVmax"] == pytest.approx(301.041, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.6743, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["Vu≤ØVmax"] is True
    assert results.iloc[1]["Vu≤ØVn"] is True


def test_shear_check_ACI_318_19_no_rebar_1(
    beam_example_imperial: RectangularBeam,
) -> None:
    # Tested with "ACI 318-19 Beam Shear 01 - Imperial.cpd" for beam that needs rebar
    f = Forces(V_z=8 * kip, N_x=0 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(2.12, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(2.12, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert beam_example_imperial._stirrup_d_b.to("mm").magnitude == 0
    assert beam_example_imperial._d_shear.to("cm").magnitude == pytest.approx(36.04, rel=1e-3)
    assert results.iloc[1]["ØVc"] == pytest.approx(35.48, rel=1e-3)
    assert results.iloc[1]["ØVs"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["ØVn"] == pytest.approx(35.48, rel=1e-3)
    assert results.iloc[1]["ØVmax"] == pytest.approx(274.96, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(1.003, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["Vu≤ØVmax"] is True
    assert results.iloc[1]["Vu≤ØVn"] is False


def test_shear_check_ACI_318_19_no_rebar_2(
    beam_example_imperial: RectangularBeam,
) -> None:
    # Tested with "ACI 318-19 Beam Shear 01 - Imperial.cpd" for beem that doesn't need rebar
    f = Forces(V_z=6 * kip, N_x=0 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert beam_example_imperial._d_shear.to("cm").magnitude == pytest.approx(36.04, rel=1e-3)
    assert results.iloc[1]["Av,min"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert beam_example_imperial._k_c_min.to("MPa").magnitude == pytest.approx(0.517, rel=1e-3)
    assert results.iloc[1]["ØVc"] == pytest.approx(35.48, rel=1e-3)
    assert results.iloc[1]["ØVs"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["ØVn"] == pytest.approx(35.48, rel=1e-3)
    assert results.iloc[1]["ØVmax"] == pytest.approx(274.96, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.752, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["Vu≤ØVmax"] is True
    assert results.iloc[1]["Vu≤ØVn"] is True


def test_shear_design_ACI_318_19(beam_example_imperial: RectangularBeam) -> None:
    # Tested with "ACI 318-19 Beam Shear 01 - Imperial.cpd" for beem that needs rebar
    f = Forces(V_z=37.727 * kip, N_x=0 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.design_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert beam_example_imperial._d_shear.to("cm").magnitude == pytest.approx(35.08, rel=1e-3)
    assert results.iloc[1]["Av,min"] == pytest.approx(2.12, rel=1e-3)
    assert beam_example_imperial._V_s_req.to("kN").magnitude == pytest.approx(109.53, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(10.06, rel=1e-3)
    assert results.iloc[1]["ØVc"] == pytest.approx(58.29, rel=1e-3)
    assert results.iloc[1]["ØVs"] == pytest.approx(122.15, rel=1e-3)
    assert results.iloc[1]["ØVn"] == pytest.approx(180.44, rel=1e-3)
    assert results.iloc[1]["ØVmax"] == pytest.approx(291.44, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(11.22, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.93, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["Vu≤ØVmax"] is True
    assert results.iloc[1]["Vu≤ØVn"] is True

    # Design result
    assert beam_example_imperial._stirrup_n == 1
    assert beam_example_imperial._stirrup_d_b.to("mm").magnitude == pytest.approx(9.525, rel=1e-3)
    assert beam_example_imperial._stirrup_s_l.to("cm").magnitude == pytest.approx(12.7, rel=1e-3)


def test_shear_design_CIRSOC_201_2025(
    beam_example_CIRSOC_201_2025: RectangularBeam,
) -> None:
    f = Forces(label="ELU_01", V_z=120 * kN, N_x=0 * kN)
    node = Node(section=beam_example_CIRSOC_201_2025, forces=f)
    results = node.design_shear()

    assert results.iloc[1]["Av,min"] == pytest.approx(1.67, rel=1e-2)
    assert results.iloc[1]["Av,req"] == pytest.approx(2.69, rel=1e-2)
    assert results.iloc[1]["Av"] == pytest.approx(4.04, rel=1e-2)
    assert results.iloc[1]["ØVc"] == pytest.approx(72.04, rel=1e-2)
    assert results.iloc[1]["ØVs"] == pytest.approx(71.89, rel=1e-2)
    assert results.iloc[1]["ØVn"] == pytest.approx(143.92, rel=1e-2)
    assert results.iloc[1]["ØVmax"] == pytest.approx(351.71, rel=1e-2)
    assert results.iloc[1]["DCR"] == pytest.approx(0.834, rel=1e-2)

    assert results.iloc[1]["Vu≤ØVmax"] is True
    assert results.iloc[1]["Vu≤ØVn"] is True

    # Design result
    assert beam_example_CIRSOC_201_2025._stirrup_n == 1
    assert beam_example_CIRSOC_201_2025._stirrup_d_b.to("mm").magnitude == 6
    assert beam_example_CIRSOC_201_2025._stirrup_s_l.to("cm").magnitude == 14


# # ------- FLEXURE TEST --------------


def test_flexure_check_EN_1992_2004_01(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    # Example from Calcpad EN 1992-1-1_2004 Beam Flexure 01 - Metric v2
    f = Forces(M_y=150 * kNm)
    beam_example_EN_1992_2004_01.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=15 * cm)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_top(n1=0, d_b1=0 * mm)
    assert beam_example_EN_1992_2004_01._A_s_bot.to(cm**2).magnitude == pytest.approx(8.042, rel=1e-2)
    assert beam_example_EN_1992_2004_01._d_bot.to(cm).magnitude == pytest.approx(56.0, rel=1e-2)
    assert beam_example_EN_1992_2004_01.width.to(cm).magnitude == pytest.approx(20.0, rel=1e-2)
    assert beam_example_EN_1992_2004_01._stirrup_d_b.to(mm).magnitude == pytest.approx(6.0, rel=1e-2)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_flexure()
    assert results.iloc[1]["Label"] == "B_Example_EN_01"
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["As,min"] == pytest.approx(1.49, rel=1e-2)
    assert results.iloc[1]["As,req bot"] == pytest.approx(6.656, rel=1e-3)
    assert results.iloc[1]["As,req top"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["As"] == pytest.approx(8.042, rel=1e-2)
    assert results.iloc[1]["MRd"] == pytest.approx(178.56, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.8400771633, rel=1e-3)


def test_flexure_check_EN_1992_2004_02(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    # Example from Calcpad Validation
    f = Forces(M_y=450 * kNm)
    beam_example_EN_1992_2004_01.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=15 * cm)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=5, d_b1=25 * mm)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_top(n1=5, d_b1=12 * mm)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_flexure()
    assert results.iloc[1]["Label"] == "B_Example_EN_01"
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["As,min"] == pytest.approx(1.482, rel=1e-2)
    assert results.iloc[1]["As,req bot"] == pytest.approx(20.77011560316, rel=1e-4)
    assert results.iloc[1]["As,req top"] == pytest.approx(9.34, rel=1e-4)
    assert results.iloc[1]["As"] == pytest.approx(24.54, rel=1e-4)
    assert results.iloc[1]["MRd"] == pytest.approx(385.156, rel=1e-2)
    assert results.iloc[1]["DCR"] == pytest.approx(1.168, rel=1e-3)


def test_flexure_check_EN_1992_2004_03(
    beam_example_EN_1992_2004_02: RectangularBeam,
) -> None:
    f = Forces(M_y=120 * kNm)
    beam_example_EN_1992_2004_02.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=15 * cm)
    beam_example_EN_1992_2004_02.set_longitudinal_rebar_bot(n1=5, d_b1=25 * mm)
    beam_example_EN_1992_2004_02.set_longitudinal_rebar_top(n1=5, d_b1=12 * mm)
    node = Node(section=beam_example_EN_1992_2004_02, forces=f)
    results = node.check_flexure()
    assert results.iloc[1]["Label"] == "B_Example_EN_02"
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["As,min"] == pytest.approx(1.72, rel=1e-3)
    assert results.iloc[1]["As,req bot"] == pytest.approx(8.69106875642, rel=1e-3)
    assert results.iloc[1]["As,req top"] == pytest.approx(0, rel=1e-3)


def test_flexure_check_EN_1992_2004_04(
    beam_example_EN_1992_2004_03: RectangularBeam,
) -> None:
    # Example from Lecture-3-Bending-and-Shear-in-Beams-Concrete Centre - Page 14
    f = Forces(M_y=-370 * kNm)
    beam_example_EN_1992_2004_03.set_longitudinal_rebar_top(n1=6, d_b1=25 * mm)
    beam_example_EN_1992_2004_03.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=20 * cm)
    assert beam_example_EN_1992_2004_03._d_top == 45.15 * cm
    node = Node(section=beam_example_EN_1992_2004_03, forces=f)
    results = node.check_flexure()
    assert results.iloc[1]["Label"] == "B_Example_EN_03"
    assert results.iloc[1]["Position"] == "Top"
    assert results.iloc[1]["As,min"] == pytest.approx(4.048, rel=1e-3)
    assert results.iloc[1]["As,req top"] == pytest.approx(25.63866, rel=1e-3)


def test_en1992_high_strength_concrete_compression_reinforcement() -> None:
    """Test EN 1992 beam with f_ck > 50 MPa requiring compression reinforcement."""
    # Use high-strength concrete (f_ck = 60 MPa)
    concrete_en = Concrete_EN_1992_2004(name="C60/75", f_c=60 * MPa)
    steel = SteelBar(name="B500", f_y=500 * MPa)

    beam = RectangularBeam(
        label="TEST_EN_HIGH",
        concrete=concrete_en,
        steel_bar=steel,
        width=20 * cm,
        height=60 * cm,
        c_c=2.5 * cm,
    )

    # Apply large moment to require compression reinforcement
    # This should exceed M_lim and trigger doubly reinforced section
    forces = Forces(M_y=500 * kNm)
    node = Node(beam, forces)

    # Design the beam
    result = node.design_flexure()

    # Should complete without error
    assert result is not None
    assert isinstance(result, pd.DataFrame)

    # Beam should have compression reinforcement (top rebar)
    assert beam._A_s_top > 0 * cm**2

    # Beam should be marked as doubly reinforced
    assert beam._doubly_reinforced is True

    # Check that x_u calculation used k_3 and k_4 (high-strength path)
    # This is implicit in the calculation, but we can verify the section was designed
    check_result = node.check_flexure()
    assert isinstance(check_result, pd.DataFrame)


def test_en1992_normal_strength_concrete_compression_reinforcement() -> None:
    """Test EN 1992 beam with f_ck <= 50 MPa requiring compression reinforcement."""
    # Use normal-strength concrete (f_ck = 30 MPa)
    concrete_en = Concrete_EN_1992_2004(name="C30/37", f_c=30 * MPa)
    steel = SteelBar(name="B500", f_y=500 * MPa)

    beam = RectangularBeam(
        label="TEST_EN_NORMAL",
        concrete=concrete_en,
        steel_bar=steel,
        width=25 * cm,
        height=50 * cm,
        c_c=3 * cm,
    )

    # Apply large moment to require compression reinforcement
    forces = Forces(M_y=200 * kNm, N_x=0 * kN, V_z=0 * kN)
    node = Node(beam, forces)

    # Design the beam
    result = node.design_flexure()

    # Should complete without error
    assert result is not None
    assert isinstance(result, pd.DataFrame)

    # If doubly reinforced, should have top rebar
    if beam._doubly_reinforced:
        assert beam._A_s_top > 0 * cm**2

    # Check that x_u calculation used k_1 and k_2 (normal-strength path)
    check_result = node.check_flexure()
    assert isinstance(check_result, pd.DataFrame)


def test_en1992_negative_moment() -> None:
    """Test EN 1992 beam with f_ck <= 50 MPa requiring compression reinforcement."""
    # Use normal-strength concrete (f_ck = 30 MPa)
    concrete_en = Concrete_EN_1992_2004(name="C30/37", f_c=30 * MPa)
    steel = SteelBar(name="B500", f_y=500 * MPa)

    beam = RectangularBeam(
        label="TEST_EN_NORMAL",
        concrete=concrete_en,
        steel_bar=steel,
        width=25 * cm,
        height=50 * cm,
        c_c=3 * cm,
    )

    # Apply large moment to require compression reinforcement
    forces = Forces(M_y=-50 * kNm, N_x=0 * kN, V_z=0 * kN)
    node = Node(beam, forces)

    # Design the beam
    result = node.design_flexure()

    # Should complete without error
    assert result is not None
    assert isinstance(result, pd.DataFrame)


def test_flexural_beam_determine_nominal_moment_simple_reinf_ACI_318_19() -> None:
    """
    Test para la funcion que determina el momento nominal de una seccion armada solo con refuerzo de traccion
    ver el ejemplo en:
    Example from https://www.google.com/search?q=ACI+318+19+nominal+moment+of+section&oq=ACI+318+19+nominal+moment+of+section&gs_lcrp=EgZjaHJvbWUyBggAEEUYOTIHCAEQIRigAdIBCTIyNzkzajBqN6gCALACAA&sourceid=chrome&ie=UTF-8#fpstate=ive&vld=cid:dab5f5e7,vid:m4H0QbGDYIg,st:0
    """
    concrete = Concrete_ACI_318_19(name="C6", f_c=6000 * psi)
    steelBar = SteelBar(name="G80", f_y=80000 * psi)
    beam = RectangularBeam(
        label="B20x30",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * inch,
        height=30 * inch,
        c_c=1 * inch,
    )
    A_s = 10.92 * inch**2
    d = 27 * inch
    result = _determine_nominal_moment_simple_reinf_ACI_318_19(beam, A_s, d)

    # Compare dictionaries with a tolerance for floating-point values, in inch**2
    # Result: 1653.84kip.ft
    assert result.to(kip * ft).magnitude == pytest.approx(1653.84, rel=1e-2)


def test_flexural_beam_determine_nominal_moment_double_reinf_ACI_318_19() -> None:
    """
    Test para la funcion que determina el momento nominal de una seccion doblemente armada
    Example from https://www.youtube.com/watch?v=7BOuV1gCcgc
    """
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steelBar = SteelBar(name="G60", f_y=60000 * psi)
    beam = RectangularBeam(
        label="B14x27",
        concrete=concrete,
        steel_bar=steelBar,
        width=14 * inch,
        height=27 * inch,
        c_c=1 * inch,
    )
    A_s = 6 * inch**2
    d = 24 * inch
    d_prime = 2.5 * inch
    A_s_prime = 1.8 * inch**2
    result = _determine_nominal_moment_double_reinf_ACI_318_19(beam, A_s, d, d_prime, A_s_prime)
    assert result.to(kip * ft).magnitude == pytest.approx(639.12, rel=1e-2)


# BEAM FROM CALPCAD "ACI 318-19 Beam Flexure 02-TEST_1"
@pytest.fixture()
def beam_example_flexure_ACI() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="fc 4000", f_c=4000 * psi)
    steelBar = SteelBar(name="fy 60000", f_y=60 * ksi)
    section = RectangularBeam(
        concrete=concrete,
        steel_bar=steelBar,
        width=12 * inch,
        height=24 * inch,
        c_c=1.5 * inch,
        label="B-12x24",
    )
    return section


def test_check_flexure_ACI_318_19_1(beam_example_flexure_ACI: RectangularBeam) -> None:
    # Testing the check of the reinforced beam with a large moment that requires compression reinforcement, the moment being positive
    # # See calcpad: ACI 318-19 Beam Flexure 02-TEST_1.cpd
    f = Forces(label="Test_01", M_y=400 * kip * ft)
    beam_example_flexure_ACI.set_longitudinal_rebar_bot(n1=2, d_b1=1.41 * inch, n3=2, d_b3=1.27 * inch)
    beam_example_flexure_ACI.set_longitudinal_rebar_top(n1=2, d_b1=0.75 * inch)
    node = Node(section=beam_example_flexure_ACI, forces=f)
    results = node.check_flexure()
    assert beam_example_flexure_ACI._d_bot.to("inch").magnitude == pytest.approx(20.3719, rel=1e-3)
    assert results.iloc[1]["Label"] == "B-12x24"
    assert results.iloc[1]["Comb."] == "Test_01"
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["As,min"] == pytest.approx(5.257, rel=1e-2)
    assert results.iloc[1]["As,req bot"] == pytest.approx(33.17, rel=1e-3)
    assert results.iloc[1]["As,req top"] == pytest.approx(4.93, rel=1e-3)
    assert results.iloc[1]["As"] == pytest.approx(36.49, rel=1e-2)
    assert results.iloc[1]["Mu"] == pytest.approx(542.33, rel=1e-3)
    assert results.iloc[1]["ØMn"] == pytest.approx(588.73, rel=1e-3)


def test_check_flexure_ACI_318_19_2(beam_example_flexure_ACI: RectangularBeam) -> None:
    # Testeo el checkeo de la viga armada con un momento grande que requiere armadura de compresión, siendo el momento NEGATIVO.
    # Ver calcpad: ACI 318-19 Beam Flexure 02-TEST_1.cpd
    # Invierto el armado respecto al test anterior para poder usar el mismo calcpad
    f = Forces(label="Test_02", M_y=-400 * kip * ft)
    beam_example_flexure_ACI.set_longitudinal_rebar_bot(n1=2, d_b1=0.75 * inch)
    beam_example_flexure_ACI.set_longitudinal_rebar_top(n1=2, d_b1=1.41 * inch, n3=2, d_b3=1.27 * inch)
    node = Node(section=beam_example_flexure_ACI, forces=f)
    results = node.check_flexure()

    assert results.iloc[1]["Label"] == "B-12x24"
    assert results.iloc[1]["Comb."] == "Test_02"
    assert results.iloc[1]["Position"] == "Top"
    assert results.iloc[1]["As,min"] == pytest.approx(5.26, rel=1e-3)
    assert results.iloc[1]["As,req bot"] == pytest.approx(4.93, rel=1e-3)
    assert results.iloc[1]["As,req top"] == pytest.approx(33.17, rel=1e-3)
    assert results.iloc[1]["As"] == pytest.approx(36.49, rel=1e-3)
    assert results.iloc[1]["Mu"] == pytest.approx(-542.33, rel=1e-5)
    assert results.iloc[1]["ØMn"] == pytest.approx(588.73, rel=1e-3)


def test_check_flexure_ACI_318_19_3(beam_example_flexure_ACI: RectangularBeam) -> None:
    # Simple bending check
    # Ver calcpad: ACI 318-19 Beam Flexure 03.cpd
    f = Forces(label="Test_03", M_y=200 * kip * ft)
    beam_example_flexure_ACI.set_longitudinal_rebar_bot(n1=2, d_b1=1.41 * inch)
    beam_example_flexure_ACI.set_longitudinal_rebar_top(n1=2, d_b1=0.75 * inch)
    node = Node(section=beam_example_flexure_ACI, forces=f)
    results = node.check_flexure()

    assert results.iloc[1]["Label"] == "B-12x24"
    assert results.iloc[1]["Comb."] == "Test_03"
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["As,min"] == pytest.approx(5.53, rel=1e-3)
    assert results.iloc[1]["As,req bot"] == pytest.approx(14.51, rel=1e-3)
    assert results.iloc[1]["As,req top"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["As"] == pytest.approx(20.15, rel=1e-3)
    assert results.iloc[1]["ØMn"] == pytest.approx(364.37, rel=1e-3)


# def test_design_flexure_ACI_318_19_1(beam_example_flexure_ACI: RectangularBeam) -> None:
#     f = Forces(label='Test_01', M_y=400*kip*ft)
#     forces=[f]
#     results = beam_example_flexure_ACI._design_flexure(forces)

#     assert results.iloc[0]['Label'] == 'B-12x24'
#     assert results.iloc[0]['Load Combo']  == 'Test_01'
#     assert results.iloc[0]['Position'] == 'Bottom'
#     assert results.iloc[0]['As,min'].to(cm**2).magnitude == pytest.approx(5.26, rel=1e-2)
#     assert results.iloc[0]['As,req bot'].to(cm**2).magnitude == pytest.approx(33.17, rel=1e-3)
#     assert results.iloc[0]['As,req top'].to(cm**2).magnitude == pytest.approx(4.93, rel=1e-2)
#     assert results.iloc[0]['As'].to(inch**2).magnitude == pytest.approx(5.66, rel=1e-2)


def testing_determine_nominal_moment_ACI_318_19(
    beam_example_flexure_ACI: RectangularBeam,
) -> None:
    beam_example_flexure_ACI.set_longitudinal_rebar_bot(
        n1=2,
        d_b1=1.128 * inch,
        n2=1,
        d_b2=1.128 * inch,
        n3=2,
        d_b3=1 * inch,
        n4=1,
        d_b4=1 * inch,
    )
    beam_example_flexure_ACI.set_longitudinal_rebar_top(
        n1=2,
        d_b1=1.128 * inch,
        n2=1,
        d_b2=1.128 * inch,
        n3=2,
        d_b3=1 * inch,
        n4=1,
        d_b4=1 * inch,
    )

    f = Forces(label="Test_01", M_y=400 * kip * ft)
    node_1 = Node(section=beam_example_flexure_ACI, forces=f)
    node_1.check_flexure()
    assert beam_example_flexure_ACI._phi_M_n_bot.to("kN*m").magnitude == pytest.approx(587.0589108678, rel=1e-2)
    assert beam_example_flexure_ACI._phi_M_n_top.to("kN*m").magnitude == pytest.approx(587.0589108678, rel=1e-2)


def test_rectangular_section_plot_components(
    beam_example_flexure_ACI: RectangularBeam,
) -> None:
    """
    Test that the plot method adds the expected matplotlib components
    (main rectangle, stirrup rectangles, annotations).
    """
    beam_example_flexure_ACI.plot()
    ax = beam_example_flexure_ACI._ax  # Access the private _ax created by plot

    assert ax is not None, "Plot method should create and assign an Axes object."

    # Check for patches (rectangles, FancyBboxPatch)
    # Expecting: 1 main Rectangle, 2 FancyBboxPatches for stirrup
    # + potentially other internal patches depending on matplotlib's rendering
    # It's safer to check for specific types and properties if possible.

    # Filter for Rectangle and FancyBboxPatch
    rectangles = [p for p in ax.patches if isinstance(p, Rectangle)]
    fancy_bboxes = [p for p in ax.patches if isinstance(p, FancyBboxPatch)]

    assert len(rectangles) >= 1, "Expected at least one Rectangle patch (the main section)."
    assert len(fancy_bboxes) == 2, "Expected exactly two FancyBboxPatch objects for the stirrup."

    # Check for annotations (dimensions)
    # Annotations are stored in ax.texts for text, ax.lines or ax.artists for arrows
    # matplotlib.axes.Axes.annotate returns an Annotation object.
    # Text labels are found in ax.texts
    assert len(ax.texts) == 6, "Expected 6 text annotations for dimensions (width, height)."

    # You could add more specific checks, e.g.:
    # - Check the coordinates of the main rectangle:
    main_rect = next((p for p in rectangles if p.get_x() == 0 and p.get_y() == 0), None)
    assert main_rect is not None
    assert np.isclose(main_rect.get_width() / 2.54, beam_example_flexure_ACI.width.magnitude)
    assert np.isclose(main_rect.get_height() / 2.54, beam_example_flexure_ACI.height.magnitude)
    assert main_rect.get_edgecolor() == to_rgba(CUSTOM_COLORS["dark_gray"])

    # - Check stirrup patch properties (example, will need exact dimensions based on calculation)
    #   This gets more complex quickly, focusing on existence and basic properties is often enough.

    plt.close()  # Close the figure created by the plot method


def test_rectangular_section_plot_metric_units_text(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    """
    Test that the plot method displays dimensions in metric units (cm)
    when the concrete's unit system is not 'imperial'.
    """
    beam_example_EN_1992_2004_01.plot()  # Ensure the plot method is called
    ax = beam_example_EN_1992_2004_01._ax
    assert ax is not None, "The plot method did not assign an Axes object to _ax."

    # Extract text labels
    text_labels = [t.get_text() for t in ax.texts]

    # Check if the labels contain "cm"
    assert any("cm" in label for label in text_labels), "Expected metric unit 'cm' in dimension text."
    assert "{:.0f~P}".format(beam_example_EN_1992_2004_01.width.to("cm")) in text_labels
    assert "{:.0f~P}".format(beam_example_EN_1992_2004_01.height.to("cm")) in text_labels

    plt.close()


def test_flexural_check_fails_min_top_with_flag() -> None:
    """Test flexural check fails minimum for top reinforcement with _A_s_bool_top=True."""
    # Create beam with insufficient top rebar
    concrete = Concrete_ACI_318_19(name="H35", f_c=35 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="TEST",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=40 * cm,
        c_c=2.5 * cm,
    )
    f = Forces(label="ELU", M_y=-20 * kNm)
    node = Node(section=beam, forces=f)
    beam.set_longitudinal_rebar_top(n1=2, d_b1=12 * mm)

    # Check flexure
    node.check_flexure()

    checks = beam._data_min_max_flexure["Ok?"]

    # Should contain the article reference "9.6.1.3" for failing minimum
    assert "9.6.1.3" in checks[0]


def test_flexural_check_fails_min_bot_with_flag() -> None:
    """Test flexural check fails minimum for bottom reinforcement with _A_s_bool_bot=True."""
    # Create beam with insufficient bottom rebar
    concrete = Concrete_ACI_318_19(name="H35", f_c=35 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="TEST",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=40 * cm,
        c_c=2.5 * cm,
    )
    f = Forces(label="ELU", M_y=20 * kNm)
    node = Node(section=beam, forces=f)
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=12 * mm)

    # Check flexure
    node.check_flexure()

    # Check that the article reference appears in the checks
    checks = beam._data_min_max_flexure["Ok?"]
    assert "9.6.1.3" in checks[2]  # Position 2 is bottom rebar check


def test_flexural_check_fails_without_flag() -> None:
    """Test flexural check fails with ❌ when flags are False."""
    conc = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)

    beam = RectangularBeam(
        label="TEST",
        concrete=conc,
        steel_bar=steel,
        width=20 * cm,
        height=50 * cm,
        c_c=2.5 * cm,
    )

    # Set minimal top rebar
    beam.set_longitudinal_rebar_top(n1=1, d_b1=8 * mm)
    beam._A_s_bool_top = False  # Flag disabled

    forces = Forces(M_y=-50 * kNm, N_x=0 * kN, V_z=0 * kN)
    node = Node(beam, forces)

    node.check_flexure()

    # Should contain "❌" instead of article reference
    checks = beam._data_min_max_flexure["Ok?"]
    assert checks[0] == "❌"


def test_low_concrete_strength_negative_sqrt() -> None:
    """Test that very low concrete strength (f_c=1MPa) triggers sqrt_value < 0 condition."""
    # Very low concrete strength will cause sqrt_value to be negative
    conc = Concrete_ACI_318_19(name="C1", f_c=35 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)

    beam = RectangularBeam(
        label="TEST",
        concrete=conc,
        steel_bar=steel,
        width=20 * cm,
        height=60 * cm,
        c_c=2.5 * cm,
    )

    # Apply large moment that will cause negative sqrt_value
    forces = Forces(M_y=1000 * kNm)
    node = Node(beam, forces)

    # Design should not crash, but DCR should be > 1
    node.design_flexure()

    # The designed reinforcement should equal A_s_max due to negative sqrt_value
    # Check that beam has bottom rebar set
    assert beam._A_s_req_bot == beam._A_s_max_bot


def test_doubly_reinforced_ignores_max_limits() -> None:
    """Test that doubly reinforced sections bypass maximum reinforcement limits."""
    conc = Concrete_ACI_318_19(name="C35", f_c=35 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)

    beam = RectangularBeam(
        label="TEST",
        concrete=conc,
        steel_bar=steel,
        width=20 * cm,
        height=60 * cm,
        c_c=2.5 * cm,
    )

    # Large moment to trigger doubly reinforced condition
    forces = Forces(M_y=600 * kNm)
    node = Node(beam, forces)

    # Check flexure
    node.design()

    # Check for "✔️ D.R." in the checks
    checks = beam._data_min_max_flexure["Ok?"]
    # Either position 0 (top) or position 2 (bottom) should have "✔️ D.R."
    assert "✔️ D.R." in checks[0] or "✔️ D.R." in checks[2]


## BEAM RESULTS TESTS


@pytest.fixture
def sample_beam() -> RectangularBeam:
    """Create a sample beam for testing."""
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)

    beam = RectangularBeam(
        label="TEST",
        concrete=concrete,
        steel_bar=steel,
        width=20 * cm,
        height=50 * cm,
        c_c=2.5 * cm,
    )
    return beam


def test_beam_data_property(sample_beam: RectangularBeam) -> None:
    """Test the data property displays beam information."""
    with patch("mento.beam.display") as mock_display:
        # Call the property
        result = sample_beam.data

        # Should return None
        assert result is None

        # Check that display was called
        mock_display.assert_called_once()

        # Check that markdown content was stored
        assert hasattr(sample_beam, "_md_data")
        assert sample_beam.label in sample_beam._md_data
        assert "Concrete" in sample_beam._md_data


def test_flexure_results_property_not_checked(sample_beam: RectangularBeam) -> None:
    """Test flexure_results property when flexure not checked."""
    with patch("mento.beam.display"):
        with pytest.warns(UserWarning, match="Flexural design has not been performed"):
            result = sample_beam.flexure_results

        assert result is None
        assert sample_beam._md_flexure_results == "Flexural results are not available."


def test_flexure_results_property_after_check() -> None:
    """Test flexure_results property after checking."""
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="TEST",
        concrete=concrete,
        steel_bar=steel,
        width=20 * cm,
        height=50 * cm,
        c_c=2.5 * cm,
    )
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    beam.set_longitudinal_rebar_top(n1=2, d_b1=12 * mm)

    forces = Forces(M_y=50 * kNm)
    node = Node(beam, forces)
    node.check_flexure()

    with patch("mento.beam.display") as mock_display:
        result = beam.flexure_results

        # Should return None
        assert result is None

        # Display should have been called
        mock_display.assert_called_once()

        # Markdown content should contain flexure info
        assert hasattr(beam, "_md_flexure_results")
        assert "longitudinal rebar" in beam._md_flexure_results.lower()


def test_shear_results_property_not_checked(sample_beam: RectangularBeam) -> None:
    """Test shear_results property when shear not checked."""
    with patch("mento.beam.display"):
        with pytest.warns(UserWarning, match="Shear design has not been performed"):
            result = sample_beam.shear_results

        assert result is None
        assert sample_beam._md_shear_results == "Shear results are not available."


def test_shear_results_property_after_check() -> None:
    """Test shear_results property after checking."""
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="TEST",
        concrete=concrete,
        steel_bar=steel,
        width=20 * cm,
        height=50 * cm,
        c_c=2.5 * cm,
    )
    beam.set_transverse_rebar(n_stirrups=2, d_b=8 * mm, s_l=20 * cm)

    forces = Forces(M_y=0 * kNm, N_x=0 * kN, V_z=50 * kN)
    node = Node(beam, forces)
    node.check_shear()

    with patch("mento.beam.display") as mock_display:
        result = beam.shear_results

        assert result is None
        mock_display.assert_called_once()
        assert hasattr(beam, "_md_shear_results")
        assert "shear" in beam._md_shear_results.lower()


def test_results_property_combined() -> None:
    """Test the combined results property."""
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="TEST",
        concrete=concrete,
        steel_bar=steel,
        width=20 * cm,
        height=50 * cm,
        c_c=2.5 * cm,
    )
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    beam.set_transverse_rebar(n_stirrups=2, d_b=8 * mm, s_l=20 * cm)

    forces = Forces(M_y=50 * kNm, N_x=0 * kN, V_z=50 * kN)
    node = Node(beam, forces)
    node.check()

    with patch("mento.beam.display") as mock_display:
        result = beam.results

        assert result is None
        # Display should be called multiple times (data, flexure, shear)
        assert mock_display.call_count >= 1


def test_flexure_results_detailed_not_checked(sample_beam: RectangularBeam) -> None:
    """Test flexure_results_detailed when not checked."""
    with pytest.warns(UserWarning, match="Flexural check has not been performed"):
        result = sample_beam.flexure_results_detailed()
        assert result is None


def test_flexure_results_detailed_after_check(capsys) -> None:
    """Test flexure_results_detailed prints output."""
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="TEST",
        concrete=concrete,
        steel_bar=steel,
        width=20 * cm,
        height=50 * cm,
        c_c=2.5 * cm,
    )
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    beam.set_longitudinal_rebar_top(n1=2, d_b1=12 * mm)

    forces = Forces(M_y=50 * kNm)
    node = Node(beam, forces)
    node.check_flexure()

    beam.flexure_results_detailed()

    # Capture printed output
    captured = capsys.readouterr()
    assert "BEAM FLEXURE DETAILED RESULTS" in captured.out
    assert "Materials" in captured.out


def test_shear_results_detailed_not_checked(sample_beam) -> None:
    """Test shear_results_detailed when not checked."""
    with pytest.warns(UserWarning, match="Shear check has not been performed"):
        result = sample_beam.shear_results_detailed()
        assert result is None


def test_shear_results_detailed_after_check(capsys) -> None:
    """Test shear_results_detailed prints output."""
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="TEST",
        concrete=concrete,
        steel_bar=steel,
        width=20 * cm,
        height=50 * cm,
        c_c=2.5 * cm,
    )
    beam.set_transverse_rebar(n_stirrups=2, d_b=8 * mm, s_l=20 * cm)

    forces = Forces(M_y=0 * kNm, N_x=0 * kN, V_z=50 * kN)
    node = Node(beam, forces)
    node.check_shear()

    beam.shear_results_detailed()

    captured = capsys.readouterr()
    assert "BEAM SHEAR DETAILED RESULTS" in captured.out
    assert "Materials" in captured.out


def test_flexure_results_detailed_with_specific_force() -> None:
    """Test flexure_results_detailed with specific Forces object."""
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="TEST",
        concrete=concrete,
        steel_bar=steel,
        width=20 * cm,
        height=50 * cm,
        c_c=2.5 * cm,
    )
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)

    forces = Forces(label="F1", M_y=50 * kNm, N_x=0 * kN, V_z=0 * kN)
    node = Node(beam, forces)
    node.check_flexure()

    # Should work with the specific force
    beam.flexure_results_detailed(force=forces)


def test_flexure_results_detailed_invalid_force() -> None:
    """Test flexure_results_detailed with invalid Forces object."""
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="TEST",
        concrete=concrete,
        steel_bar=steel,
        width=20 * cm,
        height=50 * cm,
        c_c=2.5 * cm,
    )
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)

    forces1 = Forces(M_y=50 * kNm)
    node = Node(beam, forces1)
    node.check_flexure()

    # Try with a different force that wasn't checked
    forces2 = Forces(M_y=100 * kNm)

    with pytest.raises(ValueError, match="No results found for Forces object"):
        beam.flexure_results_detailed(force=forces2)


# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()
