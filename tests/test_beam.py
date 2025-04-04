import pytest
import numpy as np

from mento.node import Node
from mento.beam import RectangularBeam
from mento.material import Concrete_ACI_318_19, SteelBar, Concrete_EN_1992_2004
from mento.units import psi, kip, inch, ksi, mm, kN, cm, MPa, ft, kNm
from mento.forces import Forces
from mento.codes.ACI_318_19_beam import (
    _determine_nominal_moment_ACI_318_19,
    _determine_nominal_moment_simple_reinf_ACI_318_19,
    _determine_nominal_moment_double_reinf_ACI_318_19,
)


@pytest.fixture()
def beam_example_imperial() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steelBar = SteelBar(name="ADN 420", f_y=60 * ksi)
    custom_settings = {"clear_cover": 1.5 * inch}
    section = RectangularBeam(
        label="V101",
        concrete=concrete,
        steel_bar=steelBar,
        width=10 * inch,
        height=16 * inch,
        settings=custom_settings,
    )
    return section


@pytest.fixture()
def beam_example_EN_1992_2004() -> RectangularBeam:
    concrete = Concrete_EN_1992_2004(name="C25", f_ck=25 * MPa)
    steelBar = SteelBar(name="B500S", f_y=500 * MPa)
    custom_settings = {"clear_cover": 2.6 * cm}
    section = RectangularBeam(
        label="V-20x60",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=60 * cm,
        settings=custom_settings,
    )
    return section


def test_shear_check_EN_1992_2004_rebar_1(
    beam_example_EN_1992_2004: RectangularBeam,
) -> None:
    f = Forces(V_z=100 * kN)
    beam_example_EN_1992_2004.set_transverse_rebar(
        n_stirrups=1, d_b=6 * mm, s_l=25 * cm
    )
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(1.825, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(2.262, rel=1e-3)
    assert results.iloc[0]["VEd,1"].magnitude == pytest.approx(100, rel=1e-3)
    assert results.iloc[0]["VEd,2"].magnitude == pytest.approx(100, rel=1e-3)
    assert results.iloc[0]["VRd,c"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["VRd,s"].magnitude == pytest.approx(123.924, rel=1e-3)
    assert results.iloc[0]["VRd"].magnitude == pytest.approx(123.924, rel=1e-3)
    assert results.iloc[0]["VRd,max"].magnitude == pytest.approx(312.811, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(0.8069, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["VEd,1<VRd,max"] is np.True_
    assert results.iloc[0]["VEd,2<VRd"] is np.True_


def test_shear_check_EN_1992_2004_rebar_2(
    beam_example_EN_1992_2004: RectangularBeam,
) -> None:
    f = Forces(V_z=350 * kN)
    beam_example_EN_1992_2004.set_transverse_rebar(
        n_stirrups=1, d_b=6 * mm, s_l=25 * cm
    )
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(7.533, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(2.262, rel=1e-3)
    assert results.iloc[0]["VEd,1"].magnitude == pytest.approx(350, rel=1e-3)
    assert results.iloc[0]["VEd,2"].magnitude == pytest.approx(350, rel=1e-3)
    assert results.iloc[0]["VRd,c"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["VRd,s"].magnitude == pytest.approx(105.099, rel=1e-3)
    assert results.iloc[0]["VRd"].magnitude == pytest.approx(105.099, rel=1e-3)
    assert results.iloc[0]["VRd,max"].magnitude == pytest.approx(350, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(3.33, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["VEd,1<VRd,max"] is np.True_
    assert results.iloc[0]["VEd,2<VRd"] is np.False_


def test_shear_check_EN_1992_2004_rebar_3(
    beam_example_EN_1992_2004: RectangularBeam,
) -> None:
    f = Forces(V_z=500 * kN)
    beam_example_EN_1992_2004.set_transverse_rebar(
        n_stirrups=1, d_b=6 * mm, s_l=25 * cm
    )
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(22.817, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(2.26, rel=1e-3)
    assert results.iloc[0]["VEd,1"].magnitude == pytest.approx(500, rel=1e-3)
    assert results.iloc[0]["VEd,2"].magnitude == pytest.approx(500, rel=1e-3)
    assert results.iloc[0]["VRd,c"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["VRd,s"].magnitude == pytest.approx(49.566, rel=1e-3)
    assert results.iloc[0]["VRd"].magnitude == pytest.approx(49.566, rel=1e-3)
    assert results.iloc[0]["VRd,max"].magnitude == pytest.approx(453.6, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(10.088, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["VEd,1<VRd,max"] is np.False_
    assert results.iloc[0]["VEd,2<VRd"] is np.False_


def test_shear_check_EN_1992_2004_no_rebar_1(
    beam_example_EN_1992_2004: RectangularBeam,
) -> None:
    f = Forces(V_z=30 * kN)
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["VEd,1"].magnitude == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]["VEd,2"].magnitude == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]["VRd,c"].magnitude == pytest.approx(56.126, rel=1e-3)
    assert results.iloc[0]["VRd,s"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["VRd"].magnitude == pytest.approx(56.126, rel=1e-3)
    assert results.iloc[0]["VRd,max"].magnitude == pytest.approx(56.126, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(0.5345, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["VEd,1<VRd,max"] is np.True_
    assert results.iloc[0]["VEd,2<VRd"] is np.True_


def test_shear_check_EN_1992_2004_no_rebar_2(
    beam_example_EN_1992_2004: RectangularBeam,
) -> None:
    f = Forces(V_z=30 * kN)
    node = Node(section=beam_example_EN_1992_2004, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["VEd,1"].magnitude == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]["VEd,2"].magnitude == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]["VRd,c"].magnitude == pytest.approx(39.681, rel=1e-3)
    assert results.iloc[0]["VRd,s"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["VRd"].magnitude == pytest.approx(39.681, rel=1e-3)
    assert results.iloc[0]["VRd,max"].magnitude == pytest.approx(39.681, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(0.756, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["VEd,1<VRd,max"] is np.True_
    assert results.iloc[0]["VEd,2<VRd"] is np.True_


def test_shear_check_EN_1992_2004_no_rebar_3(
    beam_example_EN_1992_2004: RectangularBeam,
) -> None:
    f = Forces(N_x=50 * kN, V_z=30 * kN)
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["VEd,1"].magnitude == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]["VEd,2"].magnitude == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]["VRd,c"].magnitude == pytest.approx(63.101, rel=1e-3)
    assert results.iloc[0]["VRd,s"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["VRd"].magnitude == pytest.approx(63.101, rel=1e-3)
    assert results.iloc[0]["VRd,max"].magnitude == pytest.approx(63.101, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(0.4754, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["VEd,1<VRd,max"] is np.True_
    assert results.iloc[0]["VEd,2<VRd"] is np.True_


def test_shear_design_EN_1992_2004_1(
    beam_example_EN_1992_2004: RectangularBeam,
) -> None:
    f = Forces(N_x=0 * kN, V_z=30 * kN)
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004, forces=f)
    results = node.design_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(1.616, rel=1e-3)
    assert results.iloc[0]["VEd,1"].magnitude == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]["VEd,2"].magnitude == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]["VRd,c"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["VRd,s"].magnitude == pytest.approx(88.52, rel=1e-3)
    assert results.iloc[0]["VRd"].magnitude == pytest.approx(88.52, rel=1e-3)
    assert results.iloc[0]["VRd,max"].magnitude == pytest.approx(312.81, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(0.339, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["VEd,1<VRd,max"] is np.True_
    assert results.iloc[0]["VEd,2<VRd"] is np.True_


def test_shear_check_ACI_318_19_1(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=37.727 * kip, N_x=0 * kip)
    beam_example_imperial.set_transverse_rebar(
        n_stirrups=1, d_b=0.5 * inch, s_l=6 * inch
    )
    node = Node(section=beam_example_imperial, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(10.0623, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(16.624, rel=1e-3)
    assert results.iloc[0]["ØVc"].magnitude == pytest.approx(58.288, rel=1e-3)
    assert results.iloc[0]["ØVs"].magnitude == pytest.approx(180.956, rel=1e-3)
    assert results.iloc[0]["ØVn"].magnitude == pytest.approx(239.247, rel=1e-3)
    assert results.iloc[0]["ØVmax"].magnitude == pytest.approx(291.44, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(0.70144, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["Vu<ØVmax"] is np.True_
    assert results.iloc[0]["Vu<ØVn"] is np.True_


def test_shear_check_ACI_318_19_2(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=37.727 * kip, N_x=20 * kip)
    beam_example_imperial.set_transverse_rebar(
        n_stirrups=1, d_b=0.5 * inch, s_l=6 * inch
    )
    node = Node(section=beam_example_imperial, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(9.1803, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(16.624, rel=1e-3)
    assert results.iloc[0]["ØVc"].magnitude == pytest.approx(67.888, rel=1e-3)
    assert results.iloc[0]["ØVs"].magnitude == pytest.approx(180.959, rel=1e-3)
    assert results.iloc[0]["ØVn"].magnitude == pytest.approx(248.847, rel=1e-3)
    assert results.iloc[0]["ØVmax"].magnitude == pytest.approx(301.041, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(0.6743, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["Vu<ØVmax"] is np.True_
    assert results.iloc[0]["Vu<ØVn"] is np.True_


def test_shear_check_ACI_318_19_no_rebar_1(
    beam_example_imperial: RectangularBeam,
) -> None:
    f = Forces(V_z=8 * kip, N_x=0 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["ØVc"].magnitude == pytest.approx(35.1253, rel=1e-3)
    assert results.iloc[0]["ØVs"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["ØVn"].magnitude == pytest.approx(35.125, rel=1e-3)
    assert results.iloc[0]["ØVmax"].magnitude == pytest.approx(268.278, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(1.013, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["Vu<ØVmax"] is np.True_
    assert results.iloc[0]["Vu<ØVn"] is np.False_


def test_shear_check_ACI_318_19_no_rebar_2(
    beam_example_imperial: RectangularBeam,
) -> None:
    f = Forces(V_z=6 * kip, N_x=0 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["ØVc"].magnitude == pytest.approx(58.288, rel=1e-3)
    assert results.iloc[0]["ØVs"].magnitude == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]["ØVn"].magnitude == pytest.approx(58.288, rel=1e-3)
    assert results.iloc[0]["ØVmax"].magnitude == pytest.approx(291.441, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(0.4579, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["Vu<ØVmax"] is np.True_
    assert results.iloc[0]["Vu<ØVn"] is np.True_


def test_shear_design_ACI_318_19(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=37.727 * kip, N_x=0 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.design_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[0]["Av,min"].magnitude == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]["Av,req"].magnitude == pytest.approx(10.06, rel=1e-3)
    assert results.iloc[0]["Av"].magnitude == pytest.approx(11.22, rel=1e-3)
    assert results.iloc[0]["ØVc"].magnitude == pytest.approx(58.29, rel=1e-3)
    assert results.iloc[0]["ØVs"].magnitude == pytest.approx(122.15, rel=1e-3)
    assert results.iloc[0]["ØVn"].magnitude == pytest.approx(180.44, rel=1e-3)
    assert results.iloc[0]["ØVmax"].magnitude == pytest.approx(291.44, rel=1e-3)
    assert results.iloc[0]["DCR"] == pytest.approx(0.93, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]["Vu<ØVmax"] is np.True_
    assert results.iloc[0]["Vu<ØVn"] is np.True_


# # ------- Flexural test --------------


def test_flexure_check_EN_1992_2004_1(
    beam_example_EN_1992_2004: RectangularBeam,
) -> None:
    # Example from Calcpad Validation
    f = Forces(M_y=150 * kNm)
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam_example_EN_1992_2004, forces=f)
    results = node.check_flexure()

    assert results.iloc[0]["Position"] == "Bottom"
    assert results.iloc[0]["As,min"].to(cm**2).magnitude == pytest.approx(
        1.49, rel=1e-2
    )
    assert results.iloc[0]["As,req bot"].to(cm**2).magnitude == pytest.approx(
        6.6, rel=1e-3
    )
    assert results.iloc[0]["As,req top"].to(cm**2).magnitude == pytest.approx(
        0, rel=1e-3
    )
    assert results.iloc[0]["As"].to(inch**2).magnitude == pytest.approx(
        8.04, rel=1e-2
    )
    assert results.iloc[0]["MRd"].to(kip * ft).magnitude == pytest.approx(
        182.88, rel=1e-5
    )
    assert results.iloc[0]["DCR"] == pytest.approx(0.82, rel=1e-5)


def test_flexure_check_EN_1992_2004_2(
    beam_example_EN_1992_2004: RectangularBeam,
) -> None:
    # Example from Calcpad Validation
    f = Forces(M_y=450 * kNm)
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=25 * mm)
    node = Node(section=beam_example_EN_1992_2004, forces=f)
    results = node.check_flexure()

    assert results.iloc[0]["Position"] == "Bottom"
    assert results.iloc[0]["As,min"].to(cm**2).magnitude == pytest.approx(
        1.48, rel=1e-2
    )
    assert results.iloc[0]["As,req bot"].to(cm**2).magnitude == pytest.approx(
        24.13, rel=1e-3
    )
    assert results.iloc[0]["As,req top"].to(cm**2).magnitude == pytest.approx(
        3.19, rel=1e-3
    )
    assert results.iloc[0]["As"].to(inch**2).magnitude == pytest.approx(
        19.63, rel=1e-2
    )
    assert results.iloc[0]["MRd"].to(kip * ft).magnitude == pytest.approx(
        466.16, rel=1e-5
    )
    assert results.iloc[0]["DCR"] == pytest.approx(0.965, rel=1e-5)


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
    )
    A_s = 6 * inch**2
    d = 24 * inch
    d_prime = 2.5 * inch
    A_s_prime = 1.8 * inch**2
    result = _determine_nominal_moment_double_reinf_ACI_318_19(
        beam, A_s, d, d_prime, A_s_prime
    )
    assert result.to(kip * ft).magnitude == pytest.approx(639.12, rel=1e-2)


# BEAM FROM CALPCAD "ACI 318-19 Beam Flexure 02-TEST_1"
@pytest.fixture()
def beam_example_flexure_ACI() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="fc 4000", f_c=4000 * psi)
    steelBar = SteelBar(name="fy 60000", f_y=60 * ksi)
    custom_settings = {"clear_cover": 1.5 * inch}
    section = RectangularBeam(
        label="B-12x24",
        concrete=concrete,
        steel_bar=steelBar,
        width=12 * inch,
        height=24 * inch,
        settings=custom_settings,
    )
    return section


def test_check_flexure_ACI_318_19_1(beam_example_flexure_ACI: RectangularBeam) -> None:
    # Testing the check of the reinforced beam with a large moment that requires compression reinforcement, the moment being positive
    # # See calcpad: ACI 318-19 Beam Flexure 02-TEST_1.cpd
    f = Forces(label="Test_01", M_y=400 * kip * ft)
    beam_example_flexure_ACI.set_longitudinal_rebar_bot(
        n1=2, d_b1=1.41 * inch, n3=2, d_b3=1.27 * inch
    )
    beam_example_flexure_ACI.set_longitudinal_rebar_top(n1=2, d_b1=0.75 * inch)
    node = Node(section=beam_example_flexure_ACI, forces=f)
    results = node.check_flexure()

    assert results.iloc[0]["Section Label"] == "B-12x24"
    assert results.iloc[0]["Load Combo"] == "Test_01"
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


def test_check_flexure_ACI_318_19_2(beam_example_flexure_ACI: RectangularBeam) -> None:
    # Testeo el checkeo de la viga armada con un momento grande que requiere armadura de compresión, siendo el momento NEGATIVO.
    # Ver calcpad: ACI 318-19 Beam Flexure 02-TEST_1.cpd
    # Invierto el armado respecto al test anterior para poder usar el mismo calcpad
    f = Forces(label="Test_02", M_y=-400 * kip * ft)
    beam_example_flexure_ACI.set_longitudinal_rebar_bot(n1=2, d_b1=0.75 * inch)
    beam_example_flexure_ACI.set_longitudinal_rebar_top(
        n1=2, d_b1=1.41 * inch, n3=2, d_b3=1.27 * inch
    )
    node = Node(section=beam_example_flexure_ACI, forces=f)
    results = node.check_flexure()

    assert results.iloc[0]["Section Label"] == "B-12x24"
    assert results.iloc[0]["Load Combo"] == "Test_02"
    assert results.iloc[0]["Position"] == "Top"
    assert results.iloc[0]["As,min"].to(cm**2).magnitude == pytest.approx(
        5.26, rel=1e-2
    )
    assert results.iloc[0]["As,req bot"].to(cm**2).magnitude == pytest.approx(
        4.9250577, rel=1e-3
    )
    assert results.iloc[0]["As,req top"].to(cm**2).magnitude == pytest.approx(
        33.1712108, rel=1e-3
    )
    assert results.iloc[0]["As"].to(inch**2).magnitude == pytest.approx(
        5.66, rel=1e-2
    )
    assert results.iloc[0]["Mu"].to(kip * ft).magnitude == pytest.approx(-400, rel=1e-5)


# def test_design_flexure_ACI_318_19_1(beam_example_flexure_ACI: RectangularBeam) -> None:
#     f = Forces(label='Test_01', M_y=400*kip*ft)
#     forces=[f]
#     results = beam_example_flexure_ACI._design_flexure(forces)

#     assert results.iloc[0]['Section Label'] == 'B-12x24'
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
    assert beam_example_flexure_ACI._phi_M_n_bot.to(kNm).magnitude == pytest.approx(
        587.0589108678, rel=1e-2
    )
    assert beam_example_flexure_ACI._phi_M_n_top.to(kNm).magnitude == pytest.approx(
        587.0589108678, rel=1e-2
    )


# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()
