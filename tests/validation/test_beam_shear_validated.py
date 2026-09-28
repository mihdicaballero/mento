"""Beam shear checks and designs that reproduce a case validated outside mento.

Every test here is marked ``published_example`` and says in a ``Source:`` paragraph of
its docstring where its numbers come from: the Calcpad beam-shear sheets (kept outside
the repository) and the EN 1992-1-1 shear calculators of eurocodeapplied.com.
``tests/architecture/test_published_examples.py`` enforces both.
"""

import pytest
from mento.node import Node
from mento.beam import RectangularBeam
from mento.units import kip, inch, mm, kN, cm
from mento.forces import Forces


@pytest.mark.published_example
def test_shear_check_EN_1992_2004_rebar_1(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    """Stirrups Ø6/25 under V_Ed = 100 kN.

    Source: eurocodeapplied.com, shear resistance of a reinforced concrete beam, EN 1992-1-1 §6.2.
    """
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


@pytest.mark.published_example
def test_shear_check_EN_1992_2004_rebar_2(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    """Stirrups under V_Ed = 350 kN.

    Source: eurocodeapplied.com, shear resistance of a reinforced concrete beam, EN 1992-1-1 §6.2.
    """
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


@pytest.mark.published_example
def test_shear_check_EN_1992_2004_rebar_3(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    """Stirrups under V_Ed = 500 kN.

    Source: eurocodeapplied.com, shear resistance of a reinforced concrete beam, EN 1992-1-1 §6.2.
    """
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


@pytest.mark.published_example
def test_shear_check_EN_1992_2004_no_rebar_1(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    """A beam with no stirrups at all under a shear V_Rd,c carries alone.

    Source: Calcpad "EN 1992-1-1_2004 Beam Shear 01 - Metric.cpd".
    """
    f = Forces(V_z=30 * kN)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    # The reference example has no stirrups at all, so say so: the section no
    # longer infers it from n_stirrups == 0 while still assuming the settings'
    # initial diameter for d.
    beam_example_EN_1992_2004_01.set_transverse_rebar(n_stirrups=0, d_b=0 * mm, s_l=0 * mm)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VEd,1"] == pytest.approx(30, rel=1e-3)
    assert results.iloc[1]["VEd,2"] == pytest.approx(30, rel=1e-3)
    assert beam_example_EN_1992_2004_01.shear_design.d_b.to("mm").magnitude == 0
    assert beam_example_EN_1992_2004_01._d_shear.to("cm").magnitude == pytest.approx(56.6, rel=1e-3)
    assert results.iloc[1]["VRd,c"] == pytest.approx(56.51, rel=1e-3)
    assert results.iloc[1]["VRd,s"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VRd"] == pytest.approx(56.51, rel=1e-3)
    assert results.iloc[1]["VRd,max"] == pytest.approx(56.51, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.531, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["VEd,1≤VRd,max"] is True
    assert results.iloc[1]["VEd,2≤VRd"] is True


@pytest.mark.published_example
def test_shear_check_EN_1992_2004_no_rebar_2(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    """No stirrups under V_Ed = 30 kN.

    Source: eurocodeapplied.com, shear resistance of a reinforced concrete beam, EN 1992-1-1 §6.2.
    """
    f = Forces(V_z=30 * kN)
    # The reference example has no stirrups at all, so say so: the section no
    # longer infers it from n_stirrups == 0 while still assuming the settings'
    # initial diameter for d.
    beam_example_EN_1992_2004_01.set_transverse_rebar(n_stirrups=0, d_b=0 * mm, s_l=0 * mm)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VEd,1"] == pytest.approx(30, rel=1e-3)
    assert results.iloc[1]["VEd,2"] == pytest.approx(30, rel=1e-3)
    assert beam_example_EN_1992_2004_01.shear_design.d_b.to("mm").magnitude == 0
    assert beam_example_EN_1992_2004_01._d_shear.to("cm").magnitude == pytest.approx(57, rel=1e-3)
    assert results.iloc[1]["VRd,c"] == pytest.approx(40.09, rel=1e-3)
    assert results.iloc[1]["VRd,s"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VRd"] == pytest.approx(40.09, rel=1e-3)
    assert results.iloc[1]["VRd,max"] == pytest.approx(40.09, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.748, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["VEd,1≤VRd,max"] is True
    assert results.iloc[1]["VEd,2≤VRd"] is True


@pytest.mark.published_example
def test_shear_check_EN_1992_2004_no_rebar_3(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    """No stirrups under V_Ed = 30 kN with N_Ed = 50 kN.

    Source: eurocodeapplied.com, shear resistance of a reinforced concrete beam, EN 1992-1-1 §6.2.
    """
    f = Forces(N_x=50 * kN, V_z=30 * kN)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    # The reference example has no stirrups at all, so say so: the section no
    # longer infers it from n_stirrups == 0 while still assuming the settings'
    # initial diameter for d.
    beam_example_EN_1992_2004_01.set_transverse_rebar(n_stirrups=0, d_b=0 * mm, s_l=0 * mm)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VEd,1"] == pytest.approx(30, rel=1e-3)
    assert results.iloc[1]["VEd,2"] == pytest.approx(30, rel=1e-3)
    assert beam_example_EN_1992_2004_01.shear_design.d_b.to("mm").magnitude == 0
    assert beam_example_EN_1992_2004_01._d_shear.to("cm").magnitude == pytest.approx(56.6, rel=1e-3)
    assert results.iloc[1]["VRd,c"] == pytest.approx(63.59, rel=1e-3)
    assert results.iloc[1]["VRd,s"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["VRd"] == pytest.approx(63.59, rel=1e-3)
    assert results.iloc[1]["VRd,max"] == pytest.approx(63.59, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.472, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["VEd,1≤VRd,max"] is True
    assert results.iloc[1]["VEd,2≤VRd"] is True


@pytest.mark.published_example
def test_shear_design_EN_1992_2004_1(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    """Stirrup design under V_Ed = 30 kN.

    Source: eurocodeapplied.com, shear resistance of a reinforced concrete beam, EN 1992-1-1 §6.2.
    """
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


@pytest.mark.published_example
def test_shear_check_ACI_318_19_1(beam_example_imperial: RectangularBeam) -> None:
    """Stirrups 1e#4/6 in under V_u = 37.727 kip, no axial load.

    Source: Calcpad "ACI 318-19 Beam Shear 01 - Imperial.cpd".
    """
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


@pytest.mark.published_example
def test_shear_check_ACI_318_19_2(beam_example_imperial: RectangularBeam) -> None:
    """Stirrups 1e#4/6 in under V_u = 37.727 kip with N_u = 20 kip.

    Source: Calcpad "ACI 318-19 Beam Shear 01 - Imperial.cpd".
    """
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


@pytest.mark.published_example
def test_shear_check_ACI_318_19_no_rebar_1(
    beam_example_imperial: RectangularBeam,
) -> None:
    """A beam with no stirrups under a shear that needs them (V_u > phi*V_c).

    Source: Calcpad "ACI 318-19 Beam Shear 01 - Imperial.cpd", the case that needs rebar.
    """
    f = Forces(V_z=8 * kip, N_x=0 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    # The reference example has no stirrups at all, so say so: the section no
    # longer infers it from n_stirrups == 0 while still assuming the settings'
    # initial diameter for d.
    beam_example_imperial.set_transverse_rebar(n_stirrups=0, d_b=0 * mm, s_l=0 * mm)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.check_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert results.iloc[1]["Av,min"] == pytest.approx(2.12, rel=1e-3)
    assert results.iloc[1]["Av,req"] == pytest.approx(2.12, rel=1e-3)
    assert results.iloc[1]["Av"] == pytest.approx(0, rel=1e-3)
    assert beam_example_imperial.shear_design.d_b.to("mm").magnitude == 0
    assert beam_example_imperial._d_shear.to("cm").magnitude == pytest.approx(36.04, rel=1e-3)
    assert results.iloc[1]["ØVc"] == pytest.approx(35.48, rel=1e-3)
    assert results.iloc[1]["ØVs"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["ØVn"] == pytest.approx(35.48, rel=1e-3)
    assert results.iloc[1]["ØVmax"] == pytest.approx(274.96, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(1.003, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[1]["Vu≤ØVmax"] is True
    assert results.iloc[1]["Vu≤ØVn"] is False


@pytest.mark.published_example
def test_shear_check_ACI_318_19_no_rebar_2(
    beam_example_imperial: RectangularBeam,
) -> None:
    """A beam with no stirrups under a shear that does not need them.

    Source: Calcpad "ACI 318-19 Beam Shear 01 - Imperial.cpd", the case that does not need rebar.
    """
    f = Forces(V_z=6 * kip, N_x=0 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    # The reference example has no stirrups at all, so say so: the section no
    # longer infers it from n_stirrups == 0 while still assuming the settings'
    # initial diameter for d.
    beam_example_imperial.set_transverse_rebar(n_stirrups=0, d_b=0 * mm, s_l=0 * mm)
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


@pytest.mark.published_example
def test_shear_design_ACI_318_19(beam_example_imperial: RectangularBeam) -> None:
    """Stirrup design for the beam that needs rebar.

    Source: Calcpad "ACI 318-19 Beam Shear 01 - Imperial.cpd", the case that needs rebar.
    """
    f = Forces(V_z=37.727 * kip, N_x=0 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    node = Node(section=beam_example_imperial, forces=f)
    results = node.design_shear()

    # Compare dictionaries with a tolerance for floating-point values, in m
    assert beam_example_imperial._d_shear.to("cm").magnitude == pytest.approx(35.08, rel=1e-3)
    assert results.iloc[1]["Av,min"] == pytest.approx(2.12, rel=1e-3)
    # Nominal required Vs = (Vu - phi*Vc)/phi = (167.82 - 58.29)/0.75 = 146.04 kN
    # (it held phi*Vs = 109.53 kN before the Table 9.7.6.2.2 fix; Av,req is unchanged).
    assert beam_example_imperial._V_s_req.to("kN").magnitude == pytest.approx(146.04, rel=1e-3)
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
    shear = beam_example_imperial.shear_design
    assert shear.n_stirrups == 1
    assert shear.d_b.to("mm").magnitude == pytest.approx(9.525, rel=1e-3)
    assert shear.s_l.to("cm").magnitude == pytest.approx(12.7, rel=1e-3)
