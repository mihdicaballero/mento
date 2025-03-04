import pytest
from mento.rebar import Rebar
from mento.beam import RectangularBeam
from mento.material import Concrete_ACI_318_19, SteelBar
from mento.units import psi, kip, mm, inch, ksi, cm, MPa, kN
from mento.forces import Forces
from mento.node import Node

@pytest.fixture()
def beam_example_imperial() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000*psi)  
    steelBar = SteelBar(name="ADN 420", f_y=60*ksi) 
    custom_settings = {'clear_cover': 1.5*inch, 'stirrup_diameter_ini':0.5*inch,
                       'longitudinal_diameter_ini': 1*inch} 
    section = RectangularBeam(
        label="V-10x16",
        concrete=concrete,
        steel_bar=steelBar,
        width=10*inch,  
        height=16*inch,
        settings=custom_settings  
    )
    return section

@pytest.fixture()
def beam_example_metric() -> RectangularBeam:
    concrete= Concrete_ACI_318_19(name="H30",f_c=30*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa)
    custom_settings = {'clear_cover': 30*mm, 'stirrup_diameter_ini':8*mm}
    section = RectangularBeam(
        label="V 20x50",
        concrete=concrete,
        steel_bar=steelBar,
        width=20*cm,  
        height=50*cm,
        settings=custom_settings  
    )
    return section

# Tests
def test_beam_longitudinal_rebar_ACI_318_19_metric(beam_example_metric: RectangularBeam) -> None:
    A_s_req = 5 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design['n_1'] == 2
    assert best_design['d_b1'] == 20*mm
    assert best_design['n_2'] == 0
    assert best_design['d_b2'] is None
    assert best_design['n_3'] == 0
    assert best_design['d_b3'] is None
    assert best_design['n_4'] == 0
    assert best_design['d_b3'] is None
    assert best_design['total_as'].magnitude == pytest.approx(6.28,rel=1e-3)
    assert best_design['total_bars'] == 2
    assert best_design['clear_spacing'].magnitude == pytest.approx(84,rel=1e-3)

# def test_beam_longitudinal_rebar_ACI_318_19_imperial(beam_example_imperial: RectangularBeam) -> None:
#     # This test will ensure that the longitudinal rebar design works correctly for beams defined in imperial units.
#     A_s_req = 1.5 * inch**2
#     beam_rebar = Rebar(beam_example_imperial)
#     beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
#     best_design = beam_rebar.longitudinal_rebar_design

#     assert best_design is not None
#     assert best_design['n_1'] == 2
#     assert best_design['d_b1'] == 0.625 * inch
#     assert best_design['n_2'] == 1
#     assert best_design['d_b2'] == 0.5 * inch
#     assert best_design['n_3'] == 0
#     assert best_design['d_b3'] is None
#     assert best_design['n_4'] == 0
#     assert best_design['d_b3'] is None
#     assert best_design['total_as'].magnitude == pytest.approx(1.53, rel=1e-3)
#     assert best_design['total_bars'] == 3
#     assert best_design['clear_spacing'].magnitude == pytest.approx(1.5, rel=1e-3)

def test_beam_longitudinal_rebar_min_diameter(beam_example_metric: RectangularBeam) -> None:
    # This test ensures that the minimum rebar diameter constraint is respected.
    A_s_req = 2 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design['d_b1'] >= 10 * mm
    assert best_design['d_b2'] >= 10 * mm if best_design['d_b2'] is not None else True

def test_beam_longitudinal_rebar_max_bars_per_layer(beam_example_metric: RectangularBeam) -> None:
    # This test ensures that the maximum number of bars per layer constraint is respected.
    A_s_req = 10 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design['n_1'] + best_design['n_2'] <= beam_rebar.max_bars_per_layer
    assert best_design['n_3'] + best_design['n_4'] <= beam_rebar.max_bars_per_layer

def test_beam_longitudinal_rebar_CIRSOC_201_25(beam_example_metric: RectangularBeam) -> None:
    # This test ensures that the longitudinal rebar design works correctly for different design codes.
    beam_example_metric.concrete.design_code = "CIRSOC 201-25"
    A_s_req = 5 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design['n_1'] == 2
    assert best_design['d_b1'] == 20 * mm
    assert best_design['n_2'] == 0
    assert best_design['d_b2'] is None
    assert best_design['n_3'] == 0
    assert best_design['d_b3'] is None
    assert best_design['n_4'] == 0
    assert best_design['d_b3'] is None
    assert best_design['total_as'].magnitude == pytest.approx(6.28, rel=1e-3)
    assert best_design['total_bars'] == 2
    assert best_design['clear_spacing'].magnitude == pytest.approx(84, rel=1e-3)

def test_beam_longitudinal_rebar_small_area(beam_example_metric: RectangularBeam) -> None:
    # This test ensures that the rebar design handles very small required areas correctly.
    A_s_req = 0.1 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design['total_as'] >= A_s_req

def test_beam_longitudinal_rebar_large_area(beam_example_metric: RectangularBeam) -> None:
    # This test ensures that the rebar design handles very large required areas correctly.
    A_s_req = 50 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design['n_1'] >= 2
    assert best_design['d_b1'] >= 8 * mm

def test_beam_longitudinal_rebar_zero_area_metric_ACI(beam_example_metric: RectangularBeam) -> None:
    A_s_req = 0 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design['total_as'].magnitude == pytest.approx(1.5707,rel=1e-3) # 2Ø10 minimum

# def test_beam_longitudinal_rebar_zero_area_imperial(beam_example_imperial: RectangularBeam) -> None:
#     A_s_req = 0 * inch**2
#     beam_rebar = Rebar(beam_example_imperial)
#     beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
#     best_design = beam_rebar.longitudinal_rebar_design

#     assert best_design is not None
#     assert best_design['total_as'].magnitude == pytest.approx(1.5707,rel=1e-3) # 2#3 minimum

def test_beam_transverse_rebar_ACI_318_19_imperial(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=37.727*kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625*inch)
    N1=Node(beam_example_imperial, forces=f)
    results = N1.design_shear()
    assert results is not None
    assert beam_example_imperial._stirrup_n == 1
    assert beam_example_imperial._stirrup_d_b == 3/8*inch
    assert beam_example_imperial._stirrup_s_l ==  5*inch 
    assert beam_example_imperial._A_v.to('cm**2/m').magnitude == pytest.approx(11.2214,rel=1e-3)

def test_beam_transverse_rebar_ACI_318_19_metric(beam_example_metric: RectangularBeam) -> None:
    # This test ensures that the transverse rebar design works correctly for beams defined in metric units.
    f = Forces(V_z=100 * kN)
    beam_example_metric.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)
    N1=Node(beam_example_metric, forces=f)
    results = N1.design_shear()
    
    assert results is not None
    assert beam_example_metric._stirrup_n == 1
    assert beam_example_metric._stirrup_d_b == 10 * mm
    assert beam_example_metric._stirrup_s_l == 22 * cm
    assert beam_example_metric._A_v.to('cm**2/m').magnitude == pytest.approx(7.14, rel=1e-3)

def test_beam_transverse_rebar_max_spacing_imperial(beam_example_imperial: RectangularBeam) -> None:
    # This test ensures that the maximum spacing constraints are respected.
    f = Forces(V_z=50 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    N1=Node(beam_example_imperial, forces=f)
    results = N1.design_shear()
    
    assert results is not None
    assert beam_example_imperial._stirrup_s_l <= 24 * inch
    assert beam_example_imperial._stirrup_s_w <= 24 * inch

def test_beam_transverse_rebar_min_diameter(beam_example_metric: RectangularBeam) -> None:
    # This test ensures that the minimum rebar diameter constraint is respected for transverse reinforcement.
    f = Forces(V_z=50 * kN)
    beam_example_metric.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)
    N1=Node(beam_example_metric, forces=f)
    results = N1.design_shear()
    
    assert results is not None
    assert beam_example_metric._stirrup_d_b >= 8 * mm

def test_beam_transverse_rebar_CIRSOC_201_25(beam_example_metric: RectangularBeam) -> None:
    # This test ensures that the transverse rebar design works correctly for different design codes.
    beam_example_metric.concrete.design_code = "CIRSOC 201-25"
    f = Forces(V_z=100 * kN)
    beam_example_metric.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)
    N1=Node(beam_example_metric, forces=f)
    results = N1.design_shear()
    
    assert results is not None
    assert beam_example_metric._stirrup_n == 1
    assert beam_example_metric._stirrup_d_b == 6 * mm
    assert beam_example_metric._stirrup_s_l == 22 * cm
    assert beam_example_metric._A_v.to('cm**2/m').magnitude == pytest.approx(2.57, rel=1e-3)

# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()