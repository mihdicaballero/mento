import pytest
from mento.rebar import Rebar
from mento.beam import RectangularBeam
from mento.material import Concrete_ACI_318_19, SteelBar
from mento.units import psi, kip, mm, inch, ksi, cm, MPa
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
def test_beam_longitudinal_rebar_ACI_318_19(beam_example_metric: RectangularBeam) -> None:
    A_s_req = 5 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design['n_1'] == 2
    assert best_design['d_b1'] == 16*mm
    assert best_design['n_2'] == 1
    assert best_design['d_b2'] == 12*mm
    assert best_design['n_3'] == 0
    assert best_design['d_b3'] is None
    assert best_design['n_4'] == 0
    assert best_design['d_b3'] is None
    assert best_design['total_as'].magnitude == pytest.approx(5.1522,rel=1e-3)
    assert best_design['total_bars'] == 3
    assert best_design['clear_spacing'].magnitude == pytest.approx(40,rel=1e-3)

def test_beam_longitudinal_rebar_ACI_318_19_fail(beam_example_imperial: RectangularBeam) -> None:
    as_nec = 150 * cm**2
    # Expecting a ValueError to be raised
    with pytest.raises(ValueError, match="Cannot fit the required reinforcement within the beam "
                       "width considering clear cover and spacing."):
        Rebar(beam_example_imperial).longitudinal_rebar_ACI_318_19(A_s_req=as_nec)

def test_beam_transverse_rebar_ACI_318_19(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=37.727*kip)
    Node(beam_example_imperial, forces=f)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625*inch)
    results = beam_example_imperial.design_shear()
    assert results is not None
    assert beam_example_imperial._stirrup_n == 1
    assert beam_example_imperial._stirrup_d_b == 3/8*inch
    assert beam_example_imperial._stirrup_s_l ==  5*inch 
    assert beam_example_imperial._A_v.to('cm**2/m').magnitude == pytest.approx(11.2214,rel=1e-3)

# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()