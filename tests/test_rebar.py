import pytest
from mento import material
from mento.concrete.beam import Beam
from mento.rebar import Rebar
from mento.units import psi, kip, mm, inch, ksi, cm, MPa

@pytest.fixture()
def beam_example_imperial() -> Beam:
    concrete = material.create_concrete(name="C4", f_c=4000*psi, design_code="ACI 318-19")  
    steelBar = material.SteelBar(name="ADN 420", f_y=60*ksi)  
    section = Beam(
        name="V-10x16",
        concrete=concrete,
        steel_bar=steelBar,
        width=10*inch,  
        height=16*inch,  
    )
    section.cc = 1.5*inch
    section.stirrup_d_b = 0.5*inch
    return section

@pytest.fixture()
def beam_example_metric() -> Beam:
    concrete=material.create_concrete(name="H30",f_c=30*MPa, design_code="ACI 318-19") 
    steelBar=material.SteelBar(name="ADN 420", f_y=420*MPa) 
    section = Beam(
        name="V 20x50",
        concrete=concrete,
        steel_bar=steelBar,
        width=20*cm,  
        height=50*cm,  
    )
    section.cc = 30*mm
    section.stirrup_d_b = 8*mm
    return section

# Tests
def test_beam_longitudinal_rebar_ACI_318_19(beam_example_metric: Beam) -> None:
    as_nec = 5 * cm**2
    best_combination = Rebar(beam_example_metric).beam_longitudinal_rebar_ACI_318_19(A_s_req=as_nec)

    assert best_combination is not None
    assert best_combination['layer_1'] == 1
    assert best_combination['num_bars_1'] == 3
    assert best_combination['diameter_1'] ==  16*mm 
    assert best_combination['total_as'].value == pytest.approx(0.0006033,rel=1e-3)
    assert best_combination['available_spacing_1'] == 38*mm 

def test_beam_longitudinal_rebar_ACI_318_19_fail(beam_example_imperial: Beam) -> None:
    as_nec = 150 * cm**2
    # Expecting a ValueError to be raised
    with pytest.raises(ValueError, match="Cannot fit the required reinforcement within the beam "
                       "width considering clear cover and spacing."):
        Rebar(beam_example_imperial).beam_longitudinal_rebar_ACI_318_19(A_s_req=as_nec)

def test_beam_transverse_rebar_ACI_318_19(beam_example_imperial: Beam) -> None:
    V_u=37.727*kip 
    N_u = 0*kip 
    A_s=0.847*inch**2 
    beam_example_imperial.design_shear(V_u, N_u, A_s) 
    best_combination = beam_example_imperial.transverse_rebar
    assert best_combination is not None
    assert best_combination['n_stirrups'] == 1
    assert best_combination['d_b'] == 12*mm 
    assert best_combination['spacing'] ==  7*inch 
    assert best_combination['A_v'].value == pytest.approx(0.001272,rel=1e-3)

# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()