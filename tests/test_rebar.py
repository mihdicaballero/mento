import pytest
from mento.rebar import Rebar
from mento.concrete.beam import RectangularConcreteBeam
from mento.material import Concrete_ACI_318_19, SteelBar
from mento.units import psi, kip, mm, inch, ksi, cm, MPa
from mento.forces import Forces

@pytest.fixture()
def beam_example_imperial() -> RectangularConcreteBeam:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000*psi)  
    steelBar = SteelBar(name="ADN 420", f_y=60*ksi)  
    section = RectangularConcreteBeam(
        label="V-10x16",
        concrete=concrete,
        steel_bar=steelBar,
        width=10*inch,  
        height=16*inch,  
    )
    section.cc = 1.5*inch
    section.stirrup_d_b = 0.5*inch
    return section

@pytest.fixture()
def beam_example_metric() -> RectangularConcreteBeam:
    concrete= Concrete_ACI_318_19(name="H30",f_c=30*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa) 
    section = RectangularConcreteBeam(
        label="V 20x50",
        concrete=concrete,
        steel_bar=steelBar,
        width=20*cm,  
        height=50*cm,  
    )
    section.cc = 30*mm
    section.stirrup_d_b = 8*mm
    return section

# Tests
def test_beam_longitudinal_rebar_ACI_318_19(beam_example_metric: RectangularConcreteBeam) -> None:
    A_s_req = 5 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    long_rebar_df = beam_rebar.beam_longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_combination = long_rebar_df.iloc[0]

    assert best_combination is not None
    assert best_combination['layer_1'] == 1
    assert best_combination['num_bars_1'] == 3
    assert best_combination['diameter_1'] ==  16*mm 
    assert best_combination['total_as'].magnitude == pytest.approx(6.033,rel=1e-3)
    assert best_combination['available_spacing_1'].magnitude == pytest.approx(38,rel=1e-3)

def test_beam_longitudinal_rebar_ACI_318_19_fail(beam_example_imperial: RectangularConcreteBeam) -> None:
    as_nec = 150 * cm**2
    # Expecting a ValueError to be raised
    with pytest.raises(ValueError, match="Cannot fit the required reinforcement within the beam "
                       "width considering clear cover and spacing."):
        Rebar(beam_example_imperial).beam_longitudinal_rebar_ACI_318_19(A_s_req=as_nec)

def test_beam_transverse_rebar_ACI_318_19(beam_example_imperial: RectangularConcreteBeam) -> None:
    f = Forces(V_z=37.727*kip, N_x = 0*kip )
    A_s=0.847*inch**2 
    beam_example_imperial.design_shear(f, A_s) 
    best_combination = beam_example_imperial.transverse_rebar
    assert best_combination is not None
    assert best_combination['n_stirrups'] == 1
    assert best_combination['d_b'] == 12*mm 
    assert best_combination['spacing'] ==  7*inch 
    assert best_combination['A_v'].mangnitude == pytest.approx(12.72,rel=1e-3)

# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()