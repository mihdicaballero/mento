from mento.concrete.beam import Beam
from mento import material

import pytest

@pytest.fixture()
def beam_example() -> Beam:
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

def test_shear_check_ACI_318_19(beam_example: Beam) -> None:
    V_u = 37.727*kip  
    N_u = 0*kip  
    A_s = 0.847*inch**2  
    results = beam_example.check_shear_ACI_318_19(V_u, N_u, A_s, d_b=0.5*inch, s=6*inch, n_legs=2)  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results['A_v_min'].value  == pytest.approx(0.0002116, rel=1e-3)
    assert results['A_v_req'].value  == pytest.approx(0.00094333, rel=1e-3)
    assert results['A_v'].value  == pytest.approx(0.001662, rel=1e-3)
    assert results['phi_V_c'].value  == pytest.approx(60767.33, rel=1e-3)
    assert results['phi_V_s'].value  == pytest.approx(188655, rel=1e-3)
    assert results['phi_V_n'].value  == pytest.approx(249423, rel=1e-3)
    assert results['phi_V_max'].value  == pytest.approx(303876, rel=1e-3)
    assert results["FUv"]  == pytest.approx(0.6728, rel=1e-3)

    # Assert non-numeric values directly
    assert results['shear_ok'] is True
    assert results['max_shear_ok'] is True

def test_shear_design_ACI_318_19(beam_example: Beam) -> None:
    V_u = 37.727*kip  
    N_u = 0*kip  
    A_s = 0.847*inch**2  
    results = beam_example.design_shear_ACI_318_19(V_u, N_u, A_s)  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results['A_v_min'].value  == pytest.approx(0.0002116, rel=1e-3)
    assert results['A_v_req'].value  == pytest.approx(0.00094333, rel=1e-3)
    assert results['A_v'].value  == pytest.approx(0.0012722, rel=1e-3)
    assert results['phi_V_c'].value  == pytest.approx(60767.33, rel=1e-3)
    assert results['phi_V_s'].value  == pytest.approx(144370, rel=1e-3)
    assert results['phi_V_n'].value  == pytest.approx(205137.8, rel=1e-3)
    assert results['phi_V_max'].value  == pytest.approx(303836.7, rel=1e-3)
    assert results["FUv"]  == pytest.approx(0.8181, rel=1e-3)

    # Assert non-numeric values directly
    assert results['shear_ok'] is True
    assert results['max_shear_ok'] is True

# ------- Flexural test --------------
# LO QUE ME ESTÁ FALTANDO ES SETEAR EL MEC COVER.
# Example 6.6 CRSI Design Guide
concreteFlexureTest1=material.create_concrete(name="fc4000",f_c=4000*psi, design_code="ACI 318-19") # type: ignore
steelBarFlexureTest1=material.SteelBar(name="G60", f_y=60000*psi) # type: ignore
sectionFlexureTest1 = Beam(
        name="B-12x24",
        concrete=concreteFlexureTest1,
        steelBar=steelBarFlexureTest1,
        width=12 * inch,  # type: ignore
        depth=24 * inch,  # type: ignore
    )

def test_flexural_design():
    results = sectionFlexureTest1.design_flexure(258.3*kip*ft)  # type: ignore
    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results['As_min_code'].value  == pytest.approx(0.0005548, rel=1e-3)
    assert results['As_required'].value  == pytest.approx(0.0019173, rel=1e-3)
    assert results['As_max'].value  == pytest.approx(0.0030065, rel=1e-3)
    assert results['As_adopted'].value  == pytest.approx(0.0019173, rel=1e-3)
    assert results['As_compression'].value  == pytest.approx(0, rel=1e-3)


# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()