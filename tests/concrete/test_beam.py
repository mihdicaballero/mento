from mento.concrete.beam import Beam
from mento import material
from mento.units import ksi, psi, kip, inch
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
        depth=16*inch,
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
    assert results['A_v'].value  == pytest.approx(0.0012722, rel=1e-3)
    assert results['phi_V_c'].value  == pytest.approx(56969.4, rel=1e-3)
    assert results['phi_V_s'].value  == pytest.approx(176865, rel=1e-3)
    assert results['phi_V_n'].value  == pytest.approx(233834, rel=1e-3)
    assert results['phi_V_max'].value  == pytest.approx(284847, rel=1e-3)
    assert results["FUv"]  == pytest.approx(0.7176789, rel=1e-3)

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
    assert results['phi_V_n'].value  == pytest.approx(214875, rel=1e-3)
    assert results['phi_V_max'].value  == pytest.approx(284847, rel=1e-3)
    assert results["FUv"]  == pytest.approx(0.7810047891, rel=1e-3)

    # Assert non-numeric values directly
    assert results['shear_ok'] is True
    assert results['max_shear_ok'] is True


# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()