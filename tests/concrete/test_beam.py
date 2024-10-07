from mento.concrete.beam import RectangularConcreteBeam
from mento.material import Concrete_ACI_318_19, SteelBar
from mento.units import psi, kip, inch, ksi, ft
from mento.forces import Forces
import pytest

@pytest.fixture()
def beam_example() -> RectangularConcreteBeam:
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

def test_shear_check_ACI_318_19(beam_example: RectangularConcreteBeam) -> None:
    f = Forces(V_z=37.727*kip, N_x = 0*kip)  
    A_s = 0.847*inch**2  
    results = beam_example.check_shear_ACI_318_19(f, A_s, d_b=0.5*inch, s=6*inch, n_legs=2)  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results['A_v_min'].magnitude  == pytest.approx(2.117, rel=1e-3)
    assert results['A_v_req'].magnitude  == pytest.approx(9.433, rel=1e-3)
    assert results['A_v'].magnitude  == pytest.approx(16.624, rel=1e-3)
    assert results['phi_V_c'].magnitude  == pytest.approx(60.767, rel=1e-3)
    assert results['phi_V_s'].magnitude  == pytest.approx(188.656, rel=1e-3)
    assert results['phi_V_n'].magnitude  == pytest.approx(249.423, rel=1e-3)
    assert results['phi_V_max'].magnitude  == pytest.approx(303.876, rel=1e-3)
    assert results["FUv"].magnitude  == pytest.approx(0.6728, rel=1e-3)

    # Assert non-numeric values directly
    assert results['shear_ok'] is True
    assert results['max_shear_ok'] is True

def test_shear_design_ACI_318_19(beam_example: RectangularConcreteBeam) -> None:
    f = Forces(V_z=37.727*kip, N_x = 0*kip) 
    A_s = 0.847*inch**2  
    results = beam_example.design_shear_ACI_318_19(f, A_s)  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results['A_v_min'].magnitude  == pytest.approx(2.116, rel=1e-3)
    assert results['A_v_req'].magnitude  == pytest.approx(9.4333, rel=1e-3)
    assert results['A_v'].magnitude  == pytest.approx(12.722, rel=1e-3)
    assert results['phi_V_c'].magnitude  == pytest.approx(60.767, rel=1e-3)
    assert results['phi_V_s'].magnitude  == pytest.approx(144.370, rel=1e-3)
    assert results['phi_V_n'].magnitude  == pytest.approx(205.137, rel=1e-3)
    assert results['phi_V_max'].magnitude  == pytest.approx(303.837, rel=1e-3)
    assert results["FUv"]  == pytest.approx(0.8181, rel=1e-3)

    # Assert non-numeric values directly
    assert results['shear_ok'] is True
    assert results['max_shear_ok'] is True

# ------- Flexural test --------------
# LO QUE ME ESTÁ FALTANDO ES SETEAR EL MEC COVER.
# Example 6.6 CRSI Design Guide
concreteFlexureTest1 = Concrete_ACI_318_19(name="fc4000",f_c=4000*psi)
steelBarFlexureTest1 = SteelBar(name="G60", f_y=60000*psi) 
sectionFlexureTest1 = RectangularConcreteBeam(
        label="B-12x24",
        concrete=concreteFlexureTest1,
        steel_bar=steelBarFlexureTest1,
        width=12 * inch, 
        height=24 * inch, 
    )

def test_flexural_design() -> None:
    results = sectionFlexureTest1.design_flexure(258.3*kip*ft) 
    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results['As_min_code'].magnitude  == pytest.approx(0.0005548, rel=1e-3)
    assert results['As_required'].magnitude  == pytest.approx(0.0019173, rel=1e-3)
    assert results['As_max'].magnitude  == pytest.approx(0.0030065, rel=1e-3)
    assert results['As_adopted'].magnitude  == pytest.approx(0.0019173, rel=1e-3)
    assert results['As_compression'].magnitude  == pytest.approx(0, rel=1e-3)


# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()