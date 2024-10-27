from mento.concrete.beam import RectangularConcreteBeam
from mento.material import Concrete_ACI_318_19, SteelBar
from mento.units import psi, kip, inch, ksi, ft
from mento.forces import Forces
import pytest
import numpy as np

@pytest.fixture()
def beam_example() -> RectangularConcreteBeam:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000*psi)  
    steelBar = SteelBar(name="ADN 420", f_y=60*ksi)  
    custom_settings = {'clear_cover': 1.5*inch, 'stirrup_diameter_ini':0.5*inch,
                       'longitudinal_diameter_ini': 1*inch} 
    section = RectangularConcreteBeam(
        label="V-10x16",
        concrete=concrete,
        steel_bar=steelBar,
        width=10*inch,  
        height=16*inch,
        settings=custom_settings  
    )
    return section

def test_shear_check_ACI_318_19(beam_example: RectangularConcreteBeam) -> None:
    f = Forces(V_z=37.727*kip, N_x = 0*kip)  
    A_s = 0.847*inch**2 
    beam_example.set_transverse_rebar(n_stirrups=1, d_b=0.5*inch, s_l=6*inch) 
    results = beam_example.check_shear_ACI_318_19(f, A_s)  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(10.419, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(16.624, rel=1e-3)
    assert results.iloc[0]['ØVc'].magnitude  == pytest.approx(56.969, rel=1e-3)
    assert results.iloc[0]['ØVs'].magnitude  == pytest.approx(176.865, rel=1e-3)
    assert results.iloc[0]['ØVn'].magnitude  == pytest.approx(233.834, rel=1e-3)
    assert results.iloc[0]['ØVmax'].magnitude  == pytest.approx(284.847, rel=1e-3)
    assert results.iloc[0]["FUv"].magnitude  == pytest.approx(0.7177, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['Vu<ØVmax'] is np.True_
    assert results.iloc[0]['Vu<ØVn'] is np.True_

def test_shear_design_ACI_318_19(beam_example: RectangularConcreteBeam) -> None:
    f = Forces(V_z=37.727*kip, N_x = 0*kip) 
    A_s = 0.847*inch**2  
    results = beam_example.design_shear_ACI_318_19(f, A_s)  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(10.419, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(10.472, rel=1e-3)
    assert results.iloc[0]['ØVc'].magnitude  == pytest.approx(56.969, rel=1e-3)
    assert results.iloc[0]['ØVs'].magnitude  == pytest.approx(112.288, rel=1e-3)
    assert results.iloc[0]['ØVn'].magnitude  == pytest.approx(169.258, rel=1e-3)
    assert results.iloc[0]['ØVmax'].magnitude  == pytest.approx(284.847, rel=1e-3)
    assert results.iloc[0]["FUv"].magnitude  == pytest.approx(0.9915, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['Vu<ØVmax'] is np.True_
    assert results.iloc[0]['Vu<ØVn'] is np.True_

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