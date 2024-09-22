import sys
import os

# Add the project root to PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
from structurelab.material import Concrete_ACI_318_19, Concrete_EN_1992, Concrete_EHE_08, SteelBar, SteelStrand
from structurelab.units import MPa, m, kg


def test_concrete_aci_318_19_properties():
    f_c = 25 * MPa
    concrete = Concrete_ACI_318_19(name="H30", f_c=f_c)
    properties = concrete.get_properties()

    assert properties['f_c'] == f_c
    assert properties['design_code'] == "ACI 318-19" 
    assert properties['density'] == 2500 * kg / m**3 
    assert properties['E_c'].value == pytest.approx(26875500000, rel=0.01) # Example value, replace with the correct value 
    assert properties['f_r'].value == pytest.approx(3125000, rel=0.01)   # Example value, replace with the correct value  

def test_concrete_en_1992_properties():
    f_c = 25 * MPa
    concrete = Concrete_EN_1992(name="C25/30", f_c=f_c, design_code="EN 1992")
    properties = concrete.get_properties()

    assert properties['f_c'] == f_c
    assert properties['design_code'] == "EN 1992"
    assert properties['density'] == 2500 * kg / m**3 
    assert properties['E_cm'].value == pytest.approx(31475000000, rel=0.01)
    assert properties['f_ctm'].value == pytest.approx(2564963, rel=0.01)  # Example value, replace with the correct value 

def test_concrete_ehe_08_properties():
    f_c = 25 * MPa
    concrete = Concrete_EHE_08(name="C25/30", f_c=f_c, design_code="EHE-08")
    properties = concrete.get_properties()

    assert properties['f_c'] == f_c
    assert properties['design_code'] == "EHE-08"
    assert properties['density'] == 2500 * kg / m**3 
    assert properties['f_ck'] == f_c
    assert properties['f_cm'] == f_c + 8 * MPa
    assert properties['E_cm'].value == pytest.approx(27264000000 , rel=0.01) # Example value, replace with the correct value 
    assert properties['f_ctm'].value == pytest.approx(2560000 , rel=0.01) # Example value, replace with the correct value  

def test_steel_bar_properties():
    f_y = 500 * MPa
    steelbar = SteelBar(name="ADN 500", f_y=f_y)
    properties = steelbar.get_properties()

    assert properties['name'] == "ADN 500"
    assert properties['f_y'] == f_y
    assert properties['E_s'] == 200000 * MPa
    assert steelbar._epsilon_y == pytest.approx(f_y /( properties['E_s']), rel=0.01) 

def test_steel_strand_properties():
    f_y = 1700 * MPa
    steelstrand = SteelStrand(name='Y1860', f_y=f_y)
    properties = steelstrand.get_properties()

    assert properties['name'] == "Y1860"
    assert properties['f_y'] == f_y
    assert properties['f_u'] == 1860 * MPa
    assert properties['E_s'] == 190000 * MPa

# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()