import pytest
from mento.material import Concrete_ACI_318_19, Concrete_EN_1992, Concrete_EHE_08, SteelBar, SteelStrand
from mento.units import MPa, m, kg, GPa


# Get only the value
# debug(N.magnitude, N.to('kN').magnitude)

def test_concrete_aci_318_19_properties() -> None:
    f_c = 25 * MPa
    concrete = Concrete_ACI_318_19(name="C25", f_c=f_c)

    assert concrete.f_c == f_c
    assert concrete.design_code == "ACI 318-19" 
    assert concrete.density == 2500 * kg / m**3 
    assert concrete.E_c.to('MPa').magnitude == pytest.approx(26875, rel=0.01)
    assert concrete.f_r.to('MPa').magnitude == pytest.approx(3.125, rel=0.01)   

def test_concrete_en_1992_properties() -> None:
    f_ck = 25 * MPa
    concrete = Concrete_EN_1992(name="C25/30", f_ck=f_ck)

    assert concrete.f_ck  == f_ck
    assert concrete.design_code == "EN 1992"
    assert concrete.density == 2500 * kg / m**3 
    assert concrete.E_cm.to('MPa').magnitude == pytest.approx(31475, rel=0.01)
    assert concrete.f_ctm.to('MPa').magnitude == pytest.approx(2.564, rel=0.01) 

def test_concrete_ehe_08_properties() -> None:
    f_ck = 25 * MPa
    concrete = Concrete_EHE_08(name="C25", f_ck=f_ck)

    assert concrete.f_ck == f_ck
    assert concrete.design_code == "EHE-08"
    assert concrete.density == 2500 * kg / m**3 
    assert concrete.f_ck == f_ck
    assert concrete.f_cm == f_ck + 8 * MPa
    assert concrete.E_cm.to('MPa').magnitude == pytest.approx(27264 , rel=0.01)
    assert concrete.f_ctm.to('MPa').magnitude == pytest.approx(2.56 , rel=0.01) 

def test_steel_bar_properties() -> None:
    f_y = 500 * MPa
    steelbar = SteelBar(name="ADN 500", f_y=f_y)

    assert steelbar.name == "ADN 500"
    assert steelbar.f_y == f_y
    assert steelbar.E_s == 200 * GPa
    assert steelbar.epsilon_y == pytest.approx(f_y.to('MPa') /( steelbar.E_s.to('MPa')), rel=0.01) 

def test_steel_strand_properties() -> None:
    f_y = 1700 * MPa
    steelstrand = SteelStrand(name='Y1860', f_y=f_y)

    assert steelstrand.name == "Y1860"
    assert steelstrand.f_y == f_y
    assert steelstrand.f_u == 1860 * MPa
    assert steelstrand.E_s == 190 * GPa

# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()