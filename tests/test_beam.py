from mento.beam import RectangularBeam
from mento.node import Node
from mento.material import Concrete_ACI_318_19, SteelBar, Concrete_EN_1992_2004
from mento.units import psi, kip, inch, ksi, mm, kN, cm, MPa, ft
from mento.forces import Forces
import pytest
import numpy as np

@pytest.fixture()
def beam_example_imperial() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000*psi)  
    steelBar = SteelBar(name="ADN 420", f_y=60*ksi)  
    custom_settings = {'clear_cover': 1.5*inch} 
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
def beam_example_EN_1992_2004() -> RectangularBeam:
    concrete= Concrete_EN_1992_2004(name="C25",f_ck=25*MPa) 
    steelBar= SteelBar(name="B500S", f_y=500*MPa)
    custom_settings = {'clear_cover': 2.6*cm}
    section = RectangularBeam(label="V-20x60",
                                      concrete=concrete,steel_bar=steelBar,width=20*cm, height=60*cm,
                                       settings=custom_settings)
    return section

def test_shear_check_EN_1992_2004_rebar_1(beam_example_EN_1992_2004: RectangularBeam) -> None:
    f = Forces(V_z=100*kN)
    beam_example_EN_1992_2004.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=25*cm) 
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16*mm)
    Node(section=beam_example_EN_1992_2004, forces=f)
    results = beam_example_EN_1992_2004.check_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(1.825, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(2.262, rel=1e-3)
    assert results.iloc[0]['VEd,1'].magnitude  == pytest.approx(100, rel=1e-3)
    assert results.iloc[0]['VEd,2'].magnitude  == pytest.approx(100, rel=1e-3)
    assert results.iloc[0]['VRd,c'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['VRd,s'].magnitude  == pytest.approx(123.924, rel=1e-3)
    assert results.iloc[0]['VRd'].magnitude  == pytest.approx(123.924, rel=1e-3)
    assert results.iloc[0]['VRd,max'].magnitude  == pytest.approx(312.811, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(0.8069, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['VEd,1<VRd,max'] is np.True_
    assert results.iloc[0]['VEd,2<VRd'] is np.True_

def test_shear_check_EN_1992_2004_rebar_2(beam_example_EN_1992_2004: RectangularBeam) -> None:
    f = Forces(V_z=350*kN)
    beam_example_EN_1992_2004.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=25*cm)
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16*mm)
    Node(section=beam_example_EN_1992_2004, forces=f)
    results = beam_example_EN_1992_2004.check_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(7.533, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(2.262, rel=1e-3)
    assert results.iloc[0]['VEd,1'].magnitude  == pytest.approx(350, rel=1e-3)
    assert results.iloc[0]['VEd,2'].magnitude  == pytest.approx(350, rel=1e-3)
    assert results.iloc[0]['VRd,c'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['VRd,s'].magnitude  == pytest.approx(105.099, rel=1e-3)
    assert results.iloc[0]['VRd'].magnitude  == pytest.approx(105.099, rel=1e-3)
    assert results.iloc[0]['VRd,max'].magnitude  == pytest.approx(350, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(3.33, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['VEd,1<VRd,max'] is np.True_
    assert results.iloc[0]['VEd,2<VRd'] is np.False_

def test_shear_check_EN_1992_2004_rebar_3(beam_example_EN_1992_2004: RectangularBeam) -> None:
    f = Forces(V_z=500*kN)
    beam_example_EN_1992_2004.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=25*cm)
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16*mm)
    Node(section=beam_example_EN_1992_2004, forces=f)
    results = beam_example_EN_1992_2004.check_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(22.817, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(2.26, rel=1e-3)
    assert results.iloc[0]['VEd,1'].magnitude  == pytest.approx(500, rel=1e-3)
    assert results.iloc[0]['VEd,2'].magnitude  == pytest.approx(500, rel=1e-3)
    assert results.iloc[0]['VRd,c'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['VRd,s'].magnitude  == pytest.approx(49.566, rel=1e-3)
    assert results.iloc[0]['VRd'].magnitude  == pytest.approx(49.566, rel=1e-3)
    assert results.iloc[0]['VRd,max'].magnitude  == pytest.approx(453.6, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(10.088, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['VEd,1<VRd,max'] is np.False_
    assert results.iloc[0]['VEd,2<VRd'] is np.False_

def test_shear_check_EN_1992_2004_no_rebar_1(beam_example_EN_1992_2004: RectangularBeam) -> None:
    f = Forces(V_z=30*kN)
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16*mm)
    Node(section=beam_example_EN_1992_2004, forces=f)
    results = beam_example_EN_1992_2004.check_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['VEd,1'].magnitude  == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]['VEd,2'].magnitude  == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]['VRd,c'].magnitude  == pytest.approx(56.126, rel=1e-3)
    assert results.iloc[0]['VRd,s'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['VRd'].magnitude  == pytest.approx(56.126, rel=1e-3)
    assert results.iloc[0]['VRd,max'].magnitude  == pytest.approx(56.126, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(0.5345, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['VEd,1<VRd,max'] is np.True_
    assert results.iloc[0]['VEd,2<VRd'] is np.True_

def test_shear_check_EN_1992_2004_no_rebar_2(beam_example_EN_1992_2004: RectangularBeam) -> None:
    f = Forces(V_z=30*kN)
    Node(section=beam_example_EN_1992_2004, forces=f)
    results = beam_example_EN_1992_2004.check_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['VEd,1'].magnitude  == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]['VEd,2'].magnitude  == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]['VRd,c'].magnitude  == pytest.approx(39.681, rel=1e-3)
    assert results.iloc[0]['VRd,s'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['VRd'].magnitude  == pytest.approx(39.681, rel=1e-3)
    assert results.iloc[0]['VRd,max'].magnitude  == pytest.approx(39.681, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(0.756, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['VEd,1<VRd,max'] is np.True_
    assert results.iloc[0]['VEd,2<VRd'] is np.True_

def test_shear_check_EN_1992_2004_no_rebar_3(beam_example_EN_1992_2004: RectangularBeam) -> None:
    f = Forces(N_x=50*kN, V_z=30*kN)
    beam_example_EN_1992_2004.set_longitudinal_rebar_bot(n1=4, d_b1=16*mm)
    Node(section=beam_example_EN_1992_2004, forces=f)
    results = beam_example_EN_1992_2004.check_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(1.6, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['VEd,1'].magnitude  == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]['VEd,2'].magnitude  == pytest.approx(30, rel=1e-3)
    assert results.iloc[0]['VRd,c'].magnitude  == pytest.approx(63.101, rel=1e-3)
    assert results.iloc[0]['VRd,s'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['VRd'].magnitude  == pytest.approx(63.101, rel=1e-3)
    assert results.iloc[0]['VRd,max'].magnitude  == pytest.approx(63.101, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(0.4754, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['VEd,1<VRd,max'] is np.True_
    assert results.iloc[0]['VEd,2<VRd'] is np.True_

def test_shear_check_ACI_318_19_1(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=37.727*kip, N_x = 0*kip)  
    beam_example_imperial.set_transverse_rebar(n_stirrups=1, d_b=0.5*inch, s_l=6*inch)
    Node(section=beam_example_imperial, forces=f)
    results = beam_example_imperial.check_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(10.0623, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(16.624, rel=1e-3)
    assert results.iloc[0]['ØVc'].magnitude  == pytest.approx(58.288, rel=1e-3)
    assert results.iloc[0]['ØVs'].magnitude  == pytest.approx(180.956, rel=1e-3)
    assert results.iloc[0]['ØVn'].magnitude  == pytest.approx(239.247, rel=1e-3)
    assert results.iloc[0]['ØVmax'].magnitude  == pytest.approx(291.44, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(0.70144, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['Vu<ØVmax'] is np.True_
    assert results.iloc[0]['Vu<ØVn'] is np.True_

def test_shear_check_ACI_318_19_2(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=37.727*kip, N_x = 20*kip)  
    beam_example_imperial.set_transverse_rebar(n_stirrups=1, d_b=0.5*inch, s_l=6*inch)
    Node(section=beam_example_imperial, forces=f)
    results = beam_example_imperial.check_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(9.1803, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(16.624, rel=1e-3)
    assert results.iloc[0]['ØVc'].magnitude  == pytest.approx(67.888, rel=1e-3)
    assert results.iloc[0]['ØVs'].magnitude  == pytest.approx(180.959, rel=1e-3)
    assert results.iloc[0]['ØVn'].magnitude  == pytest.approx(248.847, rel=1e-3)
    assert results.iloc[0]['ØVmax'].magnitude  == pytest.approx(301.041, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(0.6743, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['Vu<ØVmax'] is np.True_
    assert results.iloc[0]['Vu<ØVn'] is np.True_

def test_shear_check_ACI_318_19_no_rebar_1(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=8*kip, N_x = 0*kip)  
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625*inch)
    Node(section=beam_example_imperial, forces=f)
    results = beam_example_imperial.check_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['ØVc'].magnitude  == pytest.approx(35.1253, rel=1e-3)
    assert results.iloc[0]['ØVs'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['ØVn'].magnitude  == pytest.approx(35.125, rel=1e-3)
    assert results.iloc[0]['ØVmax'].magnitude  == pytest.approx(268.278, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(1.013, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['Vu<ØVmax'] is np.True_
    assert results.iloc[0]['Vu<ØVn'] is np.False_

def test_shear_check_ACI_318_19_no_rebar_2(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=6*kip, N_x = 0*kip)  
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625*inch)
    Node(section=beam_example_imperial, forces=f)
    results = beam_example_imperial.check_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['ØVc'].magnitude  == pytest.approx(58.288, rel=1e-3)
    assert results.iloc[0]['ØVs'].magnitude  == pytest.approx(0, rel=1e-3)
    assert results.iloc[0]['ØVn'].magnitude  == pytest.approx(58.288, rel=1e-3)
    assert results.iloc[0]['ØVmax'].magnitude  == pytest.approx(291.441, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(0.4579, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['Vu<ØVmax'] is np.True_
    assert results.iloc[0]['Vu<ØVn'] is np.True_

def test_shear_design_ACI_318_19(beam_example_imperial: RectangularBeam) -> None:
    f = Forces(V_z=37.727*kip, N_x = 0*kip) 
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625*inch)
    Node(section=beam_example_imperial, forces=f)
    results = beam_example_imperial.design_shear()  

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results.iloc[0]['Av,min'].magnitude  == pytest.approx(2.117, rel=1e-3)
    assert results.iloc[0]['Av,req'].magnitude  == pytest.approx(10.06, rel=1e-3)
    assert results.iloc[0]['Av'].magnitude  == pytest.approx(11.22, rel=1e-3)
    assert results.iloc[0]['ØVc'].magnitude  == pytest.approx(58.29, rel=1e-3)
    assert results.iloc[0]['ØVs'].magnitude  == pytest.approx(122.15, rel=1e-3)
    assert results.iloc[0]['ØVn'].magnitude  == pytest.approx(180.44, rel=1e-3)
    assert results.iloc[0]['ØVmax'].magnitude  == pytest.approx(291.44, rel=1e-3)
    assert results.iloc[0]["DCR"]  == pytest.approx(0.93, rel=1e-3)

    # Assert non-numeric values directly
    assert results.iloc[0]['Vu<ØVmax'] is np.True_
    assert results.iloc[0]['Vu<ØVn'] is np.True_

# # ------- Flexural test --------------

def test_flexural_beam_determine_nominal_moment_simple_reinf_ACI_318_19() -> None:
    '''
        Test para la funcion que determina el momento nominal de una seccion armada solo con refuerzo de traccion
        ver el ejemplo en:
        Example from https://www.google.com/search?q=ACI+318+19+nominal+moment+of+section&oq=ACI+318+19+nominal+moment+of+section&gs_lcrp=EgZjaHJvbWUyBggAEEUYOTIHCAEQIRigAdIBCTIyNzkzajBqN6gCALACAA&sourceid=chrome&ie=UTF-8#fpstate=ive&vld=cid:dab5f5e7,vid:m4H0QbGDYIg,st:0
    '''
    concrete = Concrete_ACI_318_19(name="C6",f_c=6000*psi) 
    steelBar = SteelBar(name="G80", f_y=80000*psi) 
    beam = RectangularBeam(
        label="B20x30",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * inch,  
        height=30 * inch,   
    )
    A_s=10.92 * inch**2
    d=27*inch
    result=beam._determine_nominal_moment_simple_reinf_ACI_318_19(A_s,d)

    # Compare dictionaries with a tolerance for floating-point values, in inch**2
    # Result: 1653.84kip.ft
    assert result.to(kip*ft).magnitude  == pytest.approx(1653.84, rel=1e-2)

def test_flexural_beam_determine_nominal_moment_double_reinf_ACI_318_19() -> None:
    '''
        Test para la funcion que determina el momento nominal de una seccion doblemente armada
        Example from https://www.youtube.com/watch?v=7BOuV1gCcgc
    '''
    concrete = Concrete_ACI_318_19(name="C4",f_c=4000*psi) 
    steelBar = SteelBar(name="G60", f_y=60000*psi) 
    beam = RectangularBeam(
       label="B14x27",
       concrete=concrete,
       steel_bar=steelBar,
       width=14 * inch,  
       height=27 * inch,   
    )
    A_s=6 * inch**2
    d=24*inch
    d_prime=2.5*inch
    A_s_prime=1.8*inch**2
    result=beam._determine_nominal_moment_double_reinf_ACI_318_19(A_s, d, d_prime, A_s_prime)
    assert result.to(kip*ft).magnitude  == pytest.approx(639.12, rel=1e-2)

#VIGA DEL CALPCAD DE JPR
@pytest.fixture()
def beam_example_flexure() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="fc 4000", f_c=4000*psi)  
    steelBar = SteelBar(name="fy 60000", f_y=60*ksi)  
    custom_settings = {'clear_cover': 1.5*inch} 
    section = RectangularBeam(
        label="B-12x24",
        concrete=concrete,
        steel_bar=steelBar,
        width=12*inch,  
        height=24*inch,
        settings=custom_settings  
    )
    return section


def test_design_flexure_ACI_318_19_1(beam_example_flexure: RectangularBeam) -> None:
    f = Forces(label='Test_01', M_y=400*kip*ft)
    Node(section=beam_example_flexure, forces=f)
    results = beam_example_flexure.design_flexure()

    assert beam_example_flexure._d_bot.to(inch).magnitude == pytest.approx(20.72, rel=1e-2)
    assert beam_example_flexure._d_b1_b.to(inch).magnitude == pytest.approx(1, rel=1e-2)
    assert beam_example_flexure._n1_b == 2
    assert beam_example_flexure._d_b2_b.to(inch).magnitude == pytest.approx(1, rel=1e-2)
    assert beam_example_flexure._n2_b == 2
    assert beam_example_flexure._d_b3_b.to(inch).magnitude == pytest.approx(1, rel=1e-2)
    assert beam_example_flexure._n3_b == 2
    assert beam_example_flexure._d_b4_b.to(inch).magnitude == pytest.approx(0.75, rel=1e-2)
    assert beam_example_flexure._n4_b == 1
    assert beam_example_flexure._d_top.to(inch).magnitude == pytest.approx(21.31, rel=1e-2)
    assert beam_example_flexure._d_b1_t.to(inch).magnitude == pytest.approx(0.5, rel=1e-2)
    assert beam_example_flexure._n1_t == 2
    assert beam_example_flexure._d_b2_t.to(inch).magnitude == pytest.approx(0.375, rel=1e-2)
    assert beam_example_flexure._d_b3_t.to(inch).magnitude == pytest.approx(0.375, rel=1e-2)
    assert beam_example_flexure._d_b4_t.to(inch).magnitude == pytest.approx(0.375, rel=1e-2)

    assert results.iloc[0]['Section Label'] == 'B-12x24'
    assert results.iloc[0]['Load Combo']  == 'Test_01'
    assert results.iloc[0]['Position'] == 'Bottom'
    assert results.iloc[0]['As,min'].to(cm**2).magnitude == pytest.approx(5.35, rel=1e-2)
    assert results.iloc[0]['As,req bot'].to(cm**2).magnitude == pytest.approx(32.65, rel=1e-3)
    assert results.iloc[0]['As,req top'].to(cm**2).magnitude == pytest.approx(3.92, rel=1e-2)
    assert results.iloc[0]['As'].to(inch**2).magnitude == pytest.approx(5.15, rel=1e-2)
    # assert results.iloc[0]['c/d'] == 
    # assert results.iloc[0]['Mu'] == 
    # assert results.iloc[0]['ØMn'] == 
    # assert results.iloc[0]['Mu<ØMn'] == 
    # assert results.iloc[0]['DCR'] == 
    # assert results.iloc[1]['Section Label'] == 'B-12x24'
    # assert results.iloc[1]['Load Combo']  == 'Test_01'
    # assert results.iloc[1]['Position'] == 'Top'
    # assert results.iloc[1]['As,min'] == 
    # assert results.iloc[1]['As,req'] == 
    # assert results.iloc[1]['As,top'].to(inch**2).magnitude == pytest.approx(0.99402, rel=1e-2)
    # assert results.iloc[1]['c/d'] == 
    # assert results.iloc[1]['Mu'] == 
    # assert results.iloc[1]['ØMn'] == 
    # assert results.iloc[1]['Mu<ØMn'] == 
    # assert results.iloc[1]['DCR'] == 





# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()