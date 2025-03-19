
from devtools import debug

from mento import cm, MPa, m, kg, kN, kNm, inch, ft, kip, lb, psi, mm, ksi
from mento.settings import Settings
from mento.material import Concrete, SteelBar, Concrete_EN_1992_2004, Concrete_ACI_318_19
from mento.section import Section
from mento.rectangular import RectangularSection
from mento import Forces, RectangularBeam

def units() -> None:
    debug(2*cm, 3*MPa, 4*kg, 1*mm, 1*m, 3*kN, 2*kNm)
    debug(1*psi, 1*lb, 1*kip, 1*psi, 1*ksi, 1*inch, 1*ft)
    N = 3*kN
    A = 3*m**2
    s = N/A
    a=3.25*m
    print(N.to_compact(), N.to('kip'),A, s)
    # Get only the value
    debug(N.magnitude, N.to('kN').magnitude)
    debug(N, A, s)

    # How to format specific output
    debug('Length is {:~P}'.format(a))

def settings() -> None:
    settings_test = Settings()
    # debug(settings_test.default_settings, settings_test.aci_318_19_settings)
    custom_settings = {'clear_cover': 50*mm, 'longitudinal_diameter_ini': 25*mm}
    settings_test.update(custom_settings)
    # debug(settings_test.default_settings)
    print(settings_test)

def section() -> None:
    concrete = Concrete('C25')
    steel_bar = SteelBar(name='B500S', f_y=500*MPa)
    section = Section(concrete= concrete, steel_bar=steel_bar) 
    debug(section.get_settings, section.id)
    custom_settings = {'clear_cover': 20}
    section.update_settings(custom_settings)
    debug(section.settings.default_settings)
    debug(section.get_settings)
    debug(section)
    debug(section.get_settings)
    custom_settings = {'clear_cover': 20}
    section.update_settings(custom_settings)
    debug(section.get_settings)

def rectangular() -> None:
    concrete = Concrete('C25', f_c=6*ksi) 
    steel_bar = SteelBar(name="ADN 420", f_y=60*ksi) 
    section = RectangularSection(concrete=concrete, steel_bar=steel_bar, width=12*inch, height=20*inch)
    # debug(section.width, section.height)
    section.plot()
    # debug(section.A_x, section.I_y, section.I_z)
    # debug(section.get_settings)

def material() -> None:
    # Test cases
    # concrete = Concrete_ACI_318_19(name="H25",f_c=4*ksi)
    concrete = Concrete_ACI_318_19(name="H25",f_c=25*MPa)
    debug(concrete.name, concrete.design_code)
    debug(concrete.get_properties())
    print(concrete)
    steelbar = SteelBar(name="ADN 500",f_y=500*MPa)
    debug(steelbar.get_properties())
    print(steelbar)
    # steelstrand = SteelStrand(name='Y1860',f_y=1700*MPa)
    # debug(steelstrand.get_properties())
    # print(concrete.f_c.to('MPa'), concrete.f_c.to('MPa').magnitude)
    # print(concrete.unit_system)
    # concrete = Concrete_EN_1992_2004(name="H25",f_ck=25*MPa)
    # debug(concrete.name, concrete.design_code)
    # debug(concrete.get_properties())
    # debug(concrete.f_ctm)

def shear_EN_1992() -> None:
    concrete= Concrete_EN_1992_2004(name="C25",f_ck=25*MPa) 
    steelBar= SteelBar(name="B500S", f_y=500*MPa)
    custom_settings = {'clear_cover': 2.5*cm}
    beam = RectangularBeam(label="101",
                                      concrete=concrete,steel_bar=steelBar,width=30*cm, height=60*cm,
                                       settings=custom_settings)
    # f = Forces(V_z=100*kN, M_y=100*kNm)
    f = Forces(V_z=30*kN, N_x=0*kN)
    forces=[f]
    beam.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=25*cm)
    beam.set_longitudinal_rebar_bot(n1=2,d_b1=32*mm, n2=3, d_b2=20*mm, n3=2, d_b3=10*mm, n4=2, d_b4=8*mm)
    beam.set_longitudinal_rebar_top(n1=2,d_b1=25*mm, n2=3, d_b2=16*mm, n3=2, d_b3=10*mm, n4=2, d_b4=8*mm)
    results = beam.check_shear(forces)
    #results = beam.design_shear(forces)
    print(results)
    beam.plot()

if __name__ == "__main__":
    # units()
    # settings()
    # section()
    # rectangular()
    # material()
    shear_EN_1992()