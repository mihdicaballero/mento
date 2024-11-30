
from devtools import debug

from mento import cm, MPa, m, kg, kN, kNm, inch, ft, kip, lb, psi, mm, ksi
from mento.settings import Settings
from mento.material import Concrete, SteelBar, Concrete_EN_1992_2004, Concrete_ACI_318_19
from mento.section import Section
from mento.rectangular import RectangularSection

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
    concrete = Concrete('C25') 
    steel_bar = SteelBar(name="ADN 420", f_y=60*ksi) 
    section = RectangularSection(concrete=concrete, steel_bar=steel_bar, width=10*inch, height=16*inch)
    debug(section.width, section.height)
    debug(section.A_x, section.I_y, section.I_z)
    debug(section.get_settings)

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

if __name__ == "__main__":
    # units()
    settings()
    # section()
    # rectangular()
    # material()