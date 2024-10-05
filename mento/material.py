from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, TYPE_CHECKING
import math
from mento.units import kg, m, MPa, ksi, GPa
from devtools import debug
from mento import ureg

# Conditional import for type checking only
if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

@dataclass
class Material:
    name: str  # Name will be provided by subclasses or instances

@dataclass
class Concrete(Material):
    f_c: PlainQuantity = field(default=25*MPa)
    density: PlainQuantity = 2500*kg/m**3
    design_code: str = field(default="ACI 318-19")

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = {
            'f_c': self.f_c.to('MPa'),
            'density': self.density.to('kg / meter ** 3')
        }
        return properties

@dataclass
class Concrete_ACI_318_19(Concrete):
    _E_c: PlainQuantity = field(init=False)
    _f_r: PlainQuantity = field(init=False)
    epsilon_c: float = field(default=0.003, init=False)
    _beta_1: float = field(init=False)

    def __post_init__(self) -> None:
        # Ensure name is properly set, either hardcode or pass during instantiation
        if not self.name:
            self.name = "Concrete ACI"  # You can set a default name here
        self.design_code = "ACI 318-19"
        self._E_c = ((self.density / (kg / m**3)) ** 1.5) * 0.043 * math.sqrt(self.f_c / MPa) * MPa
        self._f_r = 0.625 * math.sqrt(self.f_c / MPa) * MPa
        self._beta_1 = self.__beta_1()

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = super().get_properties()       
        properties['E_c'] = self._E_c.to('MPa')
        properties['f_r'] = self._f_r.to('MPa')
        properties['beta_1'] = ureg.Quantity(self._beta_1, '')
        properties['epsilon_c']=ureg.Quantity(self.epsilon_c, '')
        return properties
    
    def __beta_1(self) -> float:
        # Table 22.2.2.4.3—Values of β1 for equivalent rectangular concrete stress distribution
        # Page 399
        if 17 <= self.f_c / MPa <= 28:
            return 0.85
        elif 28 < self.f_c / MPa <= 55:
            return 0.85 - 0.05 / 7 * (self.f_c / MPa - 28)
        elif self.f_c / MPa > 55:
            return 0.65
        else:
            # Handle case where f_c / MPa < 17
            return 0.85  # or another appropriate value based on your requirements
        
    @property
    def E_c(self) -> PlainQuantity:
        return self._E_c

    @property
    def f_r(self) -> PlainQuantity:
        return self._f_r
    
    @property
    def beta_1(self) -> float:
        return self._beta_1
@dataclass
class Concrete_EN_1992(Concrete):
    _E_cm: PlainQuantity = field(init=False)  # Secant modulus of elasticity
    _f_ck: PlainQuantity = field(init=False) # Characteristic concrete strength
    _f_cm: PlainQuantity = field(init=False) # mean compressive strength
    _f_ctm: PlainQuantity = field(init=False) # Mean tensile strength

    def __init__(self, name: str, f_ck: PlainQuantity):
        super().__init__(name=name, f_c=f_ck)
        self.design_code = "EN 1992"
        self._f_ck = self.f_c
        self._f_cm = self._f_ck + 8 * MPa
        self._E_cm = 22000 * (self._f_cm / (10 * MPa)) ** 0.3 * MPa       
        self._f_ctm = 0.3 * (self._f_ck / MPa) ** (2/3) * MPa

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = super().get_properties()
        properties['E_cm'] = self._E_cm.to('MPa')
        properties['f_ctm'] = self._f_ctm.to('MPa')
        return properties

    def alpha_cc(self) -> float:
        # Example implementation for alpha_cc, as per Eurocode EN 1992-1-1
        return 1.0  # Typically, this value is taken as 1.0 for normal weight concrete
    
    @property
    def E_cm(self) -> PlainQuantity:
        return self._E_cm

    @property
    def f_ck(self) -> PlainQuantity:
        return self._f_ck
    @property
    def f_cm(self) -> PlainQuantity:
        return self._f_cm
    @property
    def f_ctm(self) -> PlainQuantity:
        return self._f_ctm

@dataclass
class Concrete_EHE_08(Concrete):
    _E_cm: PlainQuantity = field(init=False)  # Secant modulus of elasticity
    _f_ck: PlainQuantity = field(init=False) # Characteristic concrete strength
    _f_cm: PlainQuantity = field(init=False) # mean compressive strength
    _f_ctm: PlainQuantity = field(init=False) # Mean tensile strength
    _f_ctm_fl: PlainQuantity = field(init=False) # Mean flexure tensile strength

    def __init__(self, name: str, f_ck: PlainQuantity):
        super().__init__(name=name, f_c=f_ck)
        self.design_code = "EHE-08"
        # Calculate _E_cm and f_ctk based on f_c
        # Formulas are based on EHE-08 specifications
        self._f_ck = self.f_c
        self._f_cm = self._f_ck+8*MPa
        self._E_cm = 8500 *(self._f_cm / MPa)**(1/3) * MPa 
        self._f_ctm = 0.3 * (self._f_ck / MPa)**(2/3) * MPa 

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = super().get_properties()
        properties['E_cm'] = self._E_cm.to('MPa')
        properties['f_ck'] = self._f_ck.to('MPa')
        properties['f_cm'] = self._f_cm.to('MPa')
        properties['f_ctm'] = self._f_ctm.to('MPa')
        return properties

    def alpha_cc(self) -> float:
        # Example implementation for alpha_cc as per EHE-08
        return 0.85  # Example value, modify according to EHE-08 standards
    
    @property
    def E_cm(self) -> PlainQuantity:
        return self._E_cm

    @property
    def f_ck(self) -> PlainQuantity:
        return self._f_ck
    @property
    def f_cm(self) -> PlainQuantity:
        return self._f_cm
    @property
    def f_ctm(self) -> PlainQuantity:
        return self._f_ctm

# # Factory function
# def create_concrete(name: str, f_c: Quantity, design_code: str) -> Concrete:
#     if design_code == "ACI 318-19":
#         return Concrete_ACI_318_19(name=name, f_c=f_c)
#     elif design_code == "EN 1992":
#         return Concrete_EN_1992(name=name, f_c=f_c)
#     elif design_code == "EHE-08":
#         return Concrete_EHE_08(name=name, f_c=f_c)
#     else:
#         raise ValueError(f"Invalid design code: {design_code}. Options: ACI 318-19, EN 1992, EHE-08.")
@dataclass
class Steel(Material):
    _f_y: PlainQuantity = field(init=False)
    _density: PlainQuantity = field(default=7850 * kg / m**3) 

    def __init__(self, name: str, f_y: PlainQuantity, density: PlainQuantity = 7850 * kg / m**3):
        super().__init__(name)
        self._f_y = f_y
        self._density = density

    @property
    def f_y(self) -> PlainQuantity:
        return self._f_y
    
    @property
    def density(self) -> PlainQuantity:
        return self._density
    
    # Maximum f_yt for Imperial system
    @property 
    def f_yt(self) -> PlainQuantity:
        return min(self.f_y, 60*ksi)

@dataclass
class SteelBar(Steel):
    _E_s: PlainQuantity = field(default=200*GPa)
    _epsilon_y: PlainQuantity = field(init=False)

    def __init__(self, name: str, f_y: PlainQuantity, density: PlainQuantity =7850 *kg/m**3):
        super().__init__(name, f_y, density)
        self._epsilon_y = self._f_y / self._E_s # 21.2.2.1 - Page 392

    @property
    def E_s(self) -> PlainQuantity:
        return self._E_s

    @property
    def epsilon_y(self) -> PlainQuantity:
        return self._epsilon_y

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = {
            'E_s': self._E_s.to('GPa'),
            'f_y': self._f_y.to('MPa'),
            'epsilon_y': self._epsilon_y
            }
        return properties

@dataclass
class SteelStrand(Steel):
    _f_u: PlainQuantity = field(default=1860*MPa)
    _E_s: PlainQuantity = field(default=190000*MPa)
    prestress_stress: PlainQuantity = field(default=0*MPa)

    def __init__(self, name: str, f_y: PlainQuantity, density: PlainQuantity =7850 *kg/m**3):
        super().__init__(name, f_y, density)
        self._epsilon_y = self._f_y / self._E_s


    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = {
            'E_s': self._E_s.to('MPa'),
            'f_y': self._f_y.to('MPa'),
            'f_u': self._f_u.to('MPa')
        }
        return properties
    
    @property
    def f_u(self) -> PlainQuantity:
        return self._f_u
    
    @property
    def E_s(self) -> PlainQuantity:
        return self._E_s

def main() -> None:
    # Test cases
    concrete = Concrete_ACI_318_19(name="H25",f_c=25*MPa)
    debug(concrete.name, concrete.design_code)
    debug(concrete.get_properties())
    steelbar = SteelBar(name="ADN 500",f_y=500*MPa)
    debug(steelbar.get_properties())
    steelstrand = SteelStrand(name='Y1860',f_y=1700*MPa)
    debug(steelstrand.get_properties())
    print(concrete.f_c.to('MPa'), concrete.f_c.to('MPa').magnitude)

if __name__ == "__main__":
    main()


