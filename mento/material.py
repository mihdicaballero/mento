from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, TYPE_CHECKING
import math
from mento.units import kg, m, MPa, ksi, GPa, psi, Pa, lb, ft
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
    design_code: str = field(default="ACI 318-19")
    unit_system: str = field(init=False)

    def __post_init__(self) -> None:
        # Detect the unit system based on f_c
        if self.f_c.units == MPa or self.f_c.units == Pa:
            self.unit_system = "metric"
            self.density: PlainQuantity = 2500*kg/m**3
        elif self.f_c.units == psi or self.f_c.units == ksi:
            self.unit_system = "imperial"
            self.density = 155*lb/ft**3
        else:
            raise ValueError("Unsupported unit system for f_c. Please use MPa or ksi.")

    def get_properties(self) -> Dict[str, PlainQuantity]:
        # Return properties in the appropriate unit system
        properties = {
            'f_c': self.f_c,
            'density': self.density
        }
        return properties

@dataclass
class Concrete_ACI_318_19(Concrete):
    _E_c: PlainQuantity = field(init=False)
    _f_r: PlainQuantity = field(init=False)
    epsilon_c: float = field(default=0.003, init=False)
    _beta_1: float = field(init=False)

    def __post_init__(self) -> None:
        super().__post_init__()
        # Ensure name is properly set, either hardcode or pass during instantiation
        if not self.name:
            self.name = "Concrete ACI"  # You can set a default name here
        self.design_code = "ACI 318-19"
        # Adjust calculations based on unit system
        if self.unit_system == "metric":
            self._E_c = ((self.density / (kg / m**3)) ** 1.5) * 0.043 * math.sqrt(self.f_c / MPa) * MPa
            self._f_r = 0.625 * math.sqrt(self.f_c / MPa) * MPa
        else:  # imperial
            self._E_c = ((self.density / (lb / ft**3)) ** 1.5) * 33 * math.sqrt(self.f_c / psi) * psi
            self._f_r = 7.5 * math.sqrt(self.f_c / psi) * psi
        self._beta_1 = self.__beta_1()

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = super().get_properties()       
        properties['E_c'] = self._E_c
        properties['f_r'] = self._f_r
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
    _f_ck: PlainQuantity = field(init=False)  # Characteristic concrete strength
    _f_cm: PlainQuantity = field(init=False)  # Mean compressive strength
    _f_ctm: PlainQuantity = field(init=False)  # Mean tensile strength
    _f_ctk: PlainQuantity = field(init=False)  # Characteristic tensile strength
    _f_cd: PlainQuantity = field(init=False)  # Design strength of concrete
    _f_1cd: PlainQuantity = field(init=False)  # Design strength for shear
    _f_ctd: PlainQuantity = field(init=False)  # Tensile strength of concrete

    def __init__(self, name: str, f_ck: PlainQuantity):
        super().__init__(name=name, f_c=f_ck)
        gamma_c: float = 1.5
        self.design_code: str = "EHE-08"
        
        # Assign provided characteristic concrete strength
        self._f_ck = self.f_c
        
        # Calculate mean compressive strength, tensile strength, and secant modulus of elasticity
        self._f_cm = self._f_ck + 8 * MPa
        self._E_cm = 8500 * (self._f_cm / MPa) ** (1 / 3) * MPa
        self._f_ctm = 0.3 * (self._f_ck / MPa) ** (2 / 3) * MPa
        self._f_ctk = 0.7 * self._f_ctm
        
        # Calculate design strengths based on partial safety factor gamma_c
        self._f_cd = self._f_ck / gamma_c
        self._f_1cd = (
            0.6 * self._f_cd if self._f_ck <= 60 * MPa 
            else max((0.9 - self._f_ck / (200 * MPa)) * self._f_cd, 0.5 * self._f_cd)
        )
        self._f_ctd = self._f_ctk / gamma_c

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = super().get_properties()
        properties.update({
            'E_cm': self._E_cm.to('MPa'),
            'f_ck': self._f_ck.to('MPa'),
            'f_cm': self._f_cm.to('MPa'),
            'f_ctm': self._f_ctm.to('MPa'),
            'f_ctk': self._f_ctk.to('MPa'),
            'f_cd': self._f_cd.to('MPa'),
            'f_1cd': self._f_1cd.to('MPa'),
            'f_ctd': self._f_ctd.to('MPa')
        })
        return properties

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

    @property
    def f_ctk(self) -> PlainQuantity:
        return self._f_ctk

    @property
    def f_cd(self) -> PlainQuantity:
        return self._f_cd

    @property
    def f_1cd(self) -> PlainQuantity:
        return self._f_1cd

    @property
    def f_ctd(self) -> PlainQuantity:
        return self._f_ctd

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

@dataclass
class SteelBar(Steel):
    _E_s: PlainQuantity = field(default=200*GPa)
    _epsilon_y: PlainQuantity = field(init=False)

    def __init__(self, name: str, f_y: PlainQuantity, density: PlainQuantity =7850 *kg/m**3):
        super().__init__(name, f_y, density)
        self._epsilon_y = f_y.to('MPa') /(self._E_s.to('MPa')) # 21.2.2.1 - Page 392

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


