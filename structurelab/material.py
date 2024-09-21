# common/materials.py

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict
import math
import forallpeople
forallpeople.environment('structural',top_level=True)
#Useful extra units
cm = 1e-2*m # type: ignore

@dataclass
class Material:
    name: str

    def get_properties(self) -> Dict[str, float]:
        return {
            'name':self.name
        }
@dataclass
class Concrete(Material):
    f_c: float = field(default=25*MPa) # type: ignore
    design_code: str = field(default="ACI 318-19")
    density: float = 2500*kg/m**3 # type: ignore

    def get_properties(self) -> Dict[str, float]:
        properties = super().get_properties()
        properties['f_c'] = self.f_c
        properties['design_code'] = self.design_code
        properties['density'] = self.density
        return properties
@dataclass
class Concrete_ACI_318_19(Concrete):
    _E_c: float = field(init=False)
    _f_r: float = field(init=False)
    epsilon_c: float = field(default=0.003, init=False)
    _beta_1: float = field(init=False)

    def __post_init__(self):
        self._E_c =  ((self.density / (kg / m**3)) ** 1.5) * 0.043 * math.sqrt(self.f_c / MPa) * MPa # type: ignore
        self._f_r = 0.625*math.sqrt(self.f_c/MPa)*MPa # type: ignore
        self._beta_1=self.__beta_1()

    def get_properties(self) -> dict:
        properties = super().get_properties()       
        properties['E_c'] = self._E_c
        properties['f_r'] = self._f_r
        properties['beta_1'] = self._beta_1
        properties['epsilon_c']=self.epsilon_c
        return properties
    
    def __beta_1(self) -> float:
        # Table 22.2.2.4.3—Values of β1 for equivalent rectangular concrete stress distribution
        # Page 399
        if 17 <= self.f_c / MPa <= 28: # type: ignore
            return 0.85
        elif 28 < self.f_c / MPa <= 55: # type: ignore
            return 0.85 - 0.05/7 * (self.f_c / MPa - 28) # type: ignore
        elif 55 < self.f_c / MPa: # type: ignore
            return 0.65
        
    @property
    def E_c(self) -> float:
        return self._E_c

    @property
    def f_r(self) -> float:
        return self._f_r
    
    @property
    def beta_1(self) -> float:
        return self._beta_1
@dataclass
class Concrete_EN_1992(Concrete):
    _E_cm: float = field(init=False)  # Secant modulus of elasticity
    _f_ck: float = field(init=False) # Characteristic concrete strength
    _f_cm: float = field(init=False) # mean compressive strength
    _f_ctm: float = field(init=False) # Mean tensile strength

    def __post_init__(self):
        # Calculate _E_cm and f_ctm based on f_c and density
        self._f_ck = self.f_c
        self._f_cm = self._f_ck+8*MPa #type:ignore
        self._E_cm = 22000 * (self._f_cm / (10*MPa)) ** (0.3) * MPa # type: ignore       
        self._f_ctm = 0.3 * (self._f_ck / MPa) ** (2/3) * MPa # type: ignore

    def get_properties(self) -> dict:
        properties = super().get_properties()
        properties['E_cm'] = self._E_cm
        properties['f_ctm'] = self._f_ctm
        return properties

    def alpha_cc(self):
        # Example implementation for alpha_cc, as per Eurocode EN 1992-1-1
        return 1.0  # Typically, this value is taken as 1.0 for normal weight concrete

@dataclass
class Concrete_EHE_08(Concrete):
    _E_cm: float = field(init=False)  # Secant modulus of elasticity
    _f_ck: float = field(init=False) # Characteristic concrete strength
    _f_cm: float = field(init=False) # mean compressive strength
    _f_ctm: float = field(init=False) # Mean tensile strength
    _f_ctm_fl: float = field(init=False) # Mean flexure tensile strength

    def __post_init__(self):
        # Calculate _E_cm and f_ctk based on f_c
        # Formulas are based on EHE-08 specifications
        self._f_ck = self.f_c
        self._f_cm = self._f_ck+8*MPa #type:ignore
        self._E_cm = 8500 *(self._f_cm / MPa)**(1/3) * MPa  # type: ignore
        self._f_ctm = 0.3 * (self._f_ck / MPa)**(2/3) * MPa  # type: ignore

    def get_properties(self) -> dict:
        properties = super().get_properties()
        properties['E_cm'] = self._E_cm
        properties['f_ck'] = self._f_ck
        properties['f_cm'] = self._f_cm
        properties['f_ctm'] = self._f_ctm
        return properties

    def alpha_cc(self):
        # Example implementation for alpha_cc as per EHE-08
        return 0.85  # Example value, modify according to EHE-08 standards

# Factory function
def create_concrete(name: str, f_c: float, design_code: str):
    if design_code == "ACI 318-19":
        return Concrete_ACI_318_19(name=name, f_c=f_c)
    elif design_code == "EN 1992":
        return Concrete_EN_1992(name=name, f_c=f_c)
    elif design_code == "EHE-08":
        return Concrete_EHE_08(name=name, f_c=f_c)
    else:
        raise ValueError(f"Invalid design code: {design_code}")

@dataclass
class Steel(Material):
    _f_y: float = field(init=False)
    _density: float = field(default=7850 * kg / m**3)  # type: ignore

    def __init__(self, name: str, f_y: float, density: float = 7850 * kg / m**3):
        super().__init__(name)
        self._f_y = f_y
        self._density = density

@dataclass
class SteelBar(Steel):
    _E_s: float = field(default=200000*MPa) # type: ignore
    _epsilon_y: float = field(init=False)

    def __init__(self, name: str, f_y: float, density: float =7850 *kg/m**3):
        super().__init__(name, f_y, density)
        self._epsilon_y = self._f_y / self._E_s # 21.2.2.1 - Page 392

    def get_properties(self) -> dict:
        properties = super().get_properties()
        properties['E_s'] = self._E_s
        properties['f_y'] = self._f_y
        properties['epsilon_ty']=self._epsilon_y
        return properties

@dataclass
class SteelStrand(Steel):
    _f_u: float = field(default=1860*MPa) # type: ignore
    _E_s: float = field(default=190000*MPa) # type: ignore
    prestress_stress: float = field(default=0)

    def __init__(self, name: str, f_y: float, density: float =7850 *kg/m**3):
        super().__init__(name, f_y, density)
        self._epsilon_y = self._f_y / self._E_s


    def get_properties(self) -> dict:
        properties = super().get_properties()
        properties['E_s'] = self._E_s
        properties['f_y'] = self._f_y
        properties['f_u'] = self._f_u
        return properties

def main():
    # Test cases
    concrete=create_concrete(name="H25",f_c=25*MPa, design_code="ACI 318-19") # type: ignore
    print(concrete.get_properties())
    print(concrete.E_c,', ',concrete.f_r)
    steelbar = SteelBar(name="ADN 500",f_y=500*MPa) # type: ignore
    print(steelbar.get_properties())
    steelstrand = SteelStrand(name='Y1860',f_y=1700*MPa) # type: ignore
    print(steelstrand.get_properties())


if __name__ == "__main__":
    main()


