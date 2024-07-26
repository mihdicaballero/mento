# common/materials.py

from __future__ import annotations
from dataclasses import dataclass, field
import math
import forallpeople
forallpeople.environment('structural',top_level=True)
#Useful extra units
cm = 1e-2*m # type: ignore

@dataclass
class Material:
    def __init__(self, name: str):
        self.name = name

    def get_properties(self):
        return {'name': self.name}
@dataclass
class Concrete(Material):
    def __init__(self, name: str, f_c=25*MPa, design_code: str="ACI 318-19"): # type: ignore
        super().__init__(name)
        self.f_c = f_c
        self.density = 2500 * kg / m**3 # type: ignore
        self.design_code = design_code

    def get_properties(self):
        properties = super().get_properties()
        properties['f_c'] = self.f_c
        properties['density'] = self.density
        properties['design_code'] = self.design_code
        return properties
@dataclass
class ConcreteACI31819(Concrete):
    def __init__(self, name: str, f_c: float, design_code):
        super().__init__(name, f_c, design_code)
        self.E_c =  ((self.density / (kg / m**3)) ** 1.5) * 0.043 * math.sqrt(self.f_c / MPa) * MPa # type: ignore
        self.f_r = 0.625*math.sqrt(self.f_c/MPa)*MPa # type: ignore
        self.epsilon_c = 0.003

    def get_properties(self):
        properties = super().get_properties()       
        properties['E_c'] = self.E_c
        properties['f_r'] = self.f_r
        return properties
    
    def beta_1(self):
        if 17 <= self.f_c / MPa <= 28: # type: ignore
            return 0.85
        elif 28 < self.f_c / MPa <= 56: # type: ignore
            return 0.85 - 0.05/7 * (self.f_c / MPa - 20) # type: ignore
        elif 56 < self._c / MPa: # type: ignore
            return 0.65
@dataclass
class ConcreteEN1992(Concrete):
    # Example properties specific to EN 1992
    def __init__(self, name: str, f_c: float, design_code: str = "EN 1992"):
        super().__init__(name, f_c, design_code)
        # Additional EN 1992 specific initializations

# Factory function
def create_concrete(name: str, f_c: float, design_code: str):
    if design_code == "ACI 318-19":
        return ConcreteACI31819(name, f_c, design_code)
    elif design_code == "EN 1992":
        return ConcreteEN1992(name, f_c, design_code)
    else:
        raise ValueError(f"Unsupported design code: {design_code}")

@dataclass
class Steel(Material):
    def __init__(self, name: str, f_y: float):
        super().__init__(name)
        self.f_y = f_y
        self.density: float = 7850*kg/m**3 # type: ignore

@dataclass
class SteelBar(Steel):
    def __init__(self, name: str, f_y: float=420*MPa): # type: ignore
        super().__init__(name, f_y)
        self.E_s = 200*MPa # type: ignore
        self.epsilon_y = self.f_y/self.E_s
    def get_properties(self):
        properties = super().get_properties()
        properties['E_s'] = self.E_s
        properties['f_y'] = self.f_y
        return properties

@dataclass
class SteelStrand(Steel):
    def __init__(self, name: str, f_y: float=1700*MPa): # type: ignore
        super().__init__(name, f_y)
        self.f_u = 1860*MPa # type: ignore
        self.E_s = 190*MPa # type: ignore
        self_prestress_stress = 0 
    def get_properties(self):
        properties = super().get_properties()
        properties['E_s'] = self.E_s
        properties['f_y'] = self.f_y
        properties['f_u'] = self.f_u
        return properties


concrete=create_concrete(name="H30",f_c=30*MPa, design_code="ACI 318-19") # type: ignore
print(concrete.get_properties())
print(concrete.beta_1())
print(concrete.E_c,', ',concrete.f_r)
steelbar = SteelBar(name="ADN 500",f_y=500*MPa) # type: ignore
print(steelbar.get_properties())
steelstrand = SteelStrand(name='Y1860',f_y=1700*MPa) # type: ignore
print(steelstrand.get_properties())