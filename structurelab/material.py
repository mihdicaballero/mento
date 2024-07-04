# common/materials.py

from __future__ import annotations
from dataclasses import dataclass, field
import math
import forallpeople
forallpeople.environment('structural',top_level=True)
#Useful extra units
cm = 1e-2*m

@dataclass
class Material:
    def __init__(self, name: str):
        self.name = name

    def get_properties(self):
        return {'name': self.name}
@dataclass
class Concrete(Material):
    def __init__(self, name: str, f_c: float):
        super().__init__(name)
        self.f_c = f_c
        self.density = 2500 * kg / m**3 

    def get_properties(self):
        properties = super().get_properties()
        properties['f_c'] = self.f_c
        properties['density'] = self.density
        return properties
@dataclass
class ConcreteACI(Concrete):
    def __init__(self, name: str, f_c: float=25*MPa, design_code: str="318-19"):
        super().__init__(name, f_c)
        self.design_code = design_code
        self.E_c = self.calculate_E_c()
        self.f_r = self.calculate_f_r()

    def get_properties(self):
        properties = super().get_properties()
        properties['design_code'] = self.design_code
        properties['E_c'] = self.E_c
        properties['f_r'] = self.f_r
        return properties
    def calculate_E_c(self):
        return ((self.density / (kg / m**3)) ** 1.5) * 0.043 * math.sqrt(self.f_c / MPa) * MPa
    def calculate_f_r(self):
        return 0.625*math.sqrt(self.f_c/MPa)*MPa
    def beta_1(self):
        if 17 <= self.f_c / MPa <= 28:
            return 0.85
        elif 28 < self.f_c / MPa <= 56:
            return 0.85 - 0.05/7 * (self.f_c / MPa - 20)
        elif 56 < self._c / MPa:
            return 0.65


@dataclass
class Steel(Material):
    def __init__(self, name: str, f_y: float):
        super().__init__(name)
        self.f_y = f_y
        self.density: float = 7850*kg/m**3

@dataclass
class SteelBar(Steel):
    def __init__(self, name: str, f_y: float=420*MPa):
        super().__init__(name, f_y)
        self.E_s = 200*MPa
    def get_properties(self):
        properties = super().get_properties()
        properties['E_s'] = self.E_s
        properties['f_y'] = self.f_y
        return properties

@dataclass
class SteelStrand(Steel):
    def __init__(self, name: str, f_y: float=1700*MPa):
        super().__init__(name, f_y)
        self.f_u = 1860*MPa
        self.E_s = 190*MPa
        self_prestress_stress = 0
    def get_properties(self):
        properties = super().get_properties()
        properties['E_s'] = self.E_s
        properties['f_y'] = self.f_y
        properties['f_u'] = self.f_u
        return properties



concrete=ConcreteACI(name="H25",f_c=30*MPa, design_code="318-19")
print(concrete.get_properties())
print(concrete.beta_1())
print(concrete.E_c,', ',concrete.f_r)
steelbar = SteelBar(name="ADN 500",f_y=500*MPa)
print(steelbar.get_properties())
steelstrand = SteelStrand(name='Y1860',f_y=1700*MPa)
print(steelstrand.get_properties())
