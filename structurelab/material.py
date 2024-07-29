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
    name: str

    def get_properties(self):
        return {'name': self.name}
@dataclass
class Concrete(Material):
    f_c: float = field(default=25*MPa) # type: ignore
    design_code: str = field(default="ACI 318-19")
    density: float = 2500*kg/m**3 # type: ignore

    def get_properties(self) -> dict:
        properties = super().get_properties()
        properties['f_c'] = self.f_c
        properties['design_code'] = self.design_code
        properties['density'] = self.density
        return properties
@dataclass
class Concrete_ACI_318_19(Concrete):
    E_c: float = field(init=False)
    f_r: float = field(init=False)
    epsilon_c: float = field(default=0.003, init=False)

    def __post_init__(self):
        self.E_c =  ((self.density / (kg / m**3)) ** 1.5) * 0.043 * math.sqrt(self.f_c / MPa) * MPa # type: ignore
        self.f_r = 0.625*math.sqrt(self.f_c/MPa)*MPa # type: ignore

    def get_properties(self) -> dict:
        properties = super().get_properties()       
        properties['E_c'] = self.E_c
        properties['f_r'] = self.f_r
        return properties
    
    def beta_1(self):
        if 17 <= self.f_c / MPa <= 28: # type: ignore
            return 0.85
        elif 28 < self.f_c / MPa <= 56: # type: ignore
            return 0.85 - 0.05/7 * (self.f_c / MPa - 20) # type: ignore
        elif 56 < self.f_c / MPa: # type: ignore
            return 0.65
@dataclass
class Concrete_EN_1992(Concrete):
    E_cm: float = field(init=False)  # Secant modulus of elasticity
    f_ck: float = field(init=False) # Characteristic concrete strength
    f_cm: float = field(init=False) # mean compressive strength
    f_ctm: float = field(init=False) # Mean tensile strength

    def __post_init__(self):
        # Calculate E_cm and f_ctm based on f_c and density
        self.f_ck = self.f_c
        self.f_cm = self.f_ck+8*MPa #type:ignore
        self.E_cm = 22000 * (self.f_cm / (10*MPa)) ** (0.3) * MPa # type: ignore       
        self.f_ctm = 0.3 * (self.f_ck / MPa) ** (2/3) * MPa # type: ignore

    def get_properties(self) -> dict:
        properties = super().get_properties()
        properties['E_cm'] = self.E_cm
        properties['f_ctm'] = self.f_ctm
        return properties

    def alpha_cc(self):
        # Example implementation for alpha_cc, as per Eurocode EN 1992-1-1
        return 1.0  # Typically, this value is taken as 1.0 for normal weight concrete

@dataclass
class Concrete_EHE_08(Concrete):
    E_cm: float = field(init=False)  # Secant modulus of elasticity
    f_ck: float = field(init=False) # Characteristic concrete strength
    f_cm: float = field(init=False) # mean compressive strength
    f_ctm: float = field(init=False) # Mean tensile strength
    f_ctm_fl: float = field(init=False) # Mean flexure tensile strength

    def __post_init__(self):
        # Calculate E_cm and f_ctk based on f_c
        # Formulas are based on EHE-08 specifications
        self.f_ck = self.f_c
        self.f_cm = self.f_ck+8*MPa #type:ignore
        self.E_cm = 8500 *(self.f_cm / MPa)**(1/3) * MPa  # type: ignore
        self.f_ctm = 0.3 * (self.f_ck / MPa)**(2/3) * MPa  # type: ignore

    def get_properties(self) -> dict:
        properties = super().get_properties()
        properties['E_cm'] = self.E_cm
        properties['f_ck'] = self.f_ck
        properties['f_cm'] = self.f_cm
        properties['f_ctm'] = self.f_ctm
        return properties

    def alpha_cc(self):
        # Example implementation for alpha_cc as per EHE-08
        return 0.85  # Example value, modify according to EHE-08 standards

# Factory function
def create_concrete(name: str, f_c: float, design_code: str):
    concrete_classes = {
        "ACI 318-19": Concrete_ACI_318_19,
        "EN 1992": Concrete_EN_1992,
        "EHE-08": Concrete_EHE_08,
    }
    if design_code not in concrete_classes:
        raise ValueError(f"Invalid design code: {design_code}")
    return concrete_classes[design_code](name, f_c, design_code)

@dataclass
class Steel(Material):
    f_y: float
    density: float = field(default=7850*kg/m**3) # type: ignore

@dataclass
class SteelBar(Steel):
    E_s: float = field(default=200*MPa) # type: ignore
    epsilon_y: float = field(init=False)

    def __post_init__(self):
        self.epsilon_y = self.f_y / self.E_s
    def get_properties(self) -> dict:
        properties = super().get_properties()
        properties['E_s'] = self.E_s
        properties['f_y'] = self.f_y
        return properties

@dataclass
class SteelStrand(Steel):
    f_u: float = field(default=1860*MPa) # type: ignore
    E_s: float = field(default=190000*MPa) # type: ignore
    prestress_stress: float = field(default=0)

    def get_properties(self) -> dict:
        properties = super().get_properties()
        properties['E_s'] = self.E_s
        properties['f_y'] = self.f_y
        properties['f_u'] = self.f_u
        return properties

def main():
    # Test cases
    concrete=create_concrete(name="H25",f_c=25*MPa, design_code="EN 1992") # type: ignore
    print(concrete.get_properties())
    print(concrete.E_c,', ',concrete.f_r)
    steelbar = SteelBar(name="ADN 500",f_y=500*MPa) # type: ignore
    print(steelbar.get_properties())
    steelstrand = SteelStrand(name='Y1860',f_y=1700*MPa) # type: ignore
    print(steelstrand.get_properties())

if __name__ == "__main__":
    main()



