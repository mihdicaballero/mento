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
    name: str
@dataclass
class Concrete(Material):
    design_code: str ="ACI" # Attribute to specify the design code (ACI or Eurocode)
    f_c: float = 25*MPa
    density: float = 2500 * kg / m**3 

    def __post_init__(self):
        # Validate the design code
        valid_design_codes = ['ACI', 'Eurocode']
        if self.design_code not in valid_design_codes:
            raise ValueError(f"Invalid design code. Allowed values are: {', '.join(valid_design_codes)}")

    def get_concrete_properties(self):
        # Concrete properties specific to the design code
        strategy = ConcretePropertiesStrategy.get_strategy(self.design_code)
        return strategy.get_concrete_properties(self)
   
@dataclass
class ConcretePropertiesStrategy:
    @staticmethod
    def get_strategy(design_code: str) -> ConcretePropertiesStrategy:
        if design_code == 'ACI':
            return ACIConcretePropertiesStrategy()
        elif design_code == 'Eurocode':
            return EurocodeConcretePropertiesStrategy()
        else:
            raise ValueError(f"Invalid design code: {design_code}")

    def get_concrete_properties(self, concrete: Concrete) -> dict:
        raise NotImplementedError("Subclasses must implement this method")

@dataclass
class ACIConcretePropertiesStrategy(ConcretePropertiesStrategy):
    def get_concrete_properties(self, concrete: Concrete) -> dict:
        # Implement ACI-specific concrete properties calculations
        properties = {
            'E_c': ((concrete.density / (kg / m**3)) ** 1.5) * 0.043 * math.sqrt(concrete.f_c / MPa) * MPa,
            'f_r':0.625*math.sqrt(concrete.f_c/MPa)*MPa,

            # Add more properties as needed
        }
        return properties
    @staticmethod
    def beta_1(f_c):
        if 17 <= f_c / MPa <= 28:
            return 0.85
        elif 28 < f_c / MPa <= 56:
            return 0.85 - 0.05/7 * (f_c / MPa - 20)
        elif 56 < f_c / MPa:
            return 0.65

@dataclass
class EurocodeConcretePropertiesStrategy(ConcretePropertiesStrategy):
    def get_concrete_properties(self, concrete: Concrete) -> dict:
        # Implement Eurocode-specific concrete properties calculations
        properties = {
            'E_c':  31*GPa,
            # Add more properties as needed
        }
        return properties
        
@dataclass
class Steel(Material):
    density: float = 7850*kg/m**3

@dataclass
class SteelBar(Steel):
    # Default values for density and modulus of elasticity for steel
    f_y: float = 420*MPa
    E_s: float = 200*GPa

@dataclass
class SteelStrand(Steel):
    yield_strength: float = 1700*MPa
    tensile_strenght: float = 1860*MPa
    E_s: float = 190*GPa
    prestress_stress: float = 0

    def get_prestress_stress(self) -> float:
        return self.prestress_stress

    def get_prestress_strain(self) -> float:
        # You can define the strain calculation based on your specific needs
        pass

