from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, TYPE_CHECKING
import math
from mento.units import kg, m, MPa, ksi, GPa, psi, Pa, lb, ft
from mento import ureg

# Conditional import for type checking only
if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity


@dataclass
class Material:
    name: str  # Name will be provided by subclasses or instances


@dataclass
class Concrete(Material):
    f_c: PlainQuantity = field(default=25 * MPa)
    design_code: str = field(default="ACI 318-19")
    unit_system: str = field(init=False)

    def __post_init__(self) -> None:
        # Detect the unit system based on f_c
        if self.f_c.units == MPa or self.f_c.units == Pa:
            self.unit_system = "metric"
            self.density: PlainQuantity = 2500 * kg / m**3
        elif self.f_c.units == psi or self.f_c.units == ksi:
            self.unit_system = "imperial"
            self.density = 155 * lb / ft**3
        else:
            raise ValueError("Unsupported unit system for f_c. Please use MPa or ksi.")

    def get_properties(self) -> Dict[str, PlainQuantity]:
        # Return properties in the appropriate unit system
        properties = {"f_c": self.f_c, "density": self.density}
        return properties


@dataclass
class Concrete_ACI_318_19(Concrete):
    _E_c: PlainQuantity = field(init=False)
    _f_r: PlainQuantity = field(init=False)
    _epsilon_c: float = field(default=0.003, init=False)
    _beta_1: float = field(init=False)

    def __post_init__(self) -> None:
        super().__post_init__()
        # Ensure name is properly set, either hardcode or pass during instantiation
        if not self.name:
            self.name = "Concrete ACI 318-19"  # You can set a default name here
        self.design_code = "ACI 318-19"
        # Adjust calculations based on unit system
        if self.unit_system == "metric":
            self._E_c = (
                ((self.density / (kg / m**3)) ** 1.5)
                * 0.043
                * math.sqrt(self.f_c / MPa)
                * MPa
            )
            self._f_r = 0.625 * math.sqrt(self.f_c / MPa) * MPa
        else:  # imperial
            self._E_c = (
                ((self.density / (lb / ft**3)) ** 1.5)
                * 33
                * math.sqrt(self.f_c / psi)
                * psi
            )
            self._f_r = 7.5 * math.sqrt(self.f_c / psi) * psi
        self._beta_1 = self.__beta_1()

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = super().get_properties()
        properties["E_c"] = self._E_c
        properties["f_r"] = self._f_r
        properties["beta_1"] = ureg.Quantity(self._beta_1, "")
        properties["epsilon_c"] = ureg.Quantity(self._epsilon_c, "")
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

    def __str__(self) -> str:
        """Customize the string representation for user-friendly display."""
        properties = self.get_properties()
        return (
            f"Concrete Properties ({self.name}):\n"
            f"  f_c: {properties['f_c']}\n"
            f"  Density: {properties['density']}\n"
            f"  E_c: {properties['E_c']}\n"
            f"  f_r: {properties['f_r']}\n"
            f"  beta_1: {properties['beta_1']}\n"
            f"  epsilon_c: {properties['epsilon_c'].magnitude:.3f}"
        )


@dataclass
class Concrete_CIRSOC_201_25(Concrete_ACI_318_19):
    """Concrete class for a CIRSOC 201 design code with metric units."""

    def __post_init__(self) -> None:
        # Call the parent class's __post_init__ to inherit initializations
        super().__post_init__()
        # Override the design code for this specific class
        self.design_code = "CIRSOC 201-25"
        # Ensure the name is specific to this design code
        if not self.name:
            self.name = "Concrete CIRSOC 201-25"


@dataclass
class Concrete_EN_1992_2004(Concrete):
    _E_cm: PlainQuantity = field(init=False)  # Secant modulus of elasticity
    _f_ck: PlainQuantity = field(init=False)  # Characteristic concrete strength
    _f_cm: PlainQuantity = field(init=False)  # Mean compressive strength
    _f_ctm: PlainQuantity = field(init=False)  # Mean tensile strength
    _epsilon_cu3: float = field(init=False)

    def __init__(self, name: str, f_ck: PlainQuantity):
        super().__init__(name=name, f_c=f_ck)
        if not self.name:
            self.name = "Concrete EN 1992-2004"
        self.design_code = "EN 1992-2004"
        self._delta=0.85
        self._f_ck = self.f_c
        self._f_cm = self._f_ck + 8 * MPa
        self._E_cm = 22000 * (self._f_cm / (10 * MPa)) ** 0.3 * MPa
        self._f_ctm = 0.3 * (self._f_ck / MPa) ** (2 / 3) * MPa
        #TODO ARREGLAR ESTO, HAY QUE CALCULAR EL epsilon_cu3 para hormigones de mas de 50MPa:
        if (self._f_ck>50*MPa):
            print("WARNING: epsilon_cu3 is set to 0.0035, which is valid for concretes with f_ck ≤ 50 MPa.")
        self._epsilon_cu3 = 0.0035  # Ultimate strain in concrete

    def alpha_cc(self) -> float:
        # Implementation for alpha_cc, as per Eurocode EN 1992-1-1
        # Designers Guide to EN 1992-1-1, Page 62
        # Study of the data available on the behaviour of compression zones at failure suggests that the
        # use of 1.0 is unconservative. For this reason, the UK National Annex recommends a value for
        # αcc of 0.85, as is proposed in the CEB Model Codes.
        return 0.85

    def lambda_factor(self) -> float:
        """
        Calculate the effective compression zone depth factor (λ) as per EN 1992-1-1.
        """
        if self._f_ck <= 50 * MPa:
            return 0.8
        else:
            return 0.8 - (self._f_ck / MPa - 50) / 400

    def eta_factor(self) -> float:
        """
        Calculate the effective compressive strength factor (η) as per EN 1992-1-1.
        """
        if self._f_ck <= 50 * MPa:
            return 1.0
        else:
            return 1.0 - (self._f_ck / MPa - 50) / 200

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = super().get_properties()
        properties.update(
            {
                "E_cm": self._E_cm.to("MPa"),
                "f_ck": self._f_ck.to("MPa"),
                "f_cm": self._f_cm.to("MPa"),
                "f_ctm": self._f_ctm.to("MPa"),
            }
        )
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


@dataclass
class Steel(Material):
    _f_y: PlainQuantity = field(init=False)
    _density: PlainQuantity = field(default=7850 * kg / m**3)

    def __init__(
        self, name: str, f_y: PlainQuantity, density: PlainQuantity = 7850 * kg / m**3
    ):
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
    _E_s: PlainQuantity = field(default=200 * GPa)
    _epsilon_y: PlainQuantity = field(init=False)

    def __init__(
        self, name: str, f_y: PlainQuantity, density: PlainQuantity = 7850 * kg / m**3
    ):
        super().__init__(name, f_y, density)
        self._epsilon_y = f_y.to("MPa") / (self._E_s.to("MPa"))  # 21.2.2.1 - Page 392

    @property
    def E_s(self) -> PlainQuantity:
        return self._E_s

    @property
    def epsilon_y(self) -> PlainQuantity:
        return self._epsilon_y

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = {
            "E_s": self._E_s.to("GPa"),
            "f_y": self._f_y.to("MPa"),
            "epsilon_y": self._epsilon_y,
        }
        return properties

    def __str__(self) -> str:
        """Customize the string representation for user-friendly display."""
        properties = self.get_properties()
        return (
            f"SteelBar Properties ({self.name}):\n"
            f"  f_y: {properties['f_y']}\n"
            f"  E_s: {properties['E_s']}\n"
            f"  epsilon_y: {properties['epsilon_y'].magnitude:.4f}\n"
            f"  Density: {self.density}"
        )


@dataclass
class SteelStrand(Steel):
    _f_u: PlainQuantity = field(default=1860 * MPa)
    _E_s: PlainQuantity = field(default=190000 * MPa)
    prestress_stress: PlainQuantity = field(default=0 * MPa)

    def __init__(
        self, name: str, f_y: PlainQuantity, density: PlainQuantity = 7850 * kg / m**3
    ):
        super().__init__(name, f_y, density)
        self._epsilon_y = self._f_y / self._E_s

    def get_properties(self) -> Dict[str, PlainQuantity]:
        properties = {
            "E_s": self._E_s.to("MPa"),
            "f_y": self._f_y.to("MPa"),
            "f_u": self._f_u.to("MPa"),
        }
        return properties

    @property
    def f_u(self) -> PlainQuantity:
        return self._f_u

    @property
    def E_s(self) -> PlainQuantity:
        return self._E_s
