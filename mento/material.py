from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, TYPE_CHECKING, Any
import math
from mento.units import kg, m, MPa, ksi, GPa, psi, Pa, lb, ft, kPa

# Conditional import for type checking only
if TYPE_CHECKING:
    from pint import Quantity


@dataclass
class Material:
    name: str  # Name will be provided by subclasses or instances


@dataclass
class Concrete(Material):
    f_c: Quantity = field(default=25 * MPa)
    design_code: str = field(default="ACI 318-19", init=False)
    unit_system: str = field(init=False)
    density: Quantity = field(init=False)

    def __post_init__(self) -> None:
        # Detect the unit system based on f_c
        if self.f_c.units == MPa or self.f_c.units == Pa or self.f_c.units == kPa:
            self.unit_system = "metric"
            self.density: Quantity = 2500 * kg / m**3
        elif self.f_c.units == psi or self.f_c.units == ksi:
            self.unit_system = "imperial"
            self.density = 155 * lb / ft**3
        else:
            raise ValueError(
                f"Unsupported unit system for f_c ({self.f_c.units}). Please use MPa, Pa, kPa, psi, or ksi."
            )

    def get_properties(self) -> Dict[str, Any]:
        # Return properties in the appropriate unit system
        properties = {"f_c": self.f_c, "density": self.density}
        return properties


@dataclass
class Concrete_ACI_318_19(Concrete):
    """
    Concrete_ACI_318_19 represents a concrete material model based on the ACI 318-19 code provisions.
    This class extends the base `Concrete` class and implements calculation of key material properties
    according to the American Concrete Institute (ACI) 318-19 standard, including modulus of elasticity,
    modulus of rupture, and stress block parameters. It supports both metric and imperial unit systems.

    Inputs:
        name: Name of the concrete material.
        f_c: Characteristic compressive strength of concrete.

    Methods:
        get_properties() -> Dict[str, Any]: Returns a dictionary of concrete properties.
        E_c: Returns modulus of elasticity.
        f_r (property): Returns modulus of rupture.
        beta_1 (property): Returns β₁ value.
        lambda_factor (property): Returns λ factor for lightweight concrete.
        phi_v (property): Returns shear strength reduction factor.
        phi_c (property): Returns compression-controlled strength reduction factor.
        phi_y (property): Returns tension-controlled strength reduction factor.

    Usage:
        Instantiate this class to represent a concrete material with properties and code factors
        compliant with ACI 318-19, suitable for use in structural analysis and design calculations.
    """

    _E_c: Quantity = field(init=False)
    _f_r: Quantity = field(init=False)
    _epsilon_c: float = field(default=0.003, init=False)
    _beta_1: float = field(init=False)
    _lambda: float = field(init=False)
    _phi_v: float = field(init=False)
    _phi_c: float = field(init=False)
    _phi_t: float = field(init=False)
    _flexural_min_reduction: bool = field(init=False)

    def __post_init__(self) -> None:
        super().__post_init__()
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
        self._lambda = 1  # Normalweight concrete
        self._phi_v = 0.75  # Shear strength reduction factor
        self._phi_c = 0.65  # Compression controlled strength reduction factor
        self._phi_t = 0.90  # Tension controlled strength reduction factor
        self._flexural_min_reduction = (
            True  # True selects 4/3 of calculated steel if it's less than minimum
        )

    def get_properties(self) -> Dict[str, Any]:
        properties = super().get_properties()
        properties["E_c"] = self._E_c
        properties["f_r"] = self._f_r
        properties["beta_1"] = self._beta_1
        properties["epsilon_c"] = self._epsilon_c
        properties["lambda"] = self._lambda
        properties["phi_v"] = self._phi_v
        properties["phi_c"] = self._phi_c
        properties["phi_t"] = self._phi_t
        return properties

    def __beta_1(self) -> float:
        # Table 22.2.2.4.3—Values of β1 for equivalent rectangular concrete stress distribution
        # Page 399
        fc_MPa = self.f_c.to("MPa").magnitude  # Ensure comparison in MPa
        if 17 <= fc_MPa <= 28:
            return 0.85
        elif 28 < fc_MPa <= 55:
            return 0.85 - 0.05 / 7 * (fc_MPa - 28)
        elif fc_MPa > 55:
            return 0.65
        else:
            # Handle case where f_c / MPa < 17
            return 0.85

    @property
    def E_c(self) -> Quantity:
        return self._E_c

    @property
    def f_r(self) -> Quantity:
        return self._f_r

    @property
    def beta_1(self) -> float:
        return self._beta_1

    @property
    def lambda_factor(self) -> float:
        return self._lambda

    @property
    def phi_v(self) -> float:
        return self._phi_v

    @property
    def phi_c(self) -> float:
        return self._phi_c

    @property
    def phi_y(self) -> float:
        return self._phi_t

    def __str__(self) -> str:
        properties = self.get_properties()
        return (
            f"Concrete Properties ({self.name}):\n"
            f"  f_c: {properties['f_c']}\n"
            f"  Density: {properties['density']}\n"
            f"  E_c: {properties['E_c']}\n"
            f"  f_r: {properties['f_r']}\n"
            f"  beta_1: {properties['beta_1']:.3f}\n"  # Format as float
            f"  epsilon_c: {properties['epsilon_c']:.4f}\n"  # Format as float
            f"  λ: {properties['lambda']:.2f}\n"
            f"  phi_v: {properties['phi_v']:.2f}\n"
            f"  phi_c: {properties['phi_c']:.2f}\n"
            f"  phi_t: {properties['phi_t']:.2f}\n"
        )


@dataclass
class Concrete_CIRSOC_201_25(Concrete_ACI_318_19):
    """Concrete class for a CIRSOC 201-25 design code with metric units."""

    def __post_init__(self) -> None:
        # Call the parent class's __post_init__ to inherit initializations
        super().__post_init__()
        # Override the design code for this specific class
        self.design_code = "CIRSOC 201-25"


@dataclass
class Concrete_EN_1992_2004(Concrete):
    """
    Concrete_EN_1992_2004 represents concrete material properties and design parameters according to Eurocode EN 1992-1-1:2004.
    This class extends the base `Concrete` class, providing Eurocode-specific calculations for characteristic and mean strengths,
    modulus of elasticity, and other design factors. It encapsulates the following key methods:

    Inputs:
        name: Name of the concrete material.
        f_c: Characteristic compressive strength of concrete (f_ck for Eurocode).

    Methods:
        get_properties() -> Dict[str, Any]: Returns a dictionary of all relevant material properties.
        E_cm (property): Returns the secant modulus of elasticity.
        f_ck (property): Returns the characteristic compressive strength.
        f_cm (property): Returns the mean compressive strength.
        f_ctm (property): Returns the mean tensile strength.
        epsilon_cu3 (property): Returns the ultimate strain in concrete.
        gamma_c (property): Returns the partial safety factor for concrete.
        alpha_cc (property): Returns the α_cc coefficient.
        Lambda_factor (property): Returns the λ factor.
        Eta_factor (property): Returns the η factor.

    Usage:
        This class is intended for use in structural engineering applications where concrete properties must comply with EN 1992-1-1:2004.
        It provides all necessary parameters for design and verification according to the code.
    """

    def __post_init__(self) -> None:
        # Crucial: Call parent's __post_init__ first to set unit_system and density
        super().__post_init__()
        self.design_code = "EN 1992-2004"

        # The f_c passed to Concrete is the f_ck for Eurocode
        self._delta=0.85
        self._f_ck = self.f_c
        self._f_cm = self._f_ck + 8 * MPa
        self._E_cm = 22000 * (self._f_cm.to("MPa").magnitude / 10) ** 0.3 * MPa
        self._f_ctm = 0.3 * (self._f_ck.to("MPa").magnitude) ** (2 / 3) * MPa
        self._epsilon_cu1 = (2.8 + 27 * ((98 - self._f_cm.to("MPa").magnitude) / 100) ** 4 if self._f_ck >= 50 * MPa else 3.5) * 1e-3
        self._epsilon_c2  = (2.0 + 0.085 * ((self._f_ck.to("MPa").magnitude - 50)) ** 0.53 if self._f_ck >= 50 * MPa else 2.0) * 1e-3
        self._epsilon_cu2 = (2.6 + 35 * ((90 - self._f_ck.to("MPa").magnitude) / 100) ** 4 if self._f_ck >= 50 * MPa else 3.5) * 1e-3
        self._epsilon_c3  = (1.75 + 0.55 * ((self._f_ck.to("MPa").magnitude - 50) / 40) if self._f_ck >= 50 * MPa else 1.75) * 1e-3
        self._epsilon_cu3 = (2.6 + 35 * ((90 - self._f_ck.to("MPa").magnitude) / 100) ** 4 if self._f_ck >= 50 * MPa else 3.5) * 1e-3
        self._gamma_c = 1.5
        self._gamma_s = 1.15
        self._alpha_cc = self._alpha_cc_calc()
        # Default values for k values EN_1992-1-1 - ART 5.5:
        self._k_1=0.44
        self._k_2=1.25*(0.6+0.0014/self._epsilon_cu2)
        self._k_3=0.54
        self._k_4=1.25*(0.6+0.0014/self._epsilon_cu2)
        self._k_5=0.7
        self._k_6=0.8*self._epsilon_cu2

    def _alpha_cc_calc(self) -> float:
        # Implementation for alpha_cc, as per Eurocode EN 1992-1-1
        # Designers Guide to EN 1992-1-1, Page 62
        # Study of the data available on the behaviour of compression zones at failure suggests that the
        # use of 1.0 is unconservative. For this reason, the UK National Annex recommends a value for
        # αcc of 0.85, as is proposed in the CEB Model Codes.
        return 0.85
    
    def _lambda_factor(self) -> float:
        """
        Calculate the effective compression zone depth factor (λ) as per EN 1992-1-1.
        """
        f_ck_mpa = self._f_ck.to("MPa").magnitude  # Ensure comparison in MPa
        if f_ck_mpa <= 50:
            return 0.8
        else:
            return 0.8 - (f_ck_mpa - 50) / 400

    def _eta_factor(self) -> float:
        """
        Calculate the effective compressive strength factor (η) as per EN 1992-1-1.
        """
        f_ck_mpa = self._f_ck.to("MPa").magnitude  # Ensure comparison in MPa
        if f_ck_mpa <= 50:
            return 1.0
        else:
            return 1.0 - (f_ck_mpa - 50) / 200

    def get_properties(self) -> Dict[str, Any]:
        properties = super().get_properties()
        properties.update(
            {
                "E_cm": self._E_cm,
                "f_ck": self._f_ck,
                "f_cm": self._f_cm,
                "f_ctm": self._f_ctm,
                "epsilon_cu3": self._epsilon_cu3,
                "gamma_c": self._gamma_c,
                "gamma_s": self._gamma_s,
                "alpha_cc": self._alpha_cc,
                "lambda_factor": self._lambda_factor(),
                "eta_factor": self._eta_factor(),
            }
        )
        return properties

    @property
    def E_cm(self) -> Quantity:
        return self._E_cm

    @property
    def f_ck(self) -> Quantity:
        return self._f_ck

    @property
    def f_cm(self) -> Quantity:
        return self._f_cm

    @property
    def f_ctm(self) -> Quantity:
        return self._f_ctm

    @property
    def epsilon_cu3(self) -> float:
        return self._epsilon_cu3

    @property
    def gamma_c(self) -> float:
        return self._gamma_c
    
    @property
    def gamma_s(self) -> float:
        return self._gamma_s

    @property
    def alpha_cc(self) -> float:  # This property returns the calculated alpha_cc
        return self._alpha_cc

    @property
    def Lambda_factor(self) -> float:  # Property for lambda_factor
        return self._lambda_factor()

    @property
    def Eta_factor(self) -> float:  # Property for eta_factor
        return self._eta_factor()

    def __str__(self) -> str:
        """Customize the string representation for user-friendly display."""
        properties = self.get_properties()
        # Access magnitude for dimensionless quantities
        return (
            f"Concrete Properties ({self.name}):\n"
            f"  Design Code: {self.design_code}\n"
            f"  f_c (Characteristic): {properties['f_ck']}\n"
            f"  f_cm (Mean Compressive): {properties['f_cm']}\n"
            f"  f_ctm (Mean Tensile): {properties['f_ctm']}\n"
            f"  E_cm (Secant Modulus): {properties['E_cm']}\n"
            f"  Density: {properties['density']}\n"
            f"  ε_cu3: {properties['epsilon_cu3']:.4f}\n"
            f"  γ_c: {properties['gamma_c']:.2f}\n"
            f"  γ_s: {properties['gamma_s']:.2f}\n"
            f"  α_cc: {properties['alpha_cc']:.2f}\n"
            f"  λ Factor: {properties['lambda_factor']:.2f}\n"
            f"  Eta Factor: {properties['eta_factor']:.2f}"
        )


@dataclass
class Steel(Material):
    _f_y: Quantity = field(init=False)
    _density: Quantity = field(default=7850 * kg / m**3)
    
    def __init__(
        self, name: str, f_y: Quantity, density: Quantity = 7850 * kg / m**3, gamma_s: float = 1.15
    ):
        super().__init__(name)
        self._f_y = f_y
        self._density = density

    @property
    def f_y(self) -> Quantity:
        return self._f_y

    @property
    def density(self) -> Quantity:
        return self._density

    

@dataclass
class SteelBar(Steel):
    _E_s: Quantity = field(default=200 * GPa)
    _epsilon_y: Quantity = field(init=False)

    def __init__(
        self, name: str, f_y: Quantity, density: Quantity = 7850 * kg / m**3
    ):
        super().__init__(name, f_y, density)
        self._epsilon_y = f_y.to("MPa") / (self._E_s.to("MPa"))  # 21.2.2.1 - Page 392

    @property
    def E_s(self) -> Quantity:
        return self._E_s

    @property
    def epsilon_y(self) -> Quantity:
        return self._epsilon_y

    def get_properties(self) -> Dict[str, Any]:
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
    _f_u: Quantity = field(default=1860 * MPa)
    _E_s: Quantity = field(default=190000 * MPa)
    prestress_stress: Quantity = field(default=0 * MPa)
    _epsilon_y: Quantity = field(init=False)

    def __init__(
        self, name: str, f_y: Quantity, density: Quantity = 7850 * kg / m**3
    ):
        super().__init__(name, f_y, density)
        self._epsilon_y = self._f_y / self._E_s

    def get_properties(self) -> Dict[str, Any]:
        properties = {
            "E_s": self._E_s.to("MPa"),
            "f_y": self._f_y.to("MPa"),
            "f_u": self._f_u.to("MPa"),
        }
        return properties

    @property
    def f_u(self) -> Quantity:
        return self._f_u

    @property
    def E_s(self) -> Quantity:
        return self._E_s

    @property
    def epsilon_y(self) -> Quantity:
        return self._epsilon_y

    def __str__(self) -> str:
        """Customize the string representation for user-friendly display."""
        properties = self.get_properties()
        return (
            f"SteelStrand Properties ({self.name}):\n"
            f"  f_y: {properties['f_y']}\n"
            f"  f_u: {properties['f_u']}\n"
            f"  E_s: {properties['E_s']}\n"
            f"  epsilon_y: {self.epsilon_y.magnitude:.4f}\n"
            f"  Prestress Stress: {self.prestress_stress}\n"
            f"  Density: {self.density}"
        )
