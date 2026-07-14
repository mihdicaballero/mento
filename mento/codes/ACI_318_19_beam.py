import math
from pint import Quantity
import pandas as pd
import numpy as np
from typing import TYPE_CHECKING, Dict, Any, cast
import warnings
# from devtools import debug

from mento.material import Concrete_ACI_318_19
from mento.rebar import Rebar, RebarDesignInfeasibleError
from mento.units import MPa, mm, kN, inch, ksi, psi, cm, m, kNm, ft, kip, dimensionless
from mento.forces import Forces


if TYPE_CHECKING:
    from ..beam import RectangularBeam  # Import Beam for type checking only


def _initialize_variables_ACI_318_19(self: "RectangularBeam", M_y: Quantity) -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        self._M_u = M_y
        if self._M_u > 0 * kNm:
            self._M_u_bot = self._M_u
            self._M_u_top = 0 * kNm
        else:
            self._M_u_bot = 0 * kNm
            self._M_u_top = self._M_u
        self.f_yt = _calculate_f_yt_aci(self)
        # Consider bottom or top tension reinforcement
        self._A_s_tension = self._A_s_bot if self._M_u >= 0 * kNm else self._A_s_top


##########################################################
# SHEAR CHECK AND DESIGN
##########################################################


def _calculate_shear_reinforcement_aci(self: "RectangularBeam") -> None:
    V_s = self._A_v * self.f_yt * self._d_shear  # Shear contribution of reinforcement
    if isinstance(self.concrete, Concrete_ACI_318_19):
        self._phi_V_s = self.concrete.phi_v * V_s  # Reduced shear contribution of reinforcement


def _calculate_effective_shear_area_aci(self: "RectangularBeam") -> None:
    self._A_cv = self.width * self._d_shear  # Effective shear area
    self._rho_w = self._A_s_tension.to("cm**2") / self._A_cv.to("cm**2")  # Longitudinal reinforcement ratio
    if self.concrete.unit_system == "metric":
        self._lambda_s = math.sqrt(2 / (1 + 0.004 * self._d_shear / mm))
    else:
        self._lambda_s = math.sqrt(2 / (1 + self._d_shear / (10 * inch)))


def _calculate_concrete_shear_strength_aci(self: "RectangularBeam") -> None:
    f_c = self.concrete.f_c
    self._sigma_Nu = min(self._N_u / (6 * self.A_x), 0.05 * f_c)  # Axial stress influence
    if isinstance(self.concrete, Concrete_ACI_318_19):
        if self.concrete.unit_system == "metric":
            V_cmin = 0 * kN

            if self._A_v < self._A_v_min or self._A_v_min == 0 * cm**2 / m:
                if self._A_s_tension == 0 * cm**2:
                    warnings.warn(
                        "Longitudinal rebar As cannot be zero if A_v is less than A_v_min.",
                        UserWarning,
                    )
                self._k_c_min = (
                    0.66
                    * self._lambda_s
                    * self.concrete.lambda_factor
                    * self._rho_w ** (1 / 3)
                    * math.sqrt(f_c / MPa)
                    * MPa
                    + self._sigma_Nu
                )
            else:
                self._k_c_min = max(
                    0.17 * self.concrete.lambda_factor * math.sqrt(f_c / MPa) * MPa + self._sigma_Nu,
                    0.66 * self.concrete.lambda_factor * self._rho_w ** (1 / 3) * math.sqrt(f_c / MPa) * MPa
                    + self._sigma_Nu,
                )
        else:
            V_cmin = 0 * kip
            if self._A_v < self._A_v_min or self._A_v_min == 0 * inch**2 / ft:
                self._k_c_min = (
                    8
                    * self._lambda_s
                    * self.concrete.lambda_factor
                    * self._rho_w ** (1 / 3)
                    * math.sqrt(f_c.to("psi").magnitude)
                    * psi
                    + self._sigma_Nu
                )
            else:
                self._k_c_min = max(
                    2 * self.concrete.lambda_factor * math.sqrt(f_c.to("psi").magnitude) * psi + self._sigma_Nu,
                    8 * self.concrete.lambda_factor * self._rho_w ** (1 / 3) * math.sqrt(f_c / psi) * psi
                    + self._sigma_Nu,
                )
        # Maximum concrete shear strength
        if self.concrete.unit_system == "metric":
            V_cmax = (0.42 * self.concrete.lambda_factor * math.sqrt(self.concrete.f_c / MPa) * MPa) * self._A_cv
        else:
            V_cmax = (5 * self.concrete.lambda_factor * math.sqrt(self.concrete.f_c / psi) * psi) * self._A_cv
        self.V_c = min(V_cmax, max(V_cmin, self._k_c_min * self._A_cv))
        self._phi_V_c = self.concrete.phi_v * self.V_c


def _calculate_max_shear_capacity_aci(self: "RectangularBeam") -> None:
    "Formula for maximum total shear capacity (V_max)"
    if isinstance(self.concrete, Concrete_ACI_318_19):
        if self.concrete.unit_system == "metric":
            V_max = (
                self.V_c + (0.66 * self.concrete.lambda_factor * math.sqrt(self.concrete.f_c / MPa) * MPa) * self._A_cv
            )
        else:
            V_max = self.V_c + (8 * self.concrete.lambda_factor * math.sqrt(self.concrete.f_c / psi) * psi) * self._A_cv
        self._phi_V_max = self.concrete.phi_v * V_max
        self._max_shear_ok = self._V_u < self._phi_V_max


def _calculate_A_v_min_ACI(self: "RectangularBeam", f_c: Quantity) -> None:
    """Calculate the minimum shear reinforcement based on unit system."""
    # 'Minimum reinforcement should be placed if the factored shear Vu
    # is greater than half the shear capacity of the concrete,
    # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
    # Rebar needed, V_u > φ_v*V_c/2 for Imperial system
    f_yt = _calculate_f_yt_aci(self)

    if self.concrete.unit_system == "metric":
        self._A_v_min = max(
            (0.062 * math.sqrt(f_c / MPa) * MPa / f_yt) * self.width,
            (0.35 * MPa / f_yt) * self.width,
        )
    else:
        self._A_v_min = max(
            (0.75 * math.sqrt(f_c / psi) * psi / f_yt) * self.width,
            (50 * psi / f_yt) * self.width,
        )


def _calculate_f_yt_aci(self: "RectangularBeam") -> Quantity:
    """Determine the yield strength of steel based on unit system."""
    if self.concrete.unit_system == "metric":
        return min(self.steel_bar.f_y, 420 * MPa)
    else:
        return min(self.steel_bar.f_y, 60 * ksi)


def _check_minimum_reinforcement_requirement_aci(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        if self.concrete.unit_system == "metric":
            if (
                self._V_u
                < 0.083
                * self.concrete.phi_v
                * self.concrete.lambda_factor
                * math.sqrt(self.concrete.f_c / MPa)
                * MPa
                * self._A_cv
            ):
                self._A_v_req = 0 * cm**2 / m
                self._A_v_min = 0 * cm**2 / m
                self._max_shear_ok = True
            elif (
                0.083
                * self.concrete.phi_v
                * self.concrete.lambda_factor
                * math.sqrt(self.concrete.f_c / MPa)
                * MPa
                * self._A_cv
                < self._V_u
                < self._phi_V_max
            ):
                _calculate_A_v_min_ACI(self, self.concrete.f_c)
                self._max_shear_ok = True
            else:
                _calculate_A_v_min_ACI(self, self.concrete.f_c)
                self._max_shear_ok = False
        else:
            if (
                self._V_u
                < self.concrete.phi_v
                * self.concrete.lambda_factor
                * math.sqrt(self.concrete.f_c / psi)
                * psi
                * self._A_cv
            ):
                self._A_v_req = 0 * inch**2 / ft
                self._A_v_min = 0 * inch**2 / ft
                self._max_shear_ok = True
                # if self._V_u != 0*kN:
                #     self._stirrup_d_b = 0*inch
            elif (
                self.concrete.phi_v
                * self.concrete.lambda_factor
                * math.sqrt(self.concrete.f_c / psi)
                * psi
                * self._A_cv
                < self._V_u
                < self._phi_V_max
            ):
                _calculate_A_v_min_ACI(self, self.concrete.f_c)
                self._max_shear_ok = True
            else:
                _calculate_A_v_min_ACI(self, self.concrete.f_c)
                self._max_shear_ok = False


def _calculate_V_s_req(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        self._V_s_req = self._V_u - self._phi_V_c
        self._A_v_req = max(
            self._V_s_req / (self.concrete.phi_v * self.f_yt * self._d_shear),
            self._A_v_min,
        ).to("cm ** 2 / m")


def _calculate_total_shear_strength_aci(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        self._phi_V_n = self.concrete.phi_v * (self.V_c + self._A_v * self.f_yt * self._d_shear)
        V_d_max = min(self._phi_V_n, self._phi_V_max)
        self._DCRv = abs((self._V_u.to("kN").magnitude / V_d_max.to("kN").magnitude))


def _calculate_rebar_spacing_aci(self: "RectangularBeam") -> None:
    section_rebar = Rebar(self)
    n_legs_actual = self._stirrup_n * 2  # Ensure legs are even
    self._stirrup_s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs_actual - 1)
    (
        self._stirrup_s_max_l,
        self._stirrup_s_max_w,
    ) = section_rebar.calculate_max_spacing_ACI_318_19(self._V_u - self._phi_V_c, self._A_cv)
    self._stirrup_s_l = max(self._stirrup_s_l, 0 * inch)
    self._stirrup_s_w = max(self._stirrup_s_w, 0 * inch)


def _check_shear_ACI_318_19(self: "RectangularBeam", force: Forces) -> pd.DataFrame:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        # Set the initial variables
        _initialize_variables_ACI_318_19(self, force.M_y)
        self._N_u = force._N_x
        self._V_u = force._V_z

        # Minimum shear reinforcement calculation
        _calculate_A_v_min_ACI(self, self.concrete.f_c)
        if self._stirrup_n > 0:
            # Shear reinforcement calculations
            _calculate_shear_reinforcement_aci(self)
        else:
            # Set current stirrup diameter to zero
            self._stirrup_d_b = 0 * mm
            self._update_effective_heights()

        # Effective shear area and longitudinal reinforcement ratio
        _calculate_effective_shear_area_aci(self)

        # Check if minimum reinforcement is required
        _check_minimum_reinforcement_requirement_aci(self)

        # Concrete shear strength calculation
        _calculate_concrete_shear_strength_aci(self)

        # Maximum total shear capacity
        _calculate_max_shear_capacity_aci(self)

        # Calculate required shear reinforcement
        _calculate_V_s_req(self)

        # Total shear strength
        _calculate_total_shear_strength_aci(self)

        # Rebar spacing checks
        _calculate_rebar_spacing_aci(self)

        # Check results and return DataFrame
        results = _compile_results_ACI_shear(self, force)
        _initialize_dicts_ACI_318_19_shear(self)
        return results


def _design_shear_ACI_318_19(self: "RectangularBeam", force: Forces) -> None:
    # Set the initial variables
    _initialize_variables_ACI_318_19(self, force.M_y)
    self._N_u = force._N_x
    self._V_u = force._V_z
    # Minimum shear reinforcement calculation
    _calculate_A_v_min_ACI(self, self.concrete.f_c)
    # Consider that the beam has minimum reinforcement
    self._A_v = self._A_v_min
    # Effective shear area and longitudinal reinforcement ratio
    _calculate_effective_shear_area_aci(self)
    # Concrete shear strength calculation
    _calculate_concrete_shear_strength_aci(self)
    # Maximum total shear capacity
    _calculate_max_shear_capacity_aci(self)
    # Check if minimum reinforcement is required
    _check_minimum_reinforcement_requirement_aci(self)
    # Calculate required shear reinforcement
    _calculate_V_s_req(self)
    # Update spacing of longitudinal reinforcement calculation
    self._update_longitudinal_rebar_attributes()

    return None


def _compile_results_ACI_shear(self: "RectangularBeam", force: Forces) -> pd.DataFrame:
    results = {
        "Label": self.label,
        "Comb.": force.label,
        "Av,min": round(self._A_v_min.to("cm²/m").magnitude, 2),
        "Av,req": round(self._A_v_req.to("cm²/m").magnitude, 2),
        "Av": round(self._A_v.to("cm²/m").magnitude, 2),
        "Vu": self._V_u.to("kN").magnitude,
        "Nu": self._N_u.to("kN").magnitude,
        "ØVc": round(self._phi_V_c.to("kN").magnitude, 2),
        "ØVs": round(self._phi_V_s.to("kN").magnitude, 2),
        "ØVn": round(self._phi_V_n.to("kN").magnitude, 2),
        "ØVmax": round(self._phi_V_max.to("kN").magnitude, 2),
        "Vu≤ØVmax": self._max_shear_ok,
        "Vu≤ØVn": self._V_u <= self._phi_V_n,
        "DCR": round(self._DCRv, 3),
    }
    return pd.DataFrame([results], index=[0])


# TODO: Delete this method since is not used

# def _calculate_phi_ACI_318_19(self: "RectangularBeam", epsilon_most_strained: float) -> float:
#     """
#     Calculates the strength reduction factor (φ) for flexural design
#     based on ACI 318-19.
#     It is used for columns; for beams, it is not required since beams
#     are always designed to be tension-controlled, with φ=0.9.

#     Parameters:
#         epsilon_most_strained (float): Strain in the most strained steel fiber.

#     Returns:
#         float: The strength reduction factor (φ), ranging from 0.65 to 0.9.

#     Description:
#         - φ = 0.65 if ε_most_strained ≤ ε_y (yield strain).
#         - φ transitions linearly from 0.65 to 0.9 if ε_y < ε_most_strained ≤ ε_y + ε_c
#         (concrete crushing strain).
#         - φ = 0.9 if ε_most_strained > ε_y + ε_c.
#     """
#     # Retrieve concrete crushing strain (ε_c)
#     epsilon_c = self.concrete.get_properties()["epsilon_c"]

#     # Calculate φ based on ε_most_strained
#     if epsilon_most_strained <= self.steel_bar.epsilon_y:
#         return 0.65
#     elif epsilon_most_strained <= self.steel_bar.epsilon_y + epsilon_c:
#         return (0.9 - 0.65) * (epsilon_most_strained - self.steel_bar.epsilon_y) / epsilon_c + 0.65
#     else:
#         return 0.9


##########################################################
# FLEXURE CHECK AND DESIGN
##########################################################


def _maximum_flexural_reinforcement_ratio_ACI_318_19(self: "RectangularBeam") -> float:
    """
    Calculates the maximum flexural reinforcement ratio (ρ_max) according to the
    ACI 318-19 design code.

    Returns:
        float: The maximum reinforcement ratio (ρ_max) for the section,
        or 0 if the design code is not ACI 318-19.

    Description:
        This function determines the maximum reinforcement ratio (ρ_max)
        allowed by the ACI 318-19 design code to ensure ductile behavior
        of reinforced concrete sections. The calculation depends on the
        properties of the concrete (β1 and ε_c) and steel (ε_y).

        Note:
        - This function only works if the `design_code` of the concrete is
        set to "ACI 318-19". For other codes, it returns 0.

    """
    # Cast the concrete object to the specific ACI subclass
    concrete_aci = cast("Concrete_ACI_318_19", self.concrete)

    # Calculate maximum reinforcement ratio (ρ_max)
    # From CRSI Guide - Beam Theory page 6-3
    rho_max = (
        0.85
        * concrete_aci.beta_1
        * self.concrete.f_c
        / self.steel_bar.f_y
        * (concrete_aci._epsilon_c / (self.steel_bar.epsilon_y + concrete_aci._epsilon_c * 2))
    )

    return rho_max


def _c_neutral_axis_at_ductility_limit_ACI_318_19(
    self: "RectangularBeam", d: Quantity
) -> Quantity:
    """
    Neutral axis depth at the ACI 318-19 tension-controlled boundary.

    At the ductility limit, eps_t_min = eps_y + eps_cu. From strain
    compatibility (c/d = eps_cu / (eps_cu + eps_t)):
        c_t = 0.003 * d / (eps_y + 0.006)

    Sections with c < c_t are ductile (tension-controlled); with c > c_t are
    over-reinforced (brittle failure).
    """
    return 0.003 * d / (self.steel_bar.epsilon_y + 0.006)


def _f_s_prime_net_at_ductility_limit_ACI_318_19(
    self: "RectangularBeam", d: Quantity, d_prime: Quantity
) -> Quantity:
    """
    Effective compression-steel stress at the ductility limit, corrected for
    displaced concrete.

    Evaluated at c_t (neutral axis at the ductility limit):
        eps_s' = (c_t - d') / c_t * eps_cu
        f_s'   = min(eps_s' * E_s, f_y)
        f_s'_net = f_s' - 0.85 * f_c

    The `- 0.85 * f_c` accounts for the concrete displaced by the bar: that
    volume was already contributing to equilibrium via the 0.85·f_c·a·b
    block, so it cannot be counted twice.
    """
    c_t = _c_neutral_axis_at_ductility_limit_ACI_318_19(self, d)
    f_s_prime = min(0.003 * self.steel_bar.E_s * (1 - d_prime / c_t), self.steel_bar.f_y)
    return f_s_prime - 0.85 * self.concrete.f_c


def _minimum_flexural_reinforcement_ratio_ACI_318_19(self: "RectangularBeam", M_u: Quantity) -> Quantity:
    """
    Calculates the minimum flexural reinforcement ratio according to ACI 318-19
    provisions based on the factored moment, M_u.

    This method determines the minimum amount of tensile reinforcement
    (in terms of a reinforcement ratio) that should be provided in a
    reinforced concrete section according to ACI 318-19. If the factored
    moment M_u is zero, it means there is no flexural demand, and hence
    no minimum flexural reinforcement is required (the ratio is zero).
    If M_u is not zero, the method checks the unit system (metric or imperial)
    and computes the required minimum ratio accordingly. These calculations
    depend on the compressive strength of the concrete (f_c) and the yield
    strength of the reinforcing steel (f_y).

    Parameters
    ----------
    M_u : Quantity
        The factored moment for the section where the minimum flexural
        reinforcement ratio is required. The unit should be consistent with
        the chosen system (e.g., kNm in metric).

    Returns
    -------
    Quantity
        The minimum flexural reinforcement ratio (dimensionless).

    Notes
    -----
    - If M_u = 0, it indicates no moment demand, thus no minimum flexural
    reinforcement is required (resulting in a zero ratio).
    - For M_u > 0, the minimum ratio is determined using formulas involving
    the square root of f_c and the value of f_y, in accordance with
    ACI 318-19.
    - The result is a dimensionless ratio representing the minimum area
    of steel to the area of the concrete section.

    References
    ----------
    ACI Committee 318. "Building Code Requirements for Structural Concrete
    (ACI 318-19) and Commentary", American Concrete Institute, 2019.
    """
    if M_u == 0 * kNm:
        minimum_ratio = 0 * dimensionless
    else:
        if self.concrete.unit_system == "metric":
            minimum_ratio = max(
                (0.25 * np.sqrt(self.concrete.f_c / MPa) * MPa / self.steel_bar.f_y),
                (1.4 * MPa / self.steel_bar.f_y),
            )
        else:
            minimum_ratio = max(
                (3 * np.sqrt(self.concrete.f_c / psi) * psi / self.steel_bar.f_y),
                (200 * psi / self.steel_bar.f_y),
            )
    return minimum_ratio


def _calculate_flexural_reinforcement_ACI_318_19(
    self: "RectangularBeam", M_u: Quantity, d: float, d_prima: float
) -> tuple[Quantity, Quantity, Quantity, Quantity, float, bool]:
    """
    Calculates the flexural reinforcement for a given factored moment according to ACI 318-19.

    This function computes the required reinforcement areas (minimum, maximum, and final) and
    the compression reinforcement (if required) for a given factored moment. The moment M_u must
    always be provided as a positive value. For a positive moment, pass 'd' as the effective depth
    of the tensile reinforcement and 'd_prima' as the effective depth of the compression reinforcement.
    For a negative moment, reverse the roles of 'd' and 'd_prima'.

    Parameters:
        M_u (Quantity): The factored moment (always a positive value).
        d (float): Effective depth of the tensile reinforcement.
        d_prima (float): Effective depth of the compression reinforcement.

    Returns:
        tuple: A tuple containing:
            - A_s_min (Quantity): Minimum reinforcement area required by the code.
            - A_s_max (Quantity): Maximum reinforcement area allowed by the code.
            - A_s_final (Quantity): Final reinforcement area adopted for the tensile zone.
            - A_s_comp (Quantity): Compression reinforcement area (if required).
            - c_d (float): Ratio of the calculated neutral axis depth to the effective depth (c/d).
            - A_s_bool: Boolean indicating if 4/3*A_s_calc is adopted instead of A_s_min
    """
    concrete_aci = cast("Concrete_ACI_318_19", self.concrete)
    # Extract relevant properties and settings
    b = self.width

    # Determine minimum and maximum reinforcement areas
    rho_min = _minimum_flexural_reinforcement_ratio_ACI_318_19(self, M_u)  # Mechanical minimum reinforcement ratio
    A_s_min = rho_min * d * b
    rho_max = _maximum_flexural_reinforcement_ratio_ACI_318_19(
        self,
    )
    A_s_max = rho_max * d * b

    # Calculate required reinforcement based on the nominal moment capacity
    R_n = M_u / (concrete_aci._phi_t * b * d**2)
    # Verify if the value under the square root is negative
    sqrt_value = 1 - 2 * R_n / (0.85 * self.concrete.f_c)
    if sqrt_value < 0:
        # Here we assign A_s_max so that the calculation does not break,
        # resulting in a DCR greater than 1.
        A_s_calc = A_s_max
    else:
        A_s_calc = 0.85 * self.concrete.f_c * b * d / self.steel_bar.f_y * (1 - np.sqrt(sqrt_value))

    # Calculate the neutral axis depth based on equilibrium: 0.85 * f_c * c * beta_1 * b = A_s * f_y
    c = A_s_calc * self.steel_bar.f_y / (0.85 * self.concrete.f_c * b * concrete_aci._beta_1)

    # Helper function to clean near-zero values
    def clean_zero(value: float, tolerance: float = 1e-6) -> float:
        return 0.0 if abs(value) < tolerance else value

    c_d = clean_zero(c / d)

    A_s_final: Quantity = 0 * cm**2

    A_s_bool = False

    # 1.8‰ of the gross section (custom geometric minimum rule)
    if self.mode == "beam":
        A_s_geo_min = (1.8 / (1000)) * self.width * self.height
    elif self.mode == "slab":
        A_s_geo_min = (1.8 / (1000)) * self.width * self.height

    if A_s_calc >= A_s_min:
        # Case 1:
        # The required steel already exceeds the ACI minimum.
        # Do not apply the 4/3 rule, as it would increase steel unnecessarily.
        A_s_final = A_s_calc
    else:
        # Case 2:
        # The required steel is less than the ACI minimum.
        # Evaluate whether using 4/3 * A_s_calc provides less steel than A_s_min.
        A_s_4_3 = (4 * A_s_calc) / 3

        if A_s_4_3 < A_s_min:
            # The 4/3 rule is potentially beneficial → check the geometric minimum
            if A_s_4_3 >= A_s_geo_min:
                # 4/3 * A_s_calc satisfies the geometric minimum → adopt it
                A_s_final = A_s_4_3
                A_s_bool = True
            else:
                # 4/3 * A_s_calc is lower than 1.8‰ of (b·h) → enforce geometric minimum
                A_s_final = A_s_geo_min
                # A_s_bool remains False (4/3 rule not effectively used)
        else:
            # 4/3 * A_s_calc is not smaller than A_s_min → use the ACI minimum
            A_s_final = A_s_min

    # Optional cleanup and unit formatting
    A_s_final = clean_zero(A_s_final.to("cm**2").magnitude) * cm**2

    # Determine if compression reinforcement is required
    if A_s_final <= A_s_max:
        A_s_comp = 0 * cm**2
    else:
        self._doubly_reinforced = True
        rho = (
            0.85
            * concrete_aci._beta_1
            * self.concrete.f_c.to(psi).magnitude
            / self.steel_bar.f_y.to(psi).magnitude
            * (0.003 / (self.steel_bar.epsilon_y + 0.006))
        )
        # Whitney block depth at the ductility limit (exact, not approximated):
        #     a_max = beta_1 * c_t
        # This replaces the textbook shortcut d - 0.59 * rho * fy * d / fc,
        # which relies on 0.59 ≈ 1/1.7 and drops ~0.3% of precision.
        c_t = _c_neutral_axis_at_ductility_limit_ACI_318_19(self, d)
        c_d = clean_zero(c_t / d)
        a_max = concrete_aci._beta_1 * c_t
        M_n_t = rho * self.steel_bar.f_y * (d - a_max / 2) * b * d
        M_n_prima = M_u / concrete_aci._phi_t - M_n_t
        f_s_prima_net = _f_s_prime_net_at_ductility_limit_ACI_318_19(self, d, d_prima)
        A_s_comp = M_n_prima / (f_s_prima_net * (d - d_prima))
        A_s_final = rho * b * d + A_s_comp * f_s_prima_net / self.steel_bar.f_y

    return A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool


def _determine_nominal_moment_simple_reinf_ACI_318_19(self: "RectangularBeam", A_s: Quantity, d: Quantity) -> Quantity:
    """
    Determines the nominal moment for a simply reinforced section according to ACI 318-19.

    This formula is used ONLY when the provided reinforcement area (A_s) is less than or equal to A_s_max.

    The equilibrium of forces is assumed (compression equals tension):
        0.85 * f_c * a * b = A_s * f_y
    which implies:
        a = (A_s * f_y) / (0.85 * f_c * b)

    Parameters:
        A_s (Quantity): The area of reinforcement.
        d (Quantity): The effective depth of the section.

    Returns:
        Quantity: The nominal moment (M_n) calculated as A_s * f_y * (d - a/2).
    """
    # Calculate the depth of the equivalent rectangular stress block (a)
    a = A_s * self.steel_bar.f_y / (0.85 * self.concrete.f_c * self.width)
    # Calculate the nominal moment (M_n)
    M_n = A_s * self.steel_bar.f_y * (d - a / 2)
    return M_n


def _determine_nominal_moment_double_reinf_ACI_318_19(
    self: "RectangularBeam",
    A_s: Quantity,
    d: Quantity,
    d_prime: Quantity,
    A_s_prime: Quantity,
) -> Quantity:
    """
    Determines the nominal moment for a doubly reinforced beam section according to ACI 318-19.

    This method is used only when the beam has reinforcement exceeding the maximum limit
    and includes compression reinforcement.

    Equilibrium is assumed (Compression = Tension):
        C_total = Concrete compression (Cc) + Compression reinforcement (Cs)
        Tension (T) = A_s * f_y

    Initially, it is assumed that the compression reinforcement yields (i.e., f_s_prime = f_y),
    so the equilibrium equation becomes:
        0.85 * f_c * a * b + A_s_prime * f_y = A_s * f_y

    Parameters:
        A_s (Quantity): Area of the tensile reinforcement.
        d (Quantity): Effective depth of the beam section.
        d_prime (Quantity): Effective depth (cover) to the compression reinforcement.
        A_s_prime (Quantity): Area of the compression reinforcement.

    Returns:
        Quantity: The nominal moment (M_n) of the doubly reinforced section.
    """
    f_c = self.concrete.f_c
    if isinstance(self.concrete, Concrete_ACI_318_19):
        beta_1 = self.concrete._beta_1
        epsilon_c = self.concrete._epsilon_c

    f_y = self.steel_bar.f_y
    E_s = self.steel_bar._E_s
    epsilon_y = self.steel_bar._epsilon_y
    b = self.width

    # -------------------------------------------------------------------------
    # Step 1: Assume that the compression steel is yielding (f_s_prime = f_y).
    # Based on equilibrium WITH displaced-concrete correction:
    #   A_s * f_y = 0.85 * f_c * a * b + A_s_prime * (f_y - 0.85 * f_c)
    # The A_s_prime bar physically occupies a volume where the 0.85*f_c*a*b
    # block already accounts for compression; the effective contribution of the
    # compression steel is therefore (f_y - 0.85*f_c), not f_y.
    # Solving for the neutral axis depth 'c_assumed' (a = beta_1 * c):
    c_assumed = (A_s * f_y - A_s_prime * (f_y - 0.85 * f_c)) / (0.85 * f_c * b * beta_1)

    # Compute the strain in the compression reinforcement with the assumed c value:
    epsilon_s = (c_assumed - d_prime) / c_assumed * epsilon_c if c_assumed > 0 else 0
    # 0 is an arbitrary value, for example, when As top and Bottom are equal, c_assumed is 0,
    # and compression steel is not yielding, so epsilon_s must be calculated.

    # -------------------------------------------------------------------------
    # Step 2: Check if the assumed compression steel strain exceeds the yield strain.
    if epsilon_s >= epsilon_y:
        # The assumption is valid (compression reinforcement yields).
        # M_n uses (f_y - 0.85*f_c) for the compression-steel contribution to
        # avoid double-counting the displaced concrete already included in the
        # 0.85*f_c*a*b block.
        a_assumed = c_assumed * beta_1
        M_n = (
            0.85 * f_c * a_assumed * b * (d - a_assumed / 2)
            + A_s_prime * (f_y - 0.85 * f_c) * (d - d_prime)
        )
        return M_n
    else:
        # The assumption is invalid, so determine the actual neutral axis depth
        # 'c' using the quadratic equation.
        # Based on equilibrium WITH displaced-concrete correction:
        #   A_s * f_y = 0.85 * f_c * b * beta_1 * c
        #             + A_s_prime * [ epsilon_c * E_s * (c - d_prime) / c
        #                             - 0.85 * f_c ]
        # Multiplying by c and rearranging into A*c^2 + B*c + C = 0:
        A = 0.85 * f_c * b * beta_1
        B = A_s_prime * (epsilon_c * E_s - 0.85 * f_c) - A_s * f_y
        C = -d_prime * A_s_prime * epsilon_c * E_s

        # Solve for c using the quadratic formula:
        c = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)

        # Compute the corresponding depth of the equivalent rectangular stress block:
        a = c * beta_1

        # Recompute the strain and stress in the compression reinforcement,
        # then apply the displaced-concrete correction:
        epsilon_s_prime = (c - d_prime) / c * epsilon_c
        f_s_prime = epsilon_s_prime * E_s
        f_s_prime_net = f_s_prime - 0.85 * f_c

        # Determine the adjusted areas for tension reinforcement:
        # a portion of the tensile reinforcement is balanced by the (net)
        # compression reinforcement contribution.
        A_s_2 = A_s_prime * f_s_prime_net / f_y
        A_s_1 = A_s - A_s_2

        # Calculate the nominal moment contributions:
        M_n_1 = A_s_1 * f_y * (d - a / 2)
        M_n_2 = A_s_prime * f_s_prime_net * (d - d_prime)
        M_n = M_n_1 + M_n_2

        return M_n


def _determine_nominal_moment_ACI_318_19(self: "RectangularBeam", force: Forces) -> None:
    """
    Determines the nominal moment for a given section with both top and bottom reinforcement,
    calculating the nominal moment for both positive and negative moment scenarios.

    For positive moments, the tension is in the bottom reinforcement.
    For negative moments, the tension is in the top reinforcement.

    Parameters:
        force (Forces): An object containing the forces acting on the section, including the moment M_y.

    Returns:
        None
    """

    # Calculate minimum and maximum reinforcement ratios
    rho_min = _minimum_flexural_reinforcement_ratio_ACI_318_19(self, force._M_y)
    rho_max = _maximum_flexural_reinforcement_ratio_ACI_318_19(self)

    # For positive moments (tension in the bottom), set minimum reinforcement accordingly.
    if force._M_y > 0 * kNm:
        rho_min_top = 0 * dimensionless
        rho_min_bot = rho_min
    else:
        rho_min_top = rho_min
        rho_min_bot = 0 * dimensionless

    # Calculate minimum and maximum bottom reinforcement areas
    self._A_s_min_bot = rho_min_bot * self._d_bot * self.width
    self._A_s_max_bot = rho_max * self._d_bot * self.width

    # Determine the nominal moment for positive moments.
    # Three-branch dispatcher based on ductility:
    #   1. Ductile section (A_s_bot <= A_s_max_bot): simple flexure with A_s real.
    #   2. Over-reinforced but redeemed by compression steel: doubly reinforced
    #      formula with A_s real. When A_s_top = 0, A_s_max_total collapses to
    #      A_s_max_bot and this branch cannot be reached (goes to branch 3).
    #   3. Over-reinforced beyond redemption: cap A_s to A_s_max_total and use
    #      doubly reinforced formula. If A_s_top = 0, the doubly reinforced
    #      formula degenerates to simple with A_s_max_bot.
    if self._A_s_bot <= self._A_s_max_bot:
        M_n_positive = _determine_nominal_moment_simple_reinf_ACI_318_19(
            self, self._A_s_bot, self._d_bot
        )
    else:
        # Compression steel contribution at the ductility limit determines how
        # much the tension-steel cap can be extended:
        #     A_s_max_total = A_s_max_bot + A_s_top * f_s'_net / f_y
        f_s_prima_net = _f_s_prime_net_at_ductility_limit_ACI_318_19(
            self, self._d_bot, self._c_mec_top
        )
        A_s_max_total = (
            self._A_s_max_bot + self._A_s_top * f_s_prima_net / self.steel_bar.f_y
        )
        A_s_eff = self._A_s_bot if self._A_s_bot <= A_s_max_total else A_s_max_total
        M_n_positive = _determine_nominal_moment_double_reinf_ACI_318_19(
            self, A_s_eff, self._d_bot, self._c_mec_top, self._A_s_top
        )

    # Determine capacity for negative moment (tension at the top)
    self._A_s_min_top = rho_min_top * self._d_top * self.width
    self._A_s_max_top = rho_max * self._d_top * self.width

    # Determine the nominal moment for negative moments (tension on top face).
    # Mirror of the positive-moment dispatcher, plus the leading edge case of
    # zero tension steel (no top → no negative capacity).
    if self._A_s_top == 0 * cm**2:
        M_n_negative = 0 * kNm
    elif self._A_s_top <= self._A_s_max_top:
        M_n_negative = _determine_nominal_moment_simple_reinf_ACI_318_19(
            self, self._A_s_top, self._d_top
        )
    else:
        f_s_prima_net = _f_s_prime_net_at_ductility_limit_ACI_318_19(
            self, self._d_top, self._c_mec_bot
        )
        A_s_max_total = (
            self._A_s_max_top + self._A_s_bot * f_s_prima_net / self.steel_bar.f_y
        )
        A_s_eff = self._A_s_top if self._A_s_top <= A_s_max_total else A_s_max_total
        M_n_negative = _determine_nominal_moment_double_reinf_ACI_318_19(
            self, A_s_eff, self._d_top, self._c_mec_bot, self._A_s_bot
        )

    # Calculate the design moment capacities for both bottom and top reinforcement
    concrete_aci = cast("Concrete_ACI_318_19", self.concrete)
    self._phi_M_n_bot = concrete_aci._phi_t * M_n_positive
    self._phi_M_n_top = concrete_aci._phi_t * M_n_negative

    return None


def _check_flexure_ACI_318_19(self: "RectangularBeam", force: Forces) -> pd.DataFrame:
    """
    Checks the flexural capacity of the section according to ACI 318-19 guidelines.

    This function accepts a single force and performs the flexural check of the section
    following the ACI 318-19 requirements. It initializes the design variables, computes the
    nominal moments for both top and bottom reinforcement, determines the required reinforcement
    areas, and calculates the design capacity ratios. Finally, the results are compiled into a
    Pandas DataFrame.

    Parameters:
        force (Forces): The force acting on the section, which must include a single moment value.

    Returns:
        pd.DataFrame: A DataFrame containing the flexural design metrics and results.
    """

    # Initialize the design variables according to ACI 318-19 requirements using the provided force.
    _initialize_variables_ACI_318_19(self, force.M_y)

    # Calculate the nominal moments for both top and bottom reinforcement.
    _determine_nominal_moment_ACI_318_19(self, force)

    if self._M_u >= 0:
        # For positive moments, calculate the reinforcement requirements for the bottom tension side.
        (
            self._A_s_min_bot,
            self._A_s_max_bot,
            self._A_s_req_bot,
            self._A_s_req_top,
            self._c_d_bot,
            self._A_s_bool_bot,
        ) = _calculate_flexural_reinforcement_ACI_318_19(
            self,
            self._M_u_bot,
            self._d_bot,
            self._c_mec_top,
        )
        self._c_d_top = 0
        # Calculate the design capacity ratio for the bottom side.
        if self._phi_M_n_bot.to("kN*m").magnitude == 0:
            self._phi_M_n_bot = 0.01 * kNm
        self._DCRb_bot = self._M_u_bot.to("kN*m").magnitude / self._phi_M_n_bot.to("kN*m").magnitude
        self._DCRb_top = 0
    else:
        # For negative moments, calculate the reinforcement requirements for the top tension side.
        (
            self._A_s_min_top,
            self._A_s_max_top,
            self._A_s_req_top,
            self._A_s_req_bot,
            self._c_d_top,
            self._A_s_bool_top,
        ) = _calculate_flexural_reinforcement_ACI_318_19(
            self,
            abs(self._M_u_top.to("kN*m").magnitude) * kNm,
            self._d_top,
            self._c_mec_bot,
        )
        self._c_d_bot = 0
        # Calculate the design capacity ratio for the top side.
        if self._phi_M_n_top.to("kN*m").magnitude == 0:
            self._phi_M_n_top = 0.01 * kNm
        self._DCRb_top = -self._M_u_top.to("kN*m").magnitude / self._phi_M_n_top.to("kN*m").magnitude
        self._DCRb_bot = 0

    # Determine the maximum detailing cover dimensions for top and bottom.
    self._d_b_max_top = max(self._d_b1_t, self._d_b2_t, self._d_b3_t, self._d_b4_t)
    self._d_b_max_bot = max(self._d_b1_b, self._d_b2_b, self._d_b3_b, self._d_b4_b)

    # Calculate the longitudinal reinforcement ratios for both sides.
    self._rho_l_bot = self._A_s_bot / (self._d_bot * self.width)
    self._rho_l_top = self._A_s_bot / (self._d_top * self.width)

    # Compile the design results into a dictionary.
    results = _compile_results_ACI_flexure_metric(self, force)

    # Initialize any additional dictionaries required for ACI 318-19 flexural checks.
    _initialize_dicts_ACI_318_19_flexure(self)

    # Return the results as a Pandas DataFrame.
    return pd.DataFrame([results], index=[0])


# ──────────────────────────────────────────────────────────────────────────────
# Cycle-safe flexural design helpers
# ──────────────────────────────────────────────────────────────────────────────
MAX_FLEXURE_ITERATIONS = 30   # safety net for slow divergence without cycling


def _rebar_design_fingerprint(rebar_design: dict) -> tuple:
    """Canonical, hashable identifier of a discrete rebar layout.

    Two iterations producing the same fingerprint mean the Picard (fixed-point
    iteration) loop has cycled back to a previous configuration. The centroid/spacing are
    intentionally excluded — they are derived from this tuple, not part of
    its identity.
    """
    def _diam_mm(key):
        q = rebar_design.get(key, 0 * mm)
        return float(q.to("mm").magnitude) if q is not None else 0.0
    return (
        int(rebar_design.get("n_1", 0)), _diam_mm("d_b1"),
        int(rebar_design.get("n_2", 0)), _diam_mm("d_b2"),
        int(rebar_design.get("n_3", 0)), _diam_mm("d_b3"),
        int(rebar_design.get("n_4", 0)), _diam_mm("d_b4"),
    )


def _select_safe_design(
    self: "RectangularBeam",
    candidate_designs: list,
    M_demand: Quantity,
    face: str,
) -> dict:
    """Among a set of candidate rebar designs — those visited during the
    Picard (fixed-point iteration) loop — return the most appropriate one for the given face.

    Selection priority:

    1. Among candidates whose actual phi*Mn satisfies phi*Mn >= M_demand,
       return the one with the smallest As (most economical valid layout).
    2. If no candidate satisfies the check, return the one with the largest
       phi*Mn (closest to passing). Downstream `check_flexure` will surface
       DCR > 1 so the user is aware that the section is insufficient.

    Never raises: keeping `design_flexure` total preserves the public
    contract — summary, plot, shear design and check must keep working even
    when the section is underdesigned.

    Parameters
    ----------
    candidate_designs : list of dict
        Each dict is a rebar_designer payload (keys n_1..n_4, d_b1..d_b4,
        total_as, ...).
    M_demand : Quantity
        Required design moment on the face (positive magnitude).
    face : {"bot", "top"}
        Which face the layout governs.
    """
    M_demand_abs = abs(M_demand.to("kN*m"))

    evaluated: list = []   # tuples of (As_provided, phi_Mn, design_dict)
    for design in candidate_designs:
        if face == "bot":
            self._apply_longitudinal_design_bot(design)
        else:
            self._apply_longitudinal_design_top(design)
        # Recompute phi_M_n with the just-applied layout (centroid included)
        probe_force = Forces(M_y=(M_demand_abs if face == "bot" else -M_demand_abs))
        _determine_nominal_moment_ACI_318_19(self, probe_force)
        phi_M_n = self._phi_M_n_bot if face == "bot" else self._phi_M_n_top
        evaluated.append((
            design.get("total_as", 0 * (cm**2)),
            phi_M_n.to("kN*m"),
            design,
        ))

    # Primary criterion: candidates that satisfy the check
    passing = [e for e in evaluated if e[1] >= M_demand_abs]
    if passing:
        passing.sort(key=lambda e: e[0])           # smallest As first
        return passing[0][2]

    # Fallback: no candidate passes — pick the one with the largest phi*Mn
    # (closest to passing). The downstream check will report DCR > 1.
    evaluated.sort(key=lambda e: -e[1])             # largest phi*Mn first
    return evaluated[0][2]


def _design_flexure_ACI_318_19(self: "RectangularBeam", max_M_y_bot: Quantity, max_M_y_top: Quantity) -> Dict[str, Any]:
    """
    Designs the flexural reinforcement for a beam cross-section according to ACI 318-19.
    Implements 'governing-face + reconciliation' so that each face's final layout covers
    tension on that face OR compression from the opposite face, whichever is larger.

    The mechanical cover is solved by a Picard (fixed-point) iteration. The loop
    is bounded by `MAX_FLEXURE_ITERATIONS` and includes per-face cycle detection:
    if the same layout fingerprint reappears, the loop exits and a safe layout
    is picked among the visited candidates by re-evaluating phi*Mn (see
    `_select_safe_design`). This avoids infinite oscillation when two or more
    discrete layouts alternate without converging in the strict tolerance.
    """

    # --- helpers -----------------------------------------------------------------
    def _design_longitudinal_for_area(A_req: Quantity, A_max: Quantity, mech_cover: Quantity):
        """Run discrete design for a target area and return best_design dict, or
        None if the rebar designer cannot fit any combination in the section
        geometry (RebarDesignInfeasibleError). Callers must handle the None
        result — preserves the public contract that design_flexure never
        crashes, delegating the "insufficient section" report to check_flexure
        via DCR>1."""
        rebar = self._create_rebar_designer()
        _ = rebar.longitudinal_rebar_ACI_318_19(A_req, A_max, mech_cover)
        try:
            return rebar.longitudinal_rebar_design
        except RebarDesignInfeasibleError:
            return None

    # --- initial guesses ----------------------------------------------------------
    rec_mec = self.c_c + self._stirrup_d_b + 1 * cm  # bottom mechanical cover to centroid (initial)
    d_prima = self.c_c + self._stirrup_d_b + 1 * cm  # top mechanical cover to centroid (initial)

    tol = 0.01 * cm
    Err = 2 * tol

    # Cycle detection — store the layout payload for each fingerprint seen.
    # Using a dict preserves insertion order (3.7+) so we can recover the full
    # cycle later if needed for diagnostics.
    bot_visited: Dict[tuple, dict] = {}
    top_visited: Dict[tuple, dict] = {}
    cycled = False

    for iteration_count in range(1, MAX_FLEXURE_ITERATIONS + 1):

        # Effective depths for this iteration
        d = self.height - rec_mec

        # --- bottom tension case (positive moment on bottom face) ----------------
        (
            self._A_s_min_bot,
            self._A_s_max_bot,
            A_s_final_bot_Positive_M,  # tension req. on bottom
            A_s_comp_top,  # compression req. on top from bottom moment
            self._c_d_bot,
            self._A_s_bool_bot,
        ) = _calculate_flexural_reinforcement_ACI_318_19(self, max_M_y_bot, d, d_prima)

        # init in case no negative moment branch runs
        A_s_comp_bot = 0 * (cm**2)
        A_s_final_top_Negative_M = 0 * (cm**2)
        self._A_s_top = A_s_comp_top

        # --- top tension case (negative moment on top face) ----------------------
        if max_M_y_top < 0:
            (
                self._A_s_min_top,
                self._A_s_max_top,
                A_s_final_top_Negative_M,  # tension req. on top
                A_s_comp_bot,  # compression req. on bottom from top moment
                self._c_d_top,
                self._A_s_bool_top,
            ) = _calculate_flexural_reinforcement_ACI_318_19(
                self,
                abs(max_M_y_top.to("kN*m").magnitude) * kNm,
                self.height - d_prima,
                rec_mec,
            )

        # Governing areas on each face (tension on the face vs. opposite-face compression)
        A_req_bot = max(A_s_final_bot_Positive_M, A_s_comp_bot)
        A_req_top = max(A_s_comp_top, A_s_final_top_Negative_M)

        self._A_s_bot = A_req_bot
        self._A_s_top = A_req_top

        # --- Discrete design for each face (independent first pass) ---------------
        if A_req_bot >= 0 * (cm**2):
            A_cap_bot = self._A_s_max_bot if A_req_bot <= self._A_s_max_bot else None
            self.flexure_design_results_bot = _design_longitudinal_for_area(
                A_req_bot, A_cap_bot, self._c_mec_bot
            )

        self.flexure_design_results_top = None
        if A_req_top >= 0 * (cm**2):
            A_cap_top = self._A_s_max_top if A_req_top <= self._A_s_max_top else None
            self.flexure_design_results_top = _design_longitudinal_for_area(
                A_req_top, A_cap_top, self._c_mec_top
            )

        # --- Apply both faces (hard overwrite) -----------------------------------
        if self.flexure_design_results_bot is not None:
            self._apply_longitudinal_design_bot(self.flexure_design_results_bot)
        if self.flexure_design_results_top is not None:
            self._apply_longitudinal_design_top(self.flexure_design_results_top)
        else:
            self._clear_top_longitudinal()

        # --- Reconciliation (override only if opposite-face compression governs) -
        A_prov_bot = (
            self.flexure_design_results_bot.get("total_as", 0 * (cm**2))
            if self.flexure_design_results_bot is not None
            else 0 * (cm**2)
        )
        A_prov_top = (
            self.flexure_design_results_top.get("total_as", 0 * (cm**2))
            if self.flexure_design_results_top is not None
            else 0 * (cm**2)
        )

        # If compression from top (A_s_comp_bot) exceeds what bottom provides, re-upgrade bottom
        if A_s_comp_bot > A_prov_bot:
            A_cap_bot = self._A_s_max_bot if A_s_comp_bot <= self._A_s_max_bot else None
            self.flexure_design_results_bot = _design_longitudinal_for_area(
                A_s_comp_bot, A_cap_bot, self._c_mec_bot
            )
            if self.flexure_design_results_bot is not None:
                self._apply_longitudinal_design_bot(self.flexure_design_results_bot)
                A_prov_bot = self.flexure_design_results_bot.get("total_as", A_s_comp_bot)

        # If compression from bottom (A_s_comp_top) exceeds what top provides, re-upgrade top
        if A_s_comp_top > A_prov_top:
            if A_s_comp_top > 0 * (cm**2):
                A_cap_top = self._A_s_max_top if A_s_comp_top <= self._A_s_max_top else None
                self.flexure_design_results_top = _design_longitudinal_for_area(
                    A_s_comp_top, A_cap_top, self._c_mec_top
                )
                if self.flexure_design_results_top is not None:
                    self._apply_longitudinal_design_top(self.flexure_design_results_top)
                    A_prov_top = self.flexure_design_results_top.get("total_as", A_s_comp_top)
            else:
                self._clear_top_longitudinal()
                A_prov_top = 0 * (cm**2)

        # --- Update geometry (centroids) for next iteration ----------------------
        c_mec_calc = self.c_c + self._stirrup_d_b + self._bot_rebar_centroid
        # If there is any top steel, use its centroid; otherwise keep previous d_prima
        has_top = (self.flexure_design_results_top is not None) and (
            int(self.flexure_design_results_top.get("n_1", 0))
            + int(self.flexure_design_results_top.get("n_2", 0))
            + int(self.flexure_design_results_top.get("n_3", 0))
            + int(self.flexure_design_results_top.get("n_4", 0))
            > 0
        )
        d_prima_calc = self.c_c + self._stirrup_d_b + self._top_rebar_centroid if has_top else d_prima

        # --- Cycle detection ------------------------------------------------------
        # A single `cycled` flag covers both faces on purpose. Bottom and top are
        # coupled: bottom's layout drives `rec_mec`, which is fed back as the
        # compression-side depth of the top design (and vice-versa). If one face
        # repeats a layout it already visited, its centroid is oscillating
        # periodically, so `rec_mec`/`d_prima` are too — and that periodic input
        # forces the OTHER face into a limit cycle as well. Detecting recurrence
        # on either face is enough evidence that the whole system is in a limit
        # cycle; continuing would only waste iterations. We exit and let
        # `_select_safe_design` pick the best layout among those visited.
        if self.flexure_design_results_bot is not None:
            fp_bot = _rebar_design_fingerprint(self.flexure_design_results_bot)
            if fp_bot in bot_visited:
                cycled = True
            else:
                bot_visited[fp_bot] = dict(self.flexure_design_results_bot)
        if self.flexure_design_results_top is not None:
            fp_top = _rebar_design_fingerprint(self.flexure_design_results_top)
            if fp_top in top_visited:
                cycled = True
            else:
                top_visited[fp_top] = dict(self.flexure_design_results_top)

        # --- Convergence update ---------------------------------------------------
        Err = max(abs(c_mec_calc - rec_mec), abs(d_prima_calc - d_prima))
        rec_mec = c_mec_calc
        d_prima = d_prima_calc

        if Err < tol:
            break
        if cycled:
            break

    # --- Final capacity verification ----------------------------------------
    # Whether the loop converged, cycled, or hit MAX_ITER, the active layout
    # is NOT guaranteed to satisfy phi*Mn >= Mu. The Picard (fixed-point
    # iteration) only ensures centroid consistency, not flexural capacity. So
    # we always run a final check; if the active layout fails, we pick the
    # safest one among the visited layouts (see `_select_safe_design`). If no
    # visited layout passes either, we pick the closest one — the downstream
    # `check_flexure` will surface DCR > 1 to signal that the section is
    # insufficient.
    def _final_phi_M_n_for(face: str) -> Quantity:
        probe = Forces(M_y=(max_M_y_bot if face == "bot" else -abs(max_M_y_top)))
        _determine_nominal_moment_ACI_318_19(self, probe)
        return self._phi_M_n_bot if face == "bot" else self._phi_M_n_top

    if max_M_y_bot > 0 * kNm:
        phi_bot_active = _final_phi_M_n_for("bot")
        if phi_bot_active < max_M_y_bot:
            chosen_bot = _select_safe_design(
                self, list(bot_visited.values()), max_M_y_bot, face="bot"
            )
            self._apply_longitudinal_design_bot(chosen_bot)

    if max_M_y_top < 0 * kNm:
        M_demand_top = abs(max_M_y_top)
        phi_top_active = _final_phi_M_n_for("top")
        if phi_top_active < M_demand_top:
            chosen_top = _select_safe_design(
                self, list(top_visited.values()), M_demand_top, face="top"
            )
            self._apply_longitudinal_design_top(chosen_top)

    # --- Results table ------------------------------------------------------------
    results = {
        "Bottom_As_adopted": self._A_s_bot.to("cm**2"),
        "Bottom separation of bars": self._available_s_bot.to("mm"),
        "As_compression_adopted": self._A_s_top.to("cm**2"),
        "Top separation of bars": self._available_s_top.to("mm"),
    }
    return pd.DataFrame([results], index=[0])


def _compile_results_ACI_flexure_metric(self: "RectangularBeam", force: Forces) -> Dict[str, Any]:
    # Create dictionaries for bottom and top rows
    if self._M_u >= 0:
        result = {
            "Label": self.label,
            "Comb.": force.label,
            "Position": "Bottom",
            "As,min": round(self._A_s_min_bot.to("cm ** 2").magnitude, 2),
            "As,req top": round(self._A_s_req_top.to("cm ** 2").magnitude, 2),
            "As,req bot": round(self._A_s_req_bot.to("cm ** 2").magnitude, 2),
            "As": round(self._A_s_bot.to("cm ** 2").magnitude, 2),
            # 'c/d': self._c_d_bot,
            "Mu": round(self._M_u_bot.to("kN*m").magnitude, 2),
            "ØMn": round(self._phi_M_n_bot.to("kN*m").magnitude, 2),
            "Mu≤ØMn": self._M_u_bot <= self._phi_M_n_bot,
            "DCR": round(self._DCRb_bot, 3),
        }
    else:
        result = {
            "Label": self.label,
            "Comb.": force.label,
            "Position": "Top",
            "As,min": round(self._A_s_min_top.to("cm ** 2").magnitude, 2),
            "As,req top": round(self._A_s_req_top.to("cm ** 2").magnitude, 2),
            "As,req bot": round(self._A_s_req_bot.to("cm ** 2").magnitude, 2),
            "As": round(self._A_s_top.to("cm ** 2").magnitude, 2),
            # 'c/d': self._c_d_top,
            "Mu": round(self._M_u_top.to("kN*m").magnitude, 2),
            "ØMn": round(self._phi_M_n_top.to("kN*m").magnitude, 2),
            "Mu≤ØMn": -self._M_u_top <= self._phi_M_n_top,
            "DCR": round(self._DCRb_top, 3),
        }
    return result


##########################################################
# RESULTS
##########################################################


def _initialize_dicts_ACI_318_19_shear(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        """Initialize the dictionaries used in check and design methods."""
        self._materials_shear = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
                "Concrete density",
                "Normalweight concrete",
                "Safety factor for shear",
            ],
            "Variable": ["", "fc", "fy", "wc", "λ", "Øv"],
            "Value": [
                self.label,
                round(self.concrete.f_c.to("MPa").magnitude, 2),
                round(self.steel_bar.f_y.to("MPa").magnitude, 2),
                round(self.concrete.density.to("kg/m**3").magnitude, 1),
                self.concrete.lambda_factor,
                self.concrete.phi_v,
            ],
            "Unit": ["", "MPa", "MPa", "kg/m³", "", ""],
        }
        self._geometry_shear = {
            "Geometry": [
                "Section height",
                "Section width",
                "Clear cover",
                "Longitudinal tension rebar",
            ],
            "Variable": ["h", "b", "cc", "As"],
            "Value": [
                self.height.to("cm").magnitude,
                self.width.to("cm").magnitude,
                self.c_c.to("cm").magnitude,
                round(self._A_s_tension.to("cm**2").magnitude, 2),
            ],
            "Unit": ["cm", "cm", "cm", "cm²"],
        }
        self._forces_shear = {
            "Design forces": [
                "Axial, positive for compression",
                "Shear",
            ],
            "Variable": ["Nu", "Vu"],
            "Value": [
                round(self._N_u.to("kN").magnitude, 2),
                round(self._V_u.to("kN").magnitude, 2),
            ],
            "Unit": ["kN", "kN"],
        }
        # Min max lists
        if self._phi_V_s == 0 * kN:
            db_min = 0 * mm if self.concrete.unit_system == "metric" else 0 * inch
            self._stirrup_d_b = 0 * mm if self.concrete.unit_system == "metric" else 0 * inch
        else:
            if self.concrete.design_code == "ACI 318-19":
                db_min = 10 * mm if self.concrete.unit_system == "metric" else 3 / 8 * inch
            else:
                db_min = 6 * mm
        min_values = [
            None,
            None,
            self._A_v_min,
            db_min,
        ]  # Use None for items without a minimum constraint
        max_values = [
            self._stirrup_s_max_l,
            self._stirrup_s_max_w,
            None,
            None,
        ]  # Use None for items without a maximum constraint
        current_values = [
            self._stirrup_s_l,
            self._stirrup_s_w,
            self._A_v,
            self._stirrup_d_b,
        ]  # Current values to check
        # Generate check marks based on the range conditions
        checks = [
            "✔️" if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else "❌"
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_shear_checks_passed = all(check == "✔️" for check in checks)
        self._data_min_max_shear = {
            "Check": [
                "Stirrup spacing along length",
                "Stirrup spacing along width",
                "Minimum shear reinforcement",
                "Minimum rebar diameter",
            ],
            "Unit": ["cm", "cm", "cm²/m", "mm"],
            "Value": [
                round(self._stirrup_s_l.to("cm").magnitude, 2),
                round(self._stirrup_s_w.to("cm").magnitude, 2),
                round(self._A_v.to("cm**2/m").magnitude, 2),
                round(self._stirrup_d_b.to("mm").magnitude, 0),
            ],
            "Min.": [
                "",
                "",
                round(self._A_v_min.to("cm**2/m").magnitude, 2),
                round(db_min.to("mm").magnitude, 0),
            ],
            "Max.": [
                round(self._stirrup_s_max_l.to("cm").magnitude, 2),
                round(self._stirrup_s_max_w.to("cm").magnitude, 2),
                "",
                "",
            ],
            "Ok?": checks,
        }
        self._shear_reinforcement = {
            "Shear reinforcement strength": [
                "Number of stirrups",
                "Stirrup diameter",
                "Stirrup spacing",
                "Effective height",
                "Minimum shear reinforcing",
                "Required shear reinforcing",
                "Defined shear reinforcing",
                "Shear rebar strength",
            ],
            "Variable": ["ns", "db", "s", "d", "Av,min", "Av,req", "Av", "ØVs"],
            "Value": [
                self._stirrup_n,
                round(self._stirrup_d_b.to("mm").magnitude, 3),
                round(self._stirrup_s_l.to("cm").magnitude, 3),
                round(self._d_shear.to("cm").magnitude, 2),
                round(self._A_v_min.to("cm**2/m").magnitude, 2),
                round(self._A_v_req.to("cm**2/m").magnitude, 2),
                round(self._A_v.to("cm**2/m").magnitude, 2),
                round(self._phi_V_s.to("kN").magnitude, 2),
            ],
            "Unit": ["", "mm", "cm", "cm", "cm²/m", "cm²/m", "cm²/m", "kN"],
        }
        check_max = "✔️" if self._max_shear_ok else "❌"
        check_FU = "✔️" if self._DCRv < 1 else "❌"
        self._shear_concrete = {
            "Shear strength": [
                "Effective shear area",
                "Longitudinal reinforcement ratio",
                "Size modification factor",
                "Axial stress",
                "Concrete effective shear stress",
                "Concrete strength",
                "Maximum shear strength",
                "Total shear strength",
                "Max shear check",
                "Demand Capacity Ratio",
            ],
            "Variable": [
                "Acv",
                "ρw",
                "λs",
                "σNu",
                "kc",
                "ØVc",
                "ØVmax",
                "ØVn",
                "",
                "DCR",
            ],
            "Value": [
                round(self._A_cv.to("cm**2").magnitude, 2),
                round(self._rho_w.magnitude, 5),
                round(self._lambda_s, 3),
                round(self._sigma_Nu.to("MPa").magnitude, 2),
                round(self._k_c_min.to("MPa").magnitude, 2),
                round(self._phi_V_c.to("kN").magnitude, 2),
                round(self._phi_V_max.to("kN").magnitude, 2),
                round(self._phi_V_n.to("kN").magnitude, 2),
                check_max,
                round(self._DCRv, 2),
            ],
            "Unit": ["cm²", "", "", "MPa", "MPa", "kN", "kN", "kN", "", check_FU],
        }
        self._shear_all_checks = self._all_shear_checks_passed and (check_max == "✔️") and (check_FU == "✔️")


def _initialize_dicts_ACI_318_19_flexure(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        # Update longitudinal rebar attributes
        self._update_longitudinal_rebar_attributes()
        """Initialize the dictionaries used in check and design methods."""
        self._materials_flexure = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
            ],
            "Variable": ["", "fc", "fy"],
            "Value": [
                self.label,
                round(self.concrete.f_c.to("MPa").magnitude, 2),
                round(self.steel_bar.f_y.to("MPa").magnitude, 2),
            ],
            "Unit": ["", "MPa", "MPa"],
        }
        self._geometry_flexure = {
            "Geometry": [
                "Section height",
                "Section width",
                "Clear cover",
                "Mechanical top cover",
                "Mechanical bottom cover",
            ],
            "Variable": ["h", "b", "cc", "cm,top", "cm,bot"],
            "Value": [
                self.height.to("cm").magnitude,
                self.width.to("cm").magnitude,
                self.c_c.to("cm").magnitude,
                round(self._c_mec_top.to("cm").magnitude, 2),
                round(self._c_mec_bot.to("cm").magnitude, 2),
            ],
            "Unit": ["cm", "cm", "cm", "cm", "cm"],
        }
        self._forces_flexure = {
            "Design forces": [
                "Top max moment",
                "Bottom max moment",
            ],
            "Variable": ["Mu,top", "Mu,bot"],
            "Value": [
                round(self._M_u_top.to("kN*m").magnitude, 2),
                round(self._M_u_bot.to("kN*m").magnitude, 2),
            ],
            "Unit": ["kNm", "kNm"],
        }
        # Min max lists
        min_spacing_top: Quantity = max(
            self.settings.clear_spacing,
            self.settings.vibrator_size,
            self._d_b_max_top,
        )
        min_spacing_bot: Quantity = max(self.settings.clear_spacing, self._d_b_max_bot)
        min_values = [
            self._A_s_min_top,
            min_spacing_top,
            self._A_s_min_bot,
            min_spacing_bot,
        ]  # Use None for items without a minimum constraint
        max_values = [
            self._A_s_max_top,
            None,
            self._A_s_max_bot,
            None,
        ]  # Use None for items without a maximum constraint
        current_values = [
            self._A_s_top,
            self._available_s_top,
            self._A_s_bot,
            self._available_s_bot,
        ]  # Current values to check

        ARTICLE_STR = "9.6.1.3"

        checks = []
        for i, (curr, min_val, max_val) in enumerate(zip(current_values, min_values, max_values)):
            # --- EXCEPTION FOR DOUBLY REINFORCED SECTIONS ---
            # If doubly reinforced, ignore maximum limits for top (i=0) and bottom (i=2)
            if self._doubly_reinforced and i in (0, 2):
                # If it passes min, we give the special tag
                if min_val is None or curr >= min_val:
                    checks.append("✔️ D.R.")
                    continue
                # If it fails min, let the normal logic handle it (fall through)
            # -------------------------------------------------

            passed = (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val)
            if passed:
                checks.append("✔️")
                continue

            # Detect if fails by min or max
            failed_min = (min_val is not None) and (curr < min_val)

            if i == 0 and failed_min and getattr(self, "_A_s_bool_top", True):
                # Position 0 -> _A_s_top vs _A_s_min_top
                checks.append(ARTICLE_STR)
            elif i == 2 and failed_min and getattr(self, "_A_s_bool_bot", True):
                # Position 2 -> _A_s_bot vs _A_s_min_bot
                checks.append(ARTICLE_STR)
            else:
                # Any other failure (includes max or no flags)
                checks.append("❌")

        self._all_flexure_checks_passed = not any(check in ("❌") for check in checks)
        self._data_min_max_flexure = {
            "Check": [
                "Min/Max As rebar top",
                "Minimum spacing top",
                "Min/Max As rebar bottom",
                "Minimum spacing bottom",
            ],
            "Unit": ["cm²", "mm", "cm²", "mm"],
            "Value": [
                round(self._A_s_top.to("cm**2").magnitude, 2),
                round(self._available_s_top.to("mm").magnitude, 2),
                round(self._A_s_bot.to("cm**2").magnitude, 2),
                round(self._available_s_bot.to("mm").magnitude, 2),
            ],
            "Min.": [
                round(self._A_s_min_top.to("cm**2").magnitude, 2),
                round(min_spacing_top.to("mm").magnitude, 2),
                round(self._A_s_min_bot.to("cm**2").magnitude, 2),
                round(min_spacing_bot.to("mm").magnitude, 2),
            ],
            "Max.": [
                round(self._A_s_max_top.to("cm**2").magnitude, 2),
                "",
                round(self._A_s_max_bot.to("cm**2").magnitude, 2),
                "",
            ],
            "Ok?": checks,
        }
        check_DCR_top = "✔️" if self._DCRb_top < 1 else "❌"
        check_DCR_bot = "✔️" if self._DCRb_bot < 1 else "❌"
        self._flexure_capacity_top = {
            "Top reinforcement check": [
                "First layer bars",
                "Second layer bars",
                "Effective height",
                "Depth of equivalent strength block ratio",
                "Minimum rebar reinforcing",
                "Required rebar reinforcing top",
                "Required rebar reinforcing bottom",
                "Defined rebar reinforcing top",
                "Longitudinal reinforcement ratio",
                "Total flexural strength",
                "Demand Capacity Ratio",
            ],
            "Variable": [
                "n1+n2",
                "n3+n4",
                "d",
                "c/d",
                "As,min",
                "As,req top",
                "As,req bot",
                "As",
                "ρl",
                "ØMn",
                "DCR",
            ],
            "Value": [
                self._format_longitudinal_rebar_string(self._n1_t, self._d_b1_t, self._n2_t, self._d_b2_t),
                self._format_longitudinal_rebar_string(self._n3_t, self._d_b3_t, self._n4_t, self._d_b4_t),
                round(self._d_top.to("cm").magnitude, 2),
                self._c_d_top,
                round(self._A_s_min_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_bot.to("cm**2").magnitude, 2),
                round(self._A_s_top.to("cm**2").magnitude, 2),
                round(self._rho_l_top.magnitude, 5),
                round(self._phi_M_n_top.to("kN*m").magnitude, 2),
                round(self._DCRb_top, 2),
            ],
            "Unit": [
                "",
                "",
                "cm",
                "",
                "cm²",
                "cm²",
                "cm²",
                "cm²",
                "",
                "kNm",
                check_DCR_top,
            ],
        }
        self._flexure_capacity_bot = {
            "Bottom reinforcement check": [
                "First layer bars",
                "Second layer bars",
                "Effective height",
                "Depth of equivalent strength block ratio",
                "Minimum rebar reinforcing",
                "Required rebar reinforcing top",
                "Required rebar reinforcing bottom",
                "Defined rebar reinforcing bottom",
                "Longitudinal reinforcement ratio",
                "Total flexural strength",
                "Demand Capacity Ratio",
            ],
            "Variable": [
                "n1+n2",
                "n3+n4",
                "d",
                "c/d",
                "As,min",
                "As,req top",
                "As,req bot",
                "As",
                "ρl",
                "ØMn",
                "DCR",
            ],
            "Value": [
                self._format_longitudinal_rebar_string(self._n1_b, self._d_b1_b, self._n2_b, self._d_b2_b),
                self._format_longitudinal_rebar_string(self._n3_b, self._d_b3_b, self._n4_b, self._d_b4_b),
                round(self._d_bot.to("cm").magnitude, 2),
                self._c_d_bot,
                round(self._A_s_min_bot.to("cm**2").magnitude, 2),
                round(self._A_s_req_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_bot.to("cm**2").magnitude, 2),
                round(self._A_s_bot.to("cm**2").magnitude, 2),
                round(self._rho_l_bot.magnitude, 5),
                round(self._phi_M_n_bot.to("kN*m").magnitude, 2),
                round(self._DCRb_bot, 2),
            ],
            "Unit": [
                "",
                "",
                "cm",
                "",
                "cm²",
                "cm²",
                "cm²",
                "cm²",
                "",
                "kNm",
                check_DCR_bot,
            ],
        }
        self._flexure_all_checks = self._all_flexure_checks_passed and (check_DCR_bot == "✔️") and (check_DCR_top == "✔️")
