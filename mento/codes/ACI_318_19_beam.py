from pint import Quantity
import pandas as pd
import numpy as np
from typing import TYPE_CHECKING, Dict, Any, cast
import warnings
# from devtools import debug

from mento.codes.flexure_design import (
    _FaceDemand,
    _run_flexure_design,
    _select_safe_design as _select_safe_design_generic,
)
from mento.codes.aci_318_19.equations import shear as shear_eq
from mento.material import Concrete_ACI_318_19
from mento.rebar import Rebar
from mento.units import MPa, mm, N, kN, inch, psi, cm, m, kNm, ft, kip, lbf, dimensionless
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
    if isinstance(self.concrete, Concrete_ACI_318_19):
        is_imperial = self.concrete.is_imperial
        length_unit = inch if is_imperial else mm
        stress_unit = psi if is_imperial else MPa
        # Shear contribution of reinforcement. A_v is an area per unit length,
        # so it converts to a plain length.
        V_s = shear_eq.shear_strength_of_reinforcement(
            self._A_v.to(length_unit).magnitude,
            self.f_yt.to(stress_unit).magnitude,
            self._d_shear.to(length_unit).magnitude,
        ) * (lbf if is_imperial else N)
        self._phi_V_s = self.concrete.phi_v * V_s  # Reduced shear contribution of reinforcement


def _calculate_effective_shear_area_aci(self: "RectangularBeam") -> None:
    self._A_cv = self.width * self._d_shear  # Effective shear area
    self._rho_w = self._A_s_tension.to("cm**2") / self._A_cv.to("cm**2")  # Longitudinal reinforcement ratio
    is_imperial = self.concrete.is_imperial
    length_unit = inch if is_imperial else mm
    self._lambda_s = shear_eq.size_effect_factor(self._d_shear.to(length_unit).magnitude, is_imperial=is_imperial)


def _calculate_concrete_shear_strength_aci(self: "RectangularBeam") -> None:
    is_imperial = self.concrete.is_imperial
    stress_unit = psi if is_imperial else MPa
    f_c = self.concrete.f_c
    # Axial stress influence
    self._sigma_Nu = (
        shear_eq.axial_stress_influence(
            self._N_u.to("lbf" if is_imperial else "N").magnitude,
            self.A_x.to("inch**2" if is_imperial else "mm**2").magnitude,
            f_c.to(stress_unit).magnitude,
        )
        * stress_unit
    )
    if isinstance(self.concrete, Concrete_ACI_318_19):
        V_cmin = 0 * (kip if is_imperial else kN)
        has_min_rebar = not (self._A_v < self._A_v_min or self._A_v_min == 0 * cm**2 / m)

        if not has_min_rebar and not is_imperial and self._A_s_tension == 0 * cm**2:
            warnings.warn(
                "Longitudinal rebar As cannot be zero if A_v is less than A_v_min.",
                UserWarning,
            )

        self._k_c_min = (
            shear_eq.concrete_shear_stress(
                f_c.to(stress_unit).magnitude,
                self.concrete.lambda_factor,
                self._rho_w.magnitude,
                self._sigma_Nu.to(stress_unit).magnitude,
                self._lambda_s,
                has_min_rebar=has_min_rebar,
                is_imperial=is_imperial,
            )
            * stress_unit
        )
        # Maximum concrete shear strength
        V_cmax = (
            shear_eq.max_concrete_shear_stress(
                f_c.to(stress_unit).magnitude, self.concrete.lambda_factor, is_imperial=is_imperial
            )
            * stress_unit
        ) * self._A_cv
        self.V_c = min(V_cmax, max(V_cmin, self._k_c_min * self._A_cv))
        self._phi_V_c = self.concrete.phi_v * self.V_c


def _calculate_max_shear_capacity_aci(self: "RectangularBeam") -> None:
    "Formula for maximum total shear capacity (V_max)"
    if isinstance(self.concrete, Concrete_ACI_318_19):
        is_imperial = self.concrete.is_imperial
        stress_unit = psi if is_imperial else MPa
        V_max = (
            self.V_c
            + (
                shear_eq.shear_stress_capacity_increment(
                    self.concrete.f_c.to(stress_unit).magnitude,
                    self.concrete.lambda_factor,
                    is_imperial=is_imperial,
                )
                * stress_unit
            )
            * self._A_cv
        )
        self._phi_V_max = self.concrete.phi_v * V_max
        self._max_shear_ok = self._V_u < self._phi_V_max


def _calculate_A_v_min_ACI(self: "RectangularBeam", f_c: Quantity) -> None:
    """Calculate the minimum shear reinforcement based on unit system."""
    # 'Minimum reinforcement should be placed if the factored shear Vu
    # is greater than half the shear capacity of the concrete,
    # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
    # Rebar needed, V_u > φ_v*V_c/2 for Imperial system
    f_yt = _calculate_f_yt_aci(self)
    is_imperial = self.concrete.is_imperial
    length_unit = inch if is_imperial else mm
    stress_unit = psi if is_imperial else MPa

    self._A_v_min = (
        shear_eq.min_shear_reinforcement_ratio(
            f_c.to(stress_unit).magnitude,
            f_yt.to(stress_unit).magnitude,
            self.width.to(length_unit).magnitude,
            is_imperial=is_imperial,
        )
        * length_unit
    )


def _calculate_f_yt_aci(self: "RectangularBeam") -> Quantity:
    """Determine the yield strength of steel based on unit system."""
    is_imperial = self.concrete.is_imperial
    stress_unit = psi if is_imperial else MPa
    return (
        shear_eq.max_yield_strength_for_shear(self.steel_bar.f_y.to(stress_unit).magnitude, is_imperial=is_imperial)
        * stress_unit
    )


def _check_minimum_reinforcement_requirement_aci(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        is_imperial = self.concrete.is_imperial
        stress_unit = psi if is_imperial else MPa
        # Demand below which ACI 318-19 §9.6.3.1 waives shear reinforcement.
        V_threshold = (
            self.concrete.phi_v
            * (
                shear_eq.min_shear_reinforcement_threshold_stress(
                    self.concrete.f_c.to(stress_unit).magnitude,
                    self.concrete.lambda_factor,
                    is_imperial=is_imperial,
                )
                * stress_unit
            )
            * self._A_cv
        )

        if self._V_u < V_threshold:
            zero_A_v = 0 * inch**2 / ft if is_imperial else 0 * cm**2 / m
            self._A_v_req = zero_A_v
            self._A_v_min = zero_A_v
            self._max_shear_ok = True
        else:
            _calculate_A_v_min_ACI(self, self.concrete.f_c)
            self._max_shear_ok = V_threshold < self._V_u < self._phi_V_max


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
        if V_d_max.to("kN").magnitude == 0:
            # No stirrups AND no longitudinal steel on the tension face. V_c in
            # Table 22.5.5.1 scales with rho_w**(1/3), so it collapses to zero
            # and the section has no shear capacity at all -- which is what the
            # warning in _calculate_concrete_shear_strength_aci flags. Report an
            # infinite DCR rather than dividing by zero: check_shear must not
            # raise because a section is insufficient. Mirrors the guard in the
            # wall module and the phi*Mn floor in the flexure check.
            self._DCRv = float("inf")
        else:
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


def _c_neutral_axis_at_ductility_limit_ACI_318_19(self: "RectangularBeam", d: Quantity) -> Quantity:
    """
    Neutral axis depth at the ACI 318-19 tension-controlled boundary.

    At the ductility limit, eps_t_min = eps_y + eps_cu. From strain
    compatibility (c/d = eps_cu / (eps_cu + eps_t)):
        c_t = 0.003 * d / (eps_y + 0.006)

    Sections with c < c_t are ductile (tension-controlled); with c > c_t are
    over-reinforced (brittle failure).
    """
    return 0.003 * d / (self.steel_bar.epsilon_y + 0.006)


def _f_s_prime_net_at_ductility_limit_ACI_318_19(self: "RectangularBeam", d: Quantity, d_prime: Quantity) -> Quantity:
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

    # 1.8‰ of the gross section (custom geometric minimum rule), same for beams and slabs
    A_s_geo_min = (1.8 / (1000)) * self.width * self.height

    if M_u == 0 * kNm:
        # Case 0:
        # No flexural demand (e.g. a shear-only load combination). ACI does not
        # require flexural minimum steel here, so rho_min -- and therefore
        # A_s_min -- is zero, but leaving A_s = 0 is not a buildable layout: the
        # section still needs detailing steel, and rho_w = 0 collapses V_c to
        # zero in the shear provisions (Table 22.5.5.1). Adopt the geometric
        # minimum.
        A_s_final = A_s_geo_min
    elif A_s_calc >= A_s_min:
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
        M_n = 0.85 * f_c * a_assumed * b * (d - a_assumed / 2) + A_s_prime * (f_y - 0.85 * f_c) * (d - d_prime)
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
        M_n_positive = _determine_nominal_moment_simple_reinf_ACI_318_19(self, self._A_s_bot, self._d_bot)
    else:
        # Compression steel contribution at the ductility limit determines how
        # much the tension-steel cap can be extended:
        #     A_s_max_total = A_s_max_bot + A_s_top * f_s'_net / f_y
        f_s_prima_net = _f_s_prime_net_at_ductility_limit_ACI_318_19(self, self._d_bot, self._c_mec_top)
        A_s_max_total = self._A_s_max_bot + self._A_s_top * f_s_prima_net / self.steel_bar.f_y
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
        M_n_negative = _determine_nominal_moment_simple_reinf_ACI_318_19(self, self._A_s_top, self._d_top)
    else:
        f_s_prima_net = _f_s_prime_net_at_ductility_limit_ACI_318_19(self, self._d_top, self._c_mec_bot)
        A_s_max_total = self._A_s_max_top + self._A_s_bot * f_s_prima_net / self.steel_bar.f_y
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


# ─────────────────────────────────────────────────────────────────────────────
# Flexural design — ACI 318-19 hooks for the shared driver
# ─────────────────────────────────────────────────────────────────────────────
# The design *strategy* -- Picard loop on the mechanical covers, discrete rebar
# selection, face reconciliation, cycle detection and final capacity
# verification -- lives in `mento.codes.flexure_design` and is shared with
# EN 1992-2004. Only the two hooks below are ACI-specific.


def _flexure_capacity_ACI_318_19(self: "RectangularBeam", face: str, M_demand: Quantity) -> Quantity:
    """phi*Mn of the layout currently applied to the section, on `face`.

    Recomputed rather than read from cached state so that the centroid of the
    layout just applied is taken into account.
    """
    M_abs = abs(M_demand)
    probe_force = Forces(M_y=(M_abs if face == "bot" else -M_abs))
    _determine_nominal_moment_ACI_318_19(self, probe_force)
    return self._phi_M_n_bot if face == "bot" else self._phi_M_n_top


def _required_areas_ACI_318_19(
    self: "RectangularBeam", face: str, M: Quantity, d: Quantity, d_prime: Quantity
) -> _FaceDemand:
    """Steel required by ACI 318-19 on `face` for the moment `M`.

    The code-specific extras that do not belong in `_FaceDemand` (the c/d ratio
    and the 4/3-rule flag) are stored on the beam for the results tables.
    """
    (
        A_s_min,
        A_s_max,
        A_s_tension,
        A_s_compression,
        c_d,
        A_s_bool,
    ) = _calculate_flexural_reinforcement_ACI_318_19(self, M, d, d_prime)
    if face == "bot":
        self._A_s_min_bot, self._A_s_max_bot = A_s_min, A_s_max
        self._c_d_bot, self._A_s_bool_bot = c_d, A_s_bool
    else:
        self._A_s_min_top, self._A_s_max_top = A_s_min, A_s_max
        self._c_d_top, self._A_s_bool_top = c_d, A_s_bool
    return _FaceDemand(A_s_min, A_s_max, A_s_tension, A_s_compression)


def _select_safe_design(
    self: "RectangularBeam",
    candidate_designs: list,
    M_demand: Quantity,
    face: str,
) -> Dict[str, Any]:
    """ACI-flavoured safe-layout selection — same selection rules, with
    phi*Mn as the capacity measure."""

    def _capacity(f: str, M: Quantity) -> Quantity:
        return _flexure_capacity_ACI_318_19(self, f, M)

    return _select_safe_design_generic(self, candidate_designs, M_demand, face, _capacity)


def _design_flexure_ACI_318_19(self: "RectangularBeam", max_M_y_bot: Quantity, max_M_y_top: Quantity) -> None:
    """Design the longitudinal reinforcement of a beam per ACI 318-19.

    Thin wrapper: everything that is not an ACI equation lives in
    ``mento.codes.flexure_design``.
    """

    def _required(face: str, M: Quantity, d: Quantity, d_prime: Quantity) -> _FaceDemand:
        return _required_areas_ACI_318_19(self, face, M, d, d_prime)

    def _capacity(face: str, M: Quantity) -> Quantity:
        return _flexure_capacity_ACI_318_19(self, face, M)

    _run_flexure_design(self, max_M_y_bot, max_M_y_top, _required, _capacity)


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
            "✅" if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else "❌"
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_shear_checks_passed = all(check == "✅" for check in checks)
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
        check_max = "✅" if self._max_shear_ok else "❌"
        check_FU = "✅" if self._DCRv < 1 else "❌"
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
        self._shear_all_checks = self._all_shear_checks_passed and (check_max == "✅") and (check_FU == "✅")


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
                    checks.append("✅ D.R.")
                    continue
                # If it fails min, let the normal logic handle it (fall through)
            # -------------------------------------------------

            passed = (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val)
            if passed:
                checks.append("✅")
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
        check_DCR_top = "✅" if self._DCRb_top < 1 else "❌"
        check_DCR_bot = "✅" if self._DCRb_bot < 1 else "❌"
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
        self._flexure_all_checks = (
            self._all_flexure_checks_passed and (check_DCR_bot == "✅") and (check_DCR_top == "✅")
        )
