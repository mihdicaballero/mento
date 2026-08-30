from pint import Quantity
from typing import TYPE_CHECKING, Dict, Any, cast
import warnings
# from devtools import debug

from mento.codes.flexure_design import (
    _FaceDemand,
    _run_flexure_design,
    _select_safe_design as _select_safe_design_generic,
)
from mento.codes.aci_318_19.equations import flexure as flexure_eq
from mento.codes.shear_state import ShearCheckState, apply_shear_state, new_shear_state
from mento.codes.aci_318_19.equations import shear as shear_eq
from mento.material import Concrete_ACI_318_19
from mento.rebar import max_stirrup_spacing_ACI_318_19
from mento.units import MPa, mm, N, kN, inch, psi, cm, m, kNm, ft, kip, lbf, dimensionless
from mento.forces import Forces


if TYPE_CHECKING:
    from ..beam import RectangularBeam  # Import Beam for type checking only

# Composite units built once; `N * mm` is a pint Unit multiplication, not free.
_MOMENT_SI = N * mm
_MOMENT_US = lbf * inch


def _flexure_units(self: "RectangularBeam") -> tuple[Any, Any, Any, Any]:
    """(stress, length, area, moment) units for the beam's unit system.

    The equations layer works in floats, so every flexure call site needs the
    same four units to strip and re-apply. Picking them in one place keeps the
    pairing (MPa with mm, psi with inch) from drifting apart.
    """
    if self.concrete.is_imperial:
        return psi, inch, inch**2, _MOMENT_US
    return MPa, mm, mm**2, _MOMENT_SI


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


def _calculate_shear_reinforcement_aci(self: "RectangularBeam", st: ShearCheckState) -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        is_imperial = self.concrete.is_imperial
        length_unit = inch if is_imperial else mm
        stress_unit = psi if is_imperial else MPa
        # Shear contribution of reinforcement. A_v is an area per unit length,
        # so it converts to a plain length.
        V_s = shear_eq.shear_strength_of_reinforcement(
            self._A_v.to(length_unit).magnitude,
            st.f_yt.to(stress_unit).magnitude,
            self._d_shear.to(length_unit).magnitude,
        ) * (lbf if is_imperial else N)
        st.phi_V_s = self.concrete.phi_v * V_s  # Reduced shear contribution of reinforcement


def _calculate_effective_shear_area_aci(self: "RectangularBeam", st: ShearCheckState) -> None:
    st.A_cv = self.width * self._d_shear  # Effective shear area
    st.rho_w = st.A_s_tension.to("cm**2") / st.A_cv.to("cm**2")  # Longitudinal reinforcement ratio
    is_imperial = self.concrete.is_imperial
    length_unit = inch if is_imperial else mm
    st.lambda_s = shear_eq.size_effect_factor(self._d_shear.to(length_unit).magnitude, is_imperial=is_imperial)


def _calculate_concrete_shear_strength_aci(self: "RectangularBeam", st: ShearCheckState) -> None:
    is_imperial = self.concrete.is_imperial
    stress_unit = psi if is_imperial else MPa
    f_c = self.concrete.f_c
    # Axial stress influence
    st.sigma_Nu = (
        shear_eq.axial_stress_influence(
            st.N_u.to("lbf" if is_imperial else "N").magnitude,
            self.A_x.to("inch**2" if is_imperial else "mm**2").magnitude,
            f_c.to(stress_unit).magnitude,
        )
        * stress_unit
    )
    if isinstance(self.concrete, Concrete_ACI_318_19):
        V_cmin = 0 * (kip if is_imperial else kN)
        has_min_rebar = not (self._A_v < st.A_v_min or st.A_v_min == 0 * cm**2 / m)

        if not has_min_rebar and not is_imperial and st.A_s_tension == 0 * cm**2:
            warnings.warn(
                "Longitudinal rebar As cannot be zero if A_v is less than A_v_min.",
                UserWarning,
            )

        st.k_c_min = (
            shear_eq.concrete_shear_stress(
                f_c.to(stress_unit).magnitude,
                self.concrete.lambda_factor,
                st.rho_w.magnitude,
                st.sigma_Nu.to(stress_unit).magnitude,
                st.lambda_s,
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
        ) * st.A_cv
        st.V_c = min(V_cmax, max(V_cmin, st.k_c_min * st.A_cv))
        st.phi_V_c = self.concrete.phi_v * st.V_c


def _calculate_max_shear_capacity_aci(self: "RectangularBeam", st: ShearCheckState) -> None:
    "Formula for maximum total shear capacity (V_max)"
    if isinstance(self.concrete, Concrete_ACI_318_19):
        is_imperial = self.concrete.is_imperial
        stress_unit = psi if is_imperial else MPa
        V_max = (
            st.V_c
            + (
                shear_eq.shear_stress_capacity_increment(
                    self.concrete.f_c.to(stress_unit).magnitude,
                    self.concrete.lambda_factor,
                    is_imperial=is_imperial,
                )
                * stress_unit
            )
            * st.A_cv
        )
        st.phi_V_max = self.concrete.phi_v * V_max
        st.max_shear_ok = st.V_u < st.phi_V_max


def _calculate_A_v_min_ACI(self: "RectangularBeam", st: ShearCheckState, f_c: Quantity) -> None:
    """Calculate the minimum shear reinforcement based on unit system."""
    # 'Minimum reinforcement should be placed if the factored shear Vu
    # is greater than half the shear capacity of the concrete,
    # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
    # Rebar needed, V_u > φ_v*V_c/2 for Imperial system
    f_yt = _calculate_f_yt_aci(self)
    is_imperial = self.concrete.is_imperial
    length_unit = inch if is_imperial else mm
    stress_unit = psi if is_imperial else MPa

    st.A_v_min = (
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


def _check_minimum_reinforcement_requirement_aci(self: "RectangularBeam", st: ShearCheckState) -> None:
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
            * st.A_cv
        )

        if st.V_u < V_threshold:
            zero_A_v = 0 * inch**2 / ft if is_imperial else 0 * cm**2 / m
            st.A_v_req = zero_A_v
            st.A_v_min = zero_A_v
            st.max_shear_ok = True
        else:
            _calculate_A_v_min_ACI(self, st, self.concrete.f_c)
            st.max_shear_ok = V_threshold < st.V_u < st.phi_V_max


def _calculate_V_s_req(self: "RectangularBeam", st: ShearCheckState) -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        st.V_s_req = st.V_u - st.phi_V_c
        st.A_v_req = max(
            st.V_s_req / (self.concrete.phi_v * st.f_yt * self._d_shear),
            st.A_v_min,
        ).to("cm ** 2 / m")


def _calculate_total_shear_strength_aci(self: "RectangularBeam", st: ShearCheckState) -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        st.phi_V_n = self.concrete.phi_v * (st.V_c + self._A_v * st.f_yt * self._d_shear)
        V_d_max = min(st.phi_V_n, st.phi_V_max)
        if V_d_max.to("kN").magnitude == 0:
            # No stirrups AND no longitudinal steel on the tension face. V_c in
            # Table 22.5.5.1 scales with rho_w**(1/3), so it collapses to zero
            # and the section has no shear capacity at all -- which is what the
            # warning in _calculate_concrete_shear_strength_aci flags. Report an
            # infinite DCR rather than dividing by zero: check_shear must not
            # raise because a section is insufficient. Mirrors the guard in the
            # wall module and the phi*Mn floor in the flexure check.
            st.DCR = float("inf")
        else:
            st.DCR = abs((st.V_u.to("kN").magnitude / V_d_max.to("kN").magnitude))


def _calculate_rebar_spacing_aci(self: "RectangularBeam", st: ShearCheckState) -> None:
    n_legs_actual = self._stirrup_n * 2  # Ensure legs are even
    st.stirrup_s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs_actual - 1)
    (
        st.stirrup_s_max_l,
        st.stirrup_s_max_w,
    ) = max_stirrup_spacing_ACI_318_19(self, st.V_u - st.phi_V_c, st.A_cv)
    st.stirrup_s_w = max(st.stirrup_s_w, 0 * inch)


def _check_shear_ACI_318_19(self: "RectangularBeam", force: Forces) -> ShearCheckState:
    """Run the ACI shear check for one combination and return what it found.

    Nothing is written to the section: the result is a value. The reporting path
    copies it back through :func:`~mento.codes.shear_state.apply_shear_state`,
    which is the compatibility layer for the report tables; the values-only path
    does not, so a loop over many sections leaves every one untouched.
    """
    # No isinstance guard: the dispatcher on the element only routes ACI and
    # CIRSOC sections here, and a guard that can never fail is dead code.
    st = new_shear_state(self)

    # Demand, and the two material values the check needs. The moments belong to
    # the flexure check, so they are deliberately not touched here.
    st.N_u = force._N_x
    st.V_u = force._V_z
    st.f_yt = _calculate_f_yt_aci(self)
    st.A_s_tension = self._A_s_bot if force._M_y >= 0 * kNm else self._A_s_top

    # Minimum shear reinforcement calculation
    _calculate_A_v_min_ACI(self, st, self.concrete.f_c)
    if self._stirrup_n > 0:
        # Shear reinforcement calculations
        _calculate_shear_reinforcement_aci(self, st)
    # A section with no stirrups assigned keeps the diameter the settings
    # assume, exactly as the flexure check does, so d_shear is the same number
    # in both. Dropping the layer here used to make a flexure check report
    # differently depending on whether shear had run first. Once a design
    # decides there are no stirrups it assigns a zero diameter, and both
    # checks follow it.

    # Effective shear area and longitudinal reinforcement ratio
    _calculate_effective_shear_area_aci(self, st)

    # Check if minimum reinforcement is required
    _check_minimum_reinforcement_requirement_aci(self, st)

    # Concrete shear strength calculation
    _calculate_concrete_shear_strength_aci(self, st)

    # Maximum total shear capacity
    _calculate_max_shear_capacity_aci(self, st)

    # Calculate required shear reinforcement
    _calculate_V_s_req(self, st)

    # Total shear strength
    _calculate_total_shear_strength_aci(self, st)

    # Rebar spacing checks
    _calculate_rebar_spacing_aci(self, st)

    return st


def _design_shear_ACI_318_19(self: "RectangularBeam", force: Forces) -> None:
    """Size the shear reinforcement for one combination.

    Unlike the check, designing is *meant* to change the section — it assigns
    the stirrups. It runs the same helpers over a state and then applies it,
    so the two paths cannot drift apart.
    """
    # Set the initial variables
    _initialize_variables_ACI_318_19(self, force.M_y)
    st = new_shear_state(self)
    st.N_u = force._N_x
    st.V_u = force._V_z
    st.f_yt = self.f_yt
    st.A_s_tension = self._A_s_tension
    # Minimum shear reinforcement calculation
    _calculate_A_v_min_ACI(self, st, self.concrete.f_c)
    # Consider that the beam has minimum reinforcement
    self._A_v = st.A_v_min
    # Effective shear area and longitudinal reinforcement ratio
    _calculate_effective_shear_area_aci(self, st)
    # Concrete shear strength calculation
    _calculate_concrete_shear_strength_aci(self, st)
    # Maximum total shear capacity
    _calculate_max_shear_capacity_aci(self, st)
    # Check if minimum reinforcement is required
    _check_minimum_reinforcement_requirement_aci(self, st)
    # Calculate required shear reinforcement
    _calculate_V_s_req(self, st)
    apply_shear_state(self, st)
    # Update spacing of longitudinal reinforcement calculation
    self._update_longitudinal_rebar_attributes()

    return None


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
    stress_unit, _, _, _ = _flexure_units(self)

    return flexure_eq.max_reinforcement_ratio(
        self.concrete.f_c.to(stress_unit).magnitude,
        self.steel_bar.f_y.to(stress_unit).magnitude,
        concrete_aci.beta_1,
        concrete_aci._epsilon_c,
        self.steel_bar.epsilon_y,
    )


def _c_neutral_axis_at_ductility_limit_ACI_318_19(self: "RectangularBeam", d: Quantity) -> Quantity:
    """
    Neutral axis depth at the ACI 318-19 tension-controlled boundary.

    At the ductility limit, eps_t_min = eps_y + eps_cu. From strain
    compatibility (c/d = eps_cu / (eps_cu + eps_t)):
        c_t = 0.003 * d / (eps_y + 0.006)

    Sections with c < c_t are ductile (tension-controlled); with c > c_t are
    over-reinforced (brittle failure).
    """
    _, length_unit, _, _ = _flexure_units(self)
    return (
        flexure_eq.neutral_axis_at_ductility_limit(d.to(length_unit).magnitude, self.steel_bar.epsilon_y) * length_unit
    )


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
    stress_unit, length_unit, _, _ = _flexure_units(self)
    c_t = _c_neutral_axis_at_ductility_limit_ACI_318_19(self, d)
    return (
        flexure_eq.compression_steel_net_stress(
            d_prime.to(length_unit).magnitude,
            c_t.to(length_unit).magnitude,
            self.steel_bar.E_s.to(stress_unit).magnitude,
            self.steel_bar.f_y.to(stress_unit).magnitude,
            self.concrete.f_c.to(stress_unit).magnitude,
        )
        * stress_unit
    )


def _minimum_flexural_reinforcement_ratio_ACI_318_19(self: "RectangularBeam", M_u: Quantity) -> float:
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
        return 0.0

    stress_unit, _, _, _ = _flexure_units(self)
    return flexure_eq.min_reinforcement_ratio(
        self.concrete.f_c.to(stress_unit).magnitude,
        self.steel_bar.f_y.to(stress_unit).magnitude,
        is_imperial=self.concrete.is_imperial,
    )


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
    stress_unit, length_unit, area_unit, moment_unit = _flexure_units(self)
    f_c_mag = self.concrete.f_c.to(stress_unit).magnitude
    f_y_mag = self.steel_bar.f_y.to(stress_unit).magnitude
    b_mag = b.to(length_unit).magnitude
    d_mag = d.to(length_unit).magnitude

    R_n = flexure_eq.flexural_resistance_factor(M_u.to(moment_unit).magnitude, concrete_aci._phi_t, b_mag, d_mag)
    # Verify if the value under the square root is negative
    if flexure_eq.singly_reinforced_discriminant(R_n, f_c_mag) < 0:
        # Here we assign A_s_max so that the calculation does not break,
        # resulting in a DCR greater than 1.
        A_s_calc = A_s_max
    else:
        A_s_calc = flexure_eq.tension_steel_for_moment(R_n, f_c_mag, f_y_mag, b_mag, d_mag) * area_unit

    # Calculate the neutral axis depth based on equilibrium: 0.85 * f_c * c * beta_1 * b = A_s * f_y
    c = (
        flexure_eq.neutral_axis_depth(A_s_calc.to(area_unit).magnitude, f_y_mag, f_c_mag, b_mag, concrete_aci._beta_1)
        * length_unit
    )

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
        # Same ductility limit as A_s_max above: with ACI's fixed eps_c = 0.003,
        # eps_y + 2*eps_c is exactly the eps_y + 0.006 this branch used to spell
        # out on its own, so one function covers both.
        rho = flexure_eq.max_reinforcement_ratio(
            f_c_mag, f_y_mag, concrete_aci._beta_1, concrete_aci._epsilon_c, self.steel_bar.epsilon_y
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
    stress_unit, length_unit, area_unit, moment_unit = _flexure_units(self)
    return (
        flexure_eq.nominal_moment_singly_reinforced(
            A_s.to(area_unit).magnitude,
            self.steel_bar.f_y.to(stress_unit).magnitude,
            self.concrete.f_c.to(stress_unit).magnitude,
            self.width.to(length_unit).magnitude,
            d.to(length_unit).magnitude,
        )
        * moment_unit
    )


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
    concrete_aci = cast("Concrete_ACI_318_19", self.concrete)
    stress_unit, length_unit, area_unit, moment_unit = _flexure_units(self)

    return (
        flexure_eq.nominal_moment_doubly_reinforced(
            A_s.to(area_unit).magnitude,
            A_s_prime.to(area_unit).magnitude,
            self.steel_bar.f_y.to(stress_unit).magnitude,
            self.concrete.f_c.to(stress_unit).magnitude,
            self.width.to(length_unit).magnitude,
            d.to(length_unit).magnitude,
            d_prime.to(length_unit).magnitude,
            concrete_aci._beta_1,
            concrete_aci._epsilon_c,
            self.steel_bar._epsilon_y,
            self.steel_bar._E_s.to(stress_unit).magnitude,
        )
        * moment_unit
    )


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


def _check_flexure_ACI_318_19(self: "RectangularBeam", force: Forces) -> None:
    """
    Checks the flexural capacity of the section according to ACI 318-19 guidelines.

    This function accepts a single force and performs the flexural check of the section
    following the ACI 318-19 requirements. It initializes the design variables, computes the
    nominal moments for both top and bottom reinforcement, determines the required reinforcement
    areas, and calculates the design capacity ratios.

    Calculation only — the results are left on the beam, and turning them into a
    report table is the caller's decision (see :mod:`mento.reports.tables`).

    Parameters:
        force (Forces): The force acting on the section, which must include a single moment value.
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


##########################################################
# RESULTS
##########################################################
