import math
from pint import Quantity
from typing import TYPE_CHECKING, Dict, Any, Tuple, cast
import pandas as pd
from pandas import DataFrame

from devtools import debug


from mento.material import Concrete_EN_1992_2004
from mento.rebar import Rebar
from mento.units import MPa, mm, kNm, dimensionless, kN, inch, cm
from mento.forces import Forces


if TYPE_CHECKING:
    from ..beam import RectangularBeam  # Import Beam for type checking only


def _initialize_variables_EN_1992_2004(self: "RectangularBeam") -> None:
    """
    Initialize variables for EN 1992-2004 design code.
    """
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        self._f_ywd = self._f_ywk / self.concrete._gamma_s
        alpha_cc_shear = 1  # Take this as 1.00 for shear design and not 0.85, as in Eurocode Applied.
        self._f_cd_shear = alpha_cc_shear * self.concrete.f_ck / self.concrete.gamma_c
        # Flexure concrete resistance with a alpha_cc reduction
        self._f_cd = (
            self.concrete._alpha_cc * self.concrete.f_ck / self.concrete.gamma_c
        )


##########################################################
# SHEAR CHECK AND DESIGN
##########################################################


def _initialize_shear_variables_EN_1992_2004(
    self: "RectangularBeam", force: Forces
) -> None:
    """
    Initialize variables for EN 1992-2004 design code.
    """
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Set the initial variables
        self._N_Ed = force.N_x
        self._V_Ed_1 = (
            force.V_z
        )  # Consider the same shear at the edge of support and in d
        self._V_Ed_2 = (
            force.V_z
        )  # Consider the same shear at the edge of support and in d

        # Minimum shear reinforcement calculation
        rho_min = (
            0.08
            * math.sqrt(self.concrete.f_ck.to("MPa").magnitude)
            / (self._f_ywk)
            * MPa
        )
        self._A_v_min = rho_min * self.width * math.sin(self._alpha)

        # Consider bottom or top tension reinforcement
        self._A_s_tension = self._A_s_bot if force._M_y >= 0 * kNm else self._A_s_top
        # Compression stress, positive
        if force._M_y >= 0 * kNm:
            self._rho_l_bot = min(
                (self._A_s_tension.to("cm**2") + self._A_p.to("cm**2"))
                / (self.width.to("cm") * self._d_shear.to("cm")),
                0.02 * dimensionless,
            )
        else:
            self._rho_l_top = min(
                (self._A_s_tension + self._A_p) / (self.width * self._d_shear),
                0.02 * dimensionless,
            )

        # Shear calculation for sections without rebar
        self._k_value = min(1 + math.sqrt(200 / self._d_shear.to("mm").magnitude), 2)

        # Positive of compression
        self._sigma_cp = min(self._N_Ed / self.A_x, 0.2 * self._f_cd_shear)


def _shear_without_rebar_EN_1992_2004(self: "RectangularBeam") -> Quantity:
    self._stirrup_d_b = 0 * mm
    self._theta = 0
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Total shear capacity without rebar
        C_rdc = 0.18 / self.concrete.gamma_c
        v_min = (
            0.035
            * self._k_value ** (3 / 2)
            * math.sqrt(self.concrete.f_ck.to("MPa").magnitude)
        )
        k_1 = 0.15
        V_Rd_c_min = (
            (v_min + k_1 * self._sigma_cp.to("MPa").magnitude)
            * self.width
            * self._d_shear
            * MPa
        ).to("kN")
        rho_l = self._rho_l_bot if self._M_Ed >= 0 * kNm else self._rho_l_top
        V_Rd_c = (
            (
                C_rdc
                * self._k_value
                * (100 * rho_l * self.concrete.f_ck.to("MPa").magnitude) ** (1 / 3)
                * MPa
                + k_1 * self._sigma_cp.to("MPa")
            )
            * self.width
            * self._d_shear
        ).to("kN")
    return max(V_Rd_c_min, V_Rd_c)


def _calculate_max_shear_strength_EN_1992_2004(self: "RectangularBeam") -> None:
    """
    Calculate the maximum shear strength (V_Rd_max) and the corresponding strut angle (θ).
    """
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        alpha_cw = 1  # Non-prestressed members or members subject to tensile stress due to axial force
        v_1 = 0.6 * (
            1 - self.concrete.f_ck.to("MPa").magnitude / 250
        )  # Strength reduction factor for concrete struts
        self._z = 0.9 * self._d_shear  # Lever arm

        # The θ angle is limited between 21.8° ≤ θ ≤ 45° (1 ≤ cot(θ) ≤ 2.5)
        # Check the minimum strut angle θ = 21.8° (cot(θ) = 2.5)
        theta_min: float = math.radians(21.8)
        cot_theta_min: float = 1 / math.tan(theta_min)
        V_Rd_max_min_angle = (
            alpha_cw
            * self.width
            * self._z
            * v_1
            * self._f_cd_shear
            / (cot_theta_min + math.tan(theta_min))
        ).to("kN")

        if self._V_Ed_1 <= V_Rd_max_min_angle:
            # If within the minimum angle
            self._theta = theta_min
            self._cot_theta = cot_theta_min
            self._V_Rd_max = V_Rd_max_min_angle
            self._max_shear_ok = True
        else:
            # Check the maximum strut angle θ = 45° (cot(θ) = 1.0)
            theta_max: float = math.radians(45)
            cot_theta_max = 1 / math.tan(theta_max)
            V_Rd_max_max_angle: Quantity = (
                alpha_cw
                * self.width
                * self._z
                * v_1
                * self._f_cd_shear
                / (cot_theta_max + math.tan(theta_max))
            ).to("kN")

            if self._V_Ed_1 > V_Rd_max_max_angle:
                self._theta = theta_max
                self._cot_theta = 1 / math.tan(self._theta)
                self._V_Rd_max = V_Rd_max_max_angle
                self._max_shear_ok = False
            else:
                self._max_shear_ok = True
                # Determine the angle θ of the strut based on the shear force
                self._theta = 0.5 * math.asin((self._V_Ed_1 / V_Rd_max_max_angle))
                self._cot_theta = 1 / math.tan(self._theta)
                self._V_Rd_max = (
                    alpha_cw
                    * self.width
                    * self._z
                    * v_1
                    * self._f_cd_shear
                    / (self._cot_theta + math.tan(self._theta))
                ).to("kN")


def _calculate_required_shear_reinforcement_EN_1992_2004(
    self: "RectangularBeam",
) -> None:
    """
    Calculate the required shear reinforcement area (A_v_req).
    """
    # Calculate the required shear reinforcement area
    self._A_v_req = max(
        (
            self._V_Ed_2 / (self._z * self._f_ywd * self._cot_theta)
        ),  # Required area based on shear force
        self._A_v_min,  # Minimum required area
    )
    self._V_Rd_s = self._A_v * self._z * self._f_ywd * self._cot_theta
    self._V_s_req = self._V_Rd_s
    # Maximum shear capacity is the same as the steel capacity
    self._V_Rd = self._V_Rd_s


def _check_shear_EN_1992_2004(self: "RectangularBeam", force: Forces) -> DataFrame:
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        if self._stirrup_n == 0:
            # Set current stirrup diameter to zero
            self._stirrup_d_b = 0 * mm
            self._update_effective_heights()
            # Initialize all the code related variables
            _initialize_variables_EN_1992_2004(self)
            _initialize_shear_variables_EN_1992_2004(self, force)
            # Calculate V_Rd_c
            self._V_Rd_c = _shear_without_rebar_EN_1992_2004(self)
            # According to EN1992-1-1 §6.2.1(4) minimum shear reinforcement should nevertheless be provided
            # according to EN1992-1-1 §9.2.2. The minimum shear reinforcement may be omitted in members where
            # transverse redistribution of loads is possible (such as slabs) and members of minor importance
            # which do not contribute significantly to the overall resistance and stability of the structure.
            self._A_v_req = self._A_v_min
            # Maximum shear capacity is the same as the concrete capacity
            self._V_Rd = self._V_Rd_c
            self._V_Rd_max = self._V_Rd
            self._max_shear_ok = self._V_Ed_1 <= self._V_Rd_max

        else:
            # Initialize all the code related variables
            _initialize_variables_EN_1992_2004(self)
            _initialize_shear_variables_EN_1992_2004(self, force)
            # Shear reinforcement calculations
            d_bs = self._stirrup_d_b
            s_l = self._stirrup_s_l
            n_legs = self._stirrup_n * 2
            A_db = (d_bs**2) * math.pi / 4  # Area of one stirrup leg
            A_vs = n_legs * A_db  # Total area of stirrups
            self._A_v = A_vs / s_l  # Stirrup area per unit length

            # Calculate maximum shear strength
            _calculate_max_shear_strength_EN_1992_2004(self)

            # Calculate required shear reinforcement area
            _calculate_required_shear_reinforcement_EN_1992_2004(self)

            # Rebar spacing checks
            section_rebar = Rebar(self)
            n_legs_actual = self._stirrup_n * 2  # Ensure legs are even
            self._stirrup_s_l = self._stirrup_s_l
            self._stirrup_s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (
                n_legs_actual - 1
            )
            (
                self._stirrup_s_max_l,
                self._stirrup_s_max_w,
            ) = section_rebar.calculate_max_spacing_EN_1992_2004(self._alpha)

        self._DCRv = abs(
            (self._V_Ed_2.to("kN").magnitude / self._V_Rd.to("kN").magnitude)
        )
        # Design results
        results = {
            "Label": self.label,  # Beam label
            "Comb.": force.label,
            "Av,min": round(
                self._A_v_min.to("cm ** 2 / m").magnitude, 2
            ),  # Minimum shear reinforcement area
            "Av,req": round(
                self._A_v_req.to("cm ** 2 / m").magnitude, 2
            ),  # Required shear reinforcing area
            "Av": round(
                self._A_v.to("cm ** 2 / m").magnitude, 2
            ),  # Provided stirrup reinforcement per unit length
            "NEd": self._N_Ed.to("kN").magnitude,
            "VEd,1": self._V_Ed_1.to(
                "kN"
            ).magnitude,  # Max Vu for the design at the support
            "VEd,2": self._V_Ed_2.to(
                "kN"
            ).magnitude,  # Max Vu for the design at d from the support
            "VRd,c": round(
                self._V_Rd_c.to("kN").magnitude, 2
            ),  # Concrete contribution to shear capacity
            "VRd,s": round(
                self._V_Rd_s.to("kN").magnitude, 2
            ),  # Reinforcement contribution to shear capacity
            "VRd": round(self._V_Rd.to("kN").magnitude, 2),  # Total shear capacity
            "VRd,max": round(
                self._V_Rd_max.to("kN").magnitude, 2
            ),  # Maximum shear capacity
            "VEd,1≤VRd,max": self._max_shear_ok,  # Check if applied shear is within max shear capacity
            "VEd,2≤VRd": self._V_Ed_2.to("kN").magnitude
            <= self._V_Rd.to(
                "kN"
            ).magnitude,  # Check if applied shear is within total capacity
            "DCR": round(self._DCRv, 3),
        }
        _initialize_dicts_EN_1992_2004_shear(self)
        return pd.DataFrame([results], index=[0])
    else:
        raise ValueError("Concrete type is not compatible with EN 1992 shear check.")


def _design_shear_EN_1992_2004(self: "RectangularBeam", force: Forces) -> None:
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Initialize all the code related variables
        _initialize_variables_EN_1992_2004(self)
        # Split bottom and top moments
        _split_top_bot_moment(self, force)
        _initialize_shear_variables_EN_1992_2004(self, force)
        # Calculate maximum shear strength
        _calculate_max_shear_strength_EN_1992_2004(self)
        # Calculate required shear reinforcement area
        _calculate_required_shear_reinforcement_EN_1992_2004(self)

        return None
    else:
        raise ValueError("Concrete type is not compatible with EN 1992 shear design.")


##########################################################
# FLEXURE CHECK AND DESIGN
##########################################################


def _min_max_flexural_reinforcement_ratio_EN_1992_2004(
    self: "RectangularBeam",
) -> Tuple[float, float]:
    """
    Initialize variables for EN 1992-2004 design code.
    """
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Calculate the minimum tensile reinforcement ratio
        rho_min = max(
            0.26
            * self.concrete.f_ctm.to("MPa").magnitude
            / self.steel_bar.f_y.to("MPa").magnitude,
            0.0013,
        )
        # Set the maximum tensile reinforcement ratio
        rho_max = 0.04

        # # For positive moments (tension in the bottom), set minimum reinforcement accordingly.
        # if force._M_y > 0 * kNm:
        #     rho_min_top = 0 * dimensionless
        #     rho_min_bot = rho_min
        # else:
        #     rho_min_top = rho_min
        #     rho_min_bot = 0 * dimensionless

        # Calculate minimum and maximum bottom reinforcement areas
        # self._A_s_min_bot = rho_min_bot * self._d_bot * self._width
        # self._A_s_max_bot = rho_max * self._d_bot * self._width

    return rho_min, rho_max


def _calculate_flexural_reinforcement_EN_1992_2004(
    self: "RectangularBeam", M_Ed: Quantity, d: Quantity, d_prima: float
) -> tuple[Quantity, Quantity, Quantity, Quantity]:
    """
    Calculate the required top and bottom reinforcement areas for bending.
    """
    rho_min, rho_max = _min_max_flexural_reinforcement_ratio_EN_1992_2004(self)
    A_s_min = rho_min * d * self.width
    A_s_max = rho_max * d * self.width  # noqa: F841

    # Constants and material properties
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        lambda_ = (
            self.concrete._lambda_factor()
        )  # Factor for effective compression zone depth (EN 1992-1-1)
        eta = self.concrete._eta_factor()  # Factor for concrete strength (EN 1992-1-1)
        # Define f_yd
        f_yd = self.steel_bar.f_y / self.concrete._gamma_s

        # Relative depth of compression zone at yielding of bottom reinforcement
        concrete_en = cast("Concrete_EN_1992_2004", self.concrete)
        xi_eff_lim_1 = (self.concrete._delta - self.concrete._k_1) / self.concrete._k_2
        xi_eff_lim_2 = lambda_ * 0.45  # Ductility condition
        xi_eff_lim = min(xi_eff_lim_1, xi_eff_lim_2)
        # Compression zone depth limit
        x_eff_lim = xi_eff_lim * d
        # Limit moment for compressive reinforcement

        M_lim = eta * self._f_cd * self.width * x_eff_lim * (d - 0.5 * x_eff_lim)

        # Check if compressive reinforcement is required
        if M_Ed <= M_lim:
            # No compressive reinforcement required
            # Relative design bending moment
            K_value = (
                (M_Ed / (self.width * d**2 * eta * self._f_cd))
                .to("dimensionless")
                .magnitude
            )

            # Compression zone depth
            x_eff = d * (1 - math.sqrt(1 - 2 * K_value))

            # Area of required tensile reinforcement
            z = d - 0.5 * lambda_ * x_eff
            A_s1 = M_Ed / (z * f_yd)
            # Ensure the area meets the minimum requirement
            A_s1 = max(A_s1, A_s_min)
            # No compressive reinforcement required
            A_s2 = 0 * cm**2

        else:
            # Compressive reinforcement is required

            # Limit tensile reinforcement area
            z = d - 0.5 * lambda_ * x_eff_lim
            A_s1_lim = (M_lim / (z * f_yd)).to("cm^2")

            # Compression zone depth with plastic limit (EC2 §5.5)
            if self.concrete._f_ck <= 50 * MPa:
                x_u = (
                    d * (self.concrete._delta - self.concrete._k_1) / self.concrete._k_2
                )
            else:
                x_u = (
                    d * (self.concrete._delta - self.concrete._k_3) / self.concrete._k_4
                )

            # Compressive reinforcement strain'
            epsilon_s2 = (x_u - d_prima) / x_u * concrete_en._epsilon_cu2

            # Compressive reinforcement stress'
            f_sd = min(epsilon_s2 * self.steel_bar._E_s, f_yd)

            # Extra moment to take with top reinforcement
            delta_M = M_Ed - M_lim

            # Lever arm of internal forces for compressive reinforcement
            z_2 = d - d_prima
            # Required compressive reinforcement area
            A_s2 = (delta_M / ((z_2) * f_sd)).to("cm^2")

            # Required tensile reinforcement area
            A_s1 = max(A_s1_lim + A_s2, A_s_min)

    return A_s_min, A_s_max, A_s1, A_s2


def _simple_determine_nominal_moment_EN_1992_2004(
    self: "RectangularBeam",
    A_s: Quantity,
    d: Quantity,
    A_s_prime: Quantity,
    d_prime: Quantity,
) -> Quantity:
    # Constants and material properties
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # For reference see Jimenez Montoya 16 Ed. P 213
        f_yd = self.steel_bar._f_y / self.concrete._gamma_s
        f_cd = self.concrete._alpha_cc * self.concrete.f_ck / self.concrete.gamma_c
        b = self.width
        omega_1 = A_s * f_yd / (b * d * f_cd)  # tension rebar ratio
        omega_2 = A_s_prime / (b * d) * f_yd / f_cd  # Compression rebar ratio
        if A_s_prime.to(cm**2).magnitude >= A_s.to(cm**2).magnitude:
            z = d - d_prime  # Lever arm of internal forces
            M_Rd = A_s * f_yd * z

        elif omega_1 - omega_2 <= 0.36:
            A_tensioned = A_s - A_s_prime
            M_omega_1_omega_2 = (
                A_tensioned
                * f_yd
                * (
                    d
                    - 0.5
                    * self.concrete._lambda_factor()
                    * A_tensioned
                    * f_yd
                    / (f_cd * b)
                )
            )
            M_Rd = M_omega_1_omega_2 + A_s_prime * f_yd * (d - d_prime)
        else:
            M_Rd = 0.295 * f_cd * b * d**2 + A_s_prime * f_yd * (d - d_prime)
    return M_Rd


def _determine_nominal_moment_EN_1992_2004(
    self: "RectangularBeam", force: Forces
) -> None:
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
    [rho_min, rho_max] = _min_max_flexural_reinforcement_ratio_EN_1992_2004(self)

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
    # Determine the nominal moment for positive moments
    self._M_Rd_bot = _simple_determine_nominal_moment_EN_1992_2004(
        self, self._A_s_bot, self._d_bot, self._A_s_top, self._c_mec_top
    )
    # Determine capacity for negative moment (tension at the top)
    self._A_s_min_top = rho_min_top * self._d_top * self.width
    self._A_s_max_top = rho_max * self._d_top * self.width
    self._M_Rd_top = _simple_determine_nominal_moment_EN_1992_2004(
        self, self._A_s_top, self._d_top, self._A_s_bot, self._c_mec_bot
    )
    return None


def _split_top_bot_moment(self: "RectangularBeam", force: Forces) -> None:
    self._M_Ed = force._M_y
    if self._M_Ed > 0 * kNm:
        self._M_Ed_bot = self._M_Ed
        self._M_Ed_top = 0 * kNm
    else:
        self._M_Ed_bot = 0 * kNm
        self._M_Ed_top = self._M_Ed


def _design_flexure_EN_1992_2004(
    self: "RectangularBeam", max_M_y_bot: Quantity, max_M_y_top: Quantity
) -> Dict[str, Any]:
    # Initialize the design variables requirements using the provided force.
    _initialize_variables_EN_1992_2004(self)

    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Initialize all the code related variables
        # Initial assumptions for mechanical cover and compression depth
        rec_mec = self.c_c + self._stirrup_d_b + 1 * cm  # type: ignore
        d_prima = self.c_c + self._stirrup_d_b + 1 * cm  # type: ignore
        # Start the iterative process
        tol = 0.01 * cm  # Tolerance for convergence
        Err = 2 * tol
        iteration_count = 0
        max_iterations = 50

        while Err >= tol and iteration_count < max_iterations:
            iteration_count += 1
            # Update the effective depth for bottom tension reinforcement
            d = self.height - rec_mec
            # Calculate reinforcement for the positive moment, even if it is 0
            (
                self._A_s_min_bot,
                self._A_s_max_bot,
                A_s_final_bot_Positive_M,
                A_s_comp_top,
            ) = _calculate_flexural_reinforcement_EN_1992_2004(
                self, max_M_y_bot, d, d_prima  # type: ignore
            )
            # Initialize bottom and top reinforcement
            self._A_s_bot = A_s_final_bot_Positive_M
            self._A_s_top = A_s_comp_top

            # If there is a negative moment, calculate the top reinforcement
            if max_M_y_top < 0:
                (
                    self._A_s_min_top,
                    self._A_s_max_top,
                    A_s_final_top_Negative_M,
                    A_s_comp_bot,
                ) = _calculate_flexural_reinforcement_EN_1992_2004(
                    self, abs(max_M_y_top / kNm) * kNm, d, d_prima  # type: ignore
                )

                # Adjust reinforcement areas based on positive and negative moments
                self._A_s_bot = max(A_s_final_bot_Positive_M, A_s_comp_bot)
                self._A_s_top = max(A_s_comp_top, A_s_final_top_Negative_M)

            # Design bottom reinforcement
            section_rebar_bot = Rebar(self)
            section_rebar_bot.longitudinal_rebar_EN_1992_2004(self._A_s_bot)
            best_design_bot = section_rebar_bot.longitudinal_rebar_design

            if self._A_s_top > 0:
                section_rebar_top = Rebar(self)
                section_rebar_top.longitudinal_rebar_EN_1992_2004(self._A_s_top)
                best_design_top = section_rebar_top.longitudinal_rebar_design
            else:
                best_design_top = None

            # Set bottom rebar
            self.set_longitudinal_rebar_bot(
                best_design_bot["n_1"],
                best_design_bot["d_b1"],
                best_design_bot["n_2"],
                best_design_bot["d_b2"],
                best_design_bot["n_3"],
                best_design_bot["d_b3"],
                best_design_bot["n_4"],
                best_design_bot["d_b4"],
            )
            if best_design_top is not None:
                self.set_longitudinal_rebar_top(
                    best_design_top["n_1"],
                    best_design_top["d_b1"],
                    best_design_top["n_2"],
                    best_design_top["d_b2"],
                    best_design_top["n_3"],
                    best_design_top["d_b3"],
                    best_design_top["n_4"],
                    best_design_top["d_b4"],
                )
            else:
                self.set_longitudinal_rebar_top(
                    0, 0 * mm, 0, 0 * mm, 0, 0 * mm, 0, 0 * mm
                )

            # Recalculate rebar centroid
            c_mec_calc = self.c_c + self._stirrup_d_b + self._bot_rebar_centroid

            if best_design_top is not None and self._A_s_top > 0:
                d_prima_calc = self.c_c + self._stirrup_d_b + self._top_rebar_centroid
            else:
                d_prima_calc = d_prima

            # Update error for iteration
            Err = max(abs(c_mec_calc - rec_mec), abs(d_prima_calc - d_prima))
            rec_mec = c_mec_calc
            d_prima = d_prima_calc

        if iteration_count >= max_iterations:
            # Podés loguear un warning si querés:
            Warning("Flexure EN 1992 did not converge within max_iterations")
            pass

        # Return results as a DataFrame
        results = {
            "Bottom_As_adopted": self._A_s_bot.to("inch**2"),
            "Bottom separation of bars": self._available_s_bot.to("inch"),
            "As_compression_adopted": self._A_s_top.to("inch**2"),
            "Top separation of bars": self._available_s_top.to("inch"),
        }
        return pd.DataFrame([results], index=[0])
    else:
        raise ValueError("Concrete type is not compatible with EN 1992 flexure design.")


def _check_flexure_EN_1992_2004(self: "RectangularBeam", force: Forces) -> pd.DataFrame:
    """ """
    # Initialize the design variables requirements using the provided force.
    _initialize_variables_EN_1992_2004(self)

    # Split bottom and top moments
    _split_top_bot_moment(self, force)

    # Calculate the nominal moments for both top and bottom reinforcement.

    _determine_nominal_moment_EN_1992_2004(self, force)
    if self._M_Ed >= 0:
        # For positive moments, calculate the reinforcement requirements for the bottom tension side.
        (
            self._A_s_min_bot,
            self._A_s_max_bot,
            self._A_s_req_bot,
            self._A_s_req_top,
        ) = _calculate_flexural_reinforcement_EN_1992_2004(
            self, self._M_Ed_bot, self._d_bot, self._c_mec_top
        )
        self._c_d_top = 0
        # Calculate the design capacity ratio for the bottom side.
        self._DCRb_bot = round(
            self._M_Ed_bot.to("kN*m").magnitude / self._M_Rd_bot.to("kN*m").magnitude,
            3,
        )
        self._DCRb_top = 0
    else:
        # For negative moments, calculate the reinforcement requirements for the top tension side.
        (
            self._A_s_min_top,
            self._A_s_max_top,
            self._A_s_req_top,
            self._A_s_req_bot,
        ) = _calculate_flexural_reinforcement_EN_1992_2004(
            self, abs(self._M_Ed_top / kNm) * kNm, self._d_top, self._c_mec_bot
        )
        self._c_d_bot = 0
        # Calculate the design capacity ratio for the top side.
        self._DCRb_top = round(
            -self._M_Ed_top.to("kN*m").magnitude / self._M_Rd_top.to("kN*m").magnitude,
            3,
        )
        self._DCRb_bot = 0

    # Determine the maximum detailing cover dimensions for top and bottom.
    self._d_b_max_top = max(self._d_b1_t, self._d_b2_t, self._d_b3_t, self._d_b4_t)
    self._d_b_max_bot = max(self._d_b1_b, self._d_b2_b, self._d_b3_b, self._d_b4_b)

    # Calculate the longitudinal reinforcement ratios for both sides.
    self._rho_l_bot = self._A_s_bot / (self._d_bot * self.width)
    self._rho_l_top = self._A_s_bot / (self._d_top * self.width)

    # Compile the design results into a dictionary.
    results = _compile_results_EN_1992_2004_flexure_metric(self, force)

    # Initialize any additional dictionaries required for ACI 318-19 flexural checks.
    _initialize_dicts_EN_1992_2004_flexure(self)

    # Return the results as a Pandas DataFrame.
    return pd.DataFrame([results], index=[0])


def _compile_results_EN_1992_2004_flexure_metric(
    self: "RectangularBeam", force: Forces
) -> Dict[str, Any]:
    # Create dictionaries for bottom and top rows
    if self._M_Ed >= 0:
        result = {
            "Label": self.label,
            "Comb.": force.label,
            "Position": "Bottom",
            "As,min": round(self._A_s_min_bot.to("cm ** 2").magnitude, 2),
            "As,req top": round(self._A_s_req_top.to("cm ** 2").magnitude, 2),
            "As,req bot": round(self._A_s_req_bot.to("cm ** 2").magnitude, 2),
            "As": round(self._A_s_bot.to("cm ** 2").magnitude, 2),
            # 'c/d': self._c_d_bot,
            "MEd": round(self._M_Ed_bot.to("kN*m").magnitude, 2),
            "MRd": round(self._M_Rd_bot.to("kN*m").magnitude, 2),
            "MEd≤MRd": self._M_Ed_bot <= self._M_Rd_bot,
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
            "MEd": round(self._M_Ed_top.to("kN*m").magnitude, 2),
            "MRd": round(self._M_Rd_top.to("kN*m").magnitude, 2),
            "MEd≤MRd": self._M_Ed_top <= self._M_Rd_top,
            "DCR": round(self._DCRb_top, 3),
        }
    return result


##########################################################
# RESULTS
##########################################################


def _initialize_dicts_EN_1992_2004_shear(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        """Initialize the dictionaries used in check and design methods."""
        self._materials_shear = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
                "Safety factor for concrete",
                "Safety factor for steel",
                "Coefficient for long term effects and loading effects",
            ],
            "Variable": ["", "fck", "fywk", "γc", "γs", "αcc"],
            "Value": [
                self.label,
                round(self.concrete.f_ck.to("MPa").magnitude, 2),
                round(self.steel_bar.f_y.to("MPa").magnitude, 2),
                self.concrete.gamma_c,
                self.concrete._gamma_s,
                self.concrete.alpha_cc,
            ],
            "Unit": ["", "MPa", "MPa", "", "", ""],
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
            "Variable": ["NEd", "VEd,2"],
            "Value": [
                round(self._N_Ed.to("kN").magnitude, 2),
                round(self._V_Ed_2.to("kN").magnitude, 2),
            ],
            "Unit": ["kN", "kN"],
        }
        # Min max lists
        if self._V_Rd_s == 0 * kN:
            self._stirrup_d_b = (
                0 * mm if self.concrete.unit_system == "metric" else 0 * inch
            )
        # Min max lists
        min_values = [
            None,
            None,
            self._A_v_min,
        ]  # Use None for items without a minimum constraint
        max_values = [
            self._stirrup_s_max_l,
            self._stirrup_s_max_w,
            None,
        ]  # Use None for items without a maximum constraint
        current_values = [
            self._stirrup_s_l,
            self._stirrup_s_w,
            self._A_v,
        ]  # Current values to check

        # Generate check marks based on the range conditions
        checks = [
            "✔️"
            if (min_val is None or curr >= min_val)
            and (max_val is None or curr <= max_val)
            else "❌"
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_shear_checks_passed = all(check == "✔️" for check in checks)
        self._data_min_max_shear = {
            "Check": [
                "Stirrup spacing along length",
                "Stirrup spacing along width",
                "Minimum shear reinforcement",
            ],
            "Unit": ["cm", "cm", "cm²/m"],
            "Value": [
                round(self._stirrup_s_l.to("cm").magnitude, 2),
                round(self._stirrup_s_w.to("cm").magnitude, 2),
                round(self._A_v.to("cm**2/m").magnitude, 2),
            ],
            "Min.": ["", "", round(self._A_v_min.to("cm**2/m").magnitude, 2)],
            "Max.": [
                round(self._stirrup_s_max_l.to("cm").magnitude, 2),
                round(self._stirrup_s_max_w.to("cm").magnitude, 2),
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
            "Variable": ["ns", "db", "s", "d", "Asw,min", "Asw,req", "Asw", "VRd,s"],
            "Value": [
                self._stirrup_n,
                self._stirrup_d_b.to("mm").magnitude,
                self._stirrup_s_l.to("cm").magnitude,
                round(self._d_shear.to("cm").magnitude, 2),
                round(self._A_v_min.to("cm**2/m").magnitude, 2),
                round(self._A_v_req.to("cm**2/m").magnitude, 2),
                round(self._A_v.to("cm**2/m").magnitude, 2),
                round(self._V_Rd_s.to("kN").magnitude, 2),
            ],
            "Unit": ["", "mm", "cm", "cm", "cm²/m", "cm²/m", "cm²/m", "kN"],
        }
        check_max = "✔️" if self._max_shear_ok else "❌"
        check_DCR = "✔️" if self._DCRv < 1 else "❌"
        rho_l = self._rho_l_bot if self._M_Ed >= 0 * kNm else self._rho_l_top
        self._shear_concrete = {
            "Shear strength": [
                "Longitudinal reinforcement ratio",
                "k value",
                "Axial stress",
                "Concrete strut angle",
                "Concrete strength",
                "Maximum shear strength",
                "Total shear strength",
                "Max shear check",
                "Demand Capacity Ratio",
            ],
            "Variable": ["ρl", "k", "σcd", "Θ", "VRd,c", "VRd,max", "VRd", "", "DCR"],
            "Value": [
                round(rho_l.magnitude, 4),
                round(self._k_value, 2),
                round(self._sigma_cp.to("MPa").magnitude, 2),
                round(math.degrees(self._theta), 1),
                round(self._V_Rd_c.to("kN").magnitude, 2),
                round(self._V_Rd_max.to("kN").magnitude, 2),
                round(self._V_Rd.to("kN").magnitude, 2),
                check_max,
                round(self._DCRv, 3),
            ],
            "Unit": ["", "", "MPa", "deg", "kN", "kN", "kN", "", check_DCR],
        }


def _initialize_dicts_EN_1992_2004_flexure(self: "RectangularBeam") -> None:
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        """Initialize the dictionaries used in check and design methods."""
        self._materials_flexure = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
            ],
            "Variable": ["", "fck", "fyk"],
            "Value": [
                self.label,
                round(self.concrete.f_ck.to("MPa").magnitude, 2),
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
            "Variable": ["MEd,top", "MEd,bot"],
            "Value": [
                round(self._M_Ed_top.to("kN*m").magnitude, 2),
                round(self._M_Ed_bot.to("kN*m").magnitude, 2),
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

        # Generate check marks based on the range conditions
        checks = [
            "✔️"
            if (min_val is None or curr >= min_val)
            and (max_val is None or curr <= max_val)
            else "❌"
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_flexure_checks_passed = all(check == "✔️" for check in checks)
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
                "MRd",
                "DCR",
            ],
            "Value": [
                self._format_longitudinal_rebar_string(
                    self._n1_t, self._d_b1_t, self._n2_t, self._d_b2_t
                ),
                self._format_longitudinal_rebar_string(
                    self._n3_t, self._d_b3_t, self._n4_t, self._d_b4_t
                ),
                round(self._d_top.to("cm").magnitude, 2),
                self._c_d_top,
                round(self._A_s_min_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_bot.to("cm**2").magnitude, 2),
                round(self._A_s_top.to("cm**2").magnitude, 2),
                round(self._rho_l_top.magnitude, 5),
                round(self._M_Rd_top.to("kN*m").magnitude, 2),
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
                "Required rebar reinforcing bottom",
                "Required rebar reinforcing top",
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
                "MRd",
                "DCR",
            ],
            "Value": [
                self._format_longitudinal_rebar_string(
                    self._n1_b, self._d_b1_b, self._n2_b, self._d_b2_b
                ),
                self._format_longitudinal_rebar_string(
                    self._n3_b, self._d_b3_b, self._n4_b, self._d_b4_b
                ),
                round(self._d_bot.to("cm").magnitude, 2),
                self._c_d_bot,
                round(self._A_s_min_bot.to("cm**2").magnitude, 2),
                round(self._A_s_req_top.to("cm**2").magnitude, 2),
                round(self._A_s_req_bot.to("cm**2").magnitude, 2),
                round(self._A_s_bot.to("cm**2").magnitude, 2),
                round(self._rho_l_bot.magnitude, 5),
                round(self._M_Rd_bot.to("kN*m").magnitude, 2),
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
            self._all_flexure_checks_passed
            and (check_DCR_bot == "✔️")
            and (check_DCR_top == "✔️")
        )
