import math
from pint import Quantity
from typing import TYPE_CHECKING, Dict, Any, Tuple, cast
import pandas as pd
from pandas import DataFrame

# from devtools import debug


from mento.codes.en_1992_2004.equations import flexure as flexure_eq
from mento.codes.en_1992_2004.equations import shear as shear_eq
from mento.codes.flexure_design import _FaceDemand, _run_flexure_design
from mento.material import Concrete_EN_1992_2004
from mento.rebar import max_stirrup_spacing_EN_1992_2004
from mento.units import MPa, mm, kNm, dimensionless, kN, N, inch
from mento.forces import Forces


if TYPE_CHECKING:
    from ..beam import RectangularBeam  # Import Beam for type checking only

# Prebuilt unit objects for the boundary conversions of ADR-0005. The equation
# layer works in N, mm, mm², MPa and N·mm; these are the only units crossing
# into and out of it, and building them once keeps the conversion off the hot
# path (`.to(_mm2)`, never `.to("mm**2")`).
_mm2 = mm**2
_Nmm = N * mm


def _initialize_variables_EN_1992_2004(self: "RectangularBeam") -> None:
    """
    Initialize variables for EN 1992-2004 design code.
    """
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        self._f_ywd = self._f_ywk / self.concrete._gamma_s
        alpha_cc_shear = 1  # Take this as 1.00 for shear design and not 0.85, as in Eurocode Applied.
        self._f_cd_shear = alpha_cc_shear * self.concrete.f_ck / self.concrete.gamma_c
        # Flexure concrete resistance with a alpha_cc reduction
        self._f_cd = self.concrete._alpha_cc * self.concrete.f_ck / self.concrete.gamma_c


##########################################################
# SHEAR CHECK AND DESIGN
##########################################################


def _initialize_shear_variables_EN_1992_2004(self: "RectangularBeam", force: Forces) -> None:
    """
    Initialize variables for EN 1992-2004 design code.
    """
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Set the initial variables
        self._N_Ed = force.N_x
        self._V_Ed_1 = force.V_z  # Consider the same shear at the edge of support and in d
        self._V_Ed_2 = force.V_z  # Consider the same shear at the edge of support and in d

        # Minimum shear reinforcement calculation. Eq. (9.5N) gives the ratio;
        # §9.2.2(5) defines it as A_sw/(s*b_w*sin(alpha)), hence the geometry here.
        rho_min = shear_eq.min_shear_reinforcement_ratio(
            self.concrete.f_ck.to(MPa).magnitude,
            self._f_ywk.to(MPa).magnitude,
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
        self._k_value = shear_eq.size_effect_factor(self._d_shear.to(mm).magnitude)

        # Positive of compression
        self._sigma_cp = (
            shear_eq.axial_stress(
                self._N_Ed.to(N).magnitude,
                self.A_x.to(_mm2).magnitude,
                self._f_cd_shear.to(MPa).magnitude,
            )
            * MPa
        )


def _shear_without_rebar_EN_1992_2004(self: "RectangularBeam") -> Quantity:
    self._stirrup_d_b = 0 * mm
    self._theta = 0
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Total shear capacity without rebar: Eq. (6.2.a), floored by Eq. (6.2.b).
        f_ck = self.concrete.f_ck.to(MPa).magnitude
        b_w = self.width.to(mm).magnitude
        d = self._d_shear.to(mm).magnitude
        sigma_cp = self._sigma_cp.to(MPa).magnitude
        rho_l = self._rho_l_bot if self._M_Ed >= 0 * kNm else self._rho_l_top
        V_Rd_c_min = (
            shear_eq.min_shear_resistance_without_reinforcement(f_ck, self._k_value, sigma_cp, b_w, d) * N
        ).to(kN)
        V_Rd_c = (
            shear_eq.shear_resistance_without_reinforcement(
                f_ck,
                self.concrete.gamma_c,
                self._k_value,
                rho_l.to(dimensionless).magnitude,
                sigma_cp,
                b_w,
                d,
            )
            * N
        ).to(kN)
    return max(V_Rd_c_min, V_Rd_c)


def _calculate_max_shear_strength_EN_1992_2004(self: "RectangularBeam") -> None:
    """
    Calculate the maximum shear strength (V_Rd_max) and the corresponding strut angle (θ).
    """
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        alpha_cw = 1.0  # Non-prestressed members or members subject to tensile stress due to axial force
        nu_1 = shear_eq.strut_strength_reduction_factor(self.concrete.f_ck.to(MPa).magnitude)
        b_w = self.width.to(mm).magnitude
        f_cd = self._f_cd_shear.to(MPa).magnitude
        z = shear_eq.lever_arm(self._d_shear.to(mm).magnitude)
        self._z = z * mm  # Lever arm

        # The θ angle is limited between 21.8° ≤ θ ≤ 45° (1 ≤ cot(θ) ≤ 2.5)
        # Check the minimum strut angle θ = 21.8° (cot(θ) = 2.5)
        theta_min: float = math.radians(21.8)
        cot_theta_min: float = 1 / math.tan(theta_min)
        V_Rd_max_min_angle = (shear_eq.max_shear_resistance(alpha_cw, b_w, z, nu_1, f_cd, theta_min) * N).to(kN)

        if self._V_Ed_1 <= V_Rd_max_min_angle:
            # If within the minimum angle
            self._theta = theta_min
            self._cot_theta = cot_theta_min
            self._V_Rd_max = V_Rd_max_min_angle
            self._max_shear_ok = True
        else:
            # Check the maximum strut angle θ = 45° (cot(θ) = 1.0)
            theta_max: float = math.radians(45)
            V_Rd_max_max_angle: Quantity = (
                shear_eq.max_shear_resistance(alpha_cw, b_w, z, nu_1, f_cd, theta_max) * N
            ).to(kN)

            if self._V_Ed_1 > V_Rd_max_max_angle:
                self._theta = theta_max
                self._cot_theta = 1 / math.tan(self._theta)
                self._V_Rd_max = V_Rd_max_max_angle
                self._max_shear_ok = False
            else:
                self._max_shear_ok = True
                # Determine the angle θ of the strut based on the shear force
                self._theta = shear_eq.strut_angle(self._V_Ed_1.to(kN).magnitude, V_Rd_max_max_angle.to(kN).magnitude)
                self._cot_theta = 1 / math.tan(self._theta)
                self._V_Rd_max = (shear_eq.max_shear_resistance(alpha_cw, b_w, z, nu_1, f_cd, self._theta) * N).to(kN)


def _calculate_required_shear_reinforcement_EN_1992_2004(
    self: "RectangularBeam",
) -> None:
    """
    Calculate the required shear reinforcement area (A_v_req).
    """
    z = self._z.to(mm).magnitude
    f_ywd = self._f_ywd.to(MPa).magnitude
    # Calculate the required shear reinforcement area
    self._A_v_req = max(
        shear_eq.required_shear_reinforcement(self._V_Ed_2.to(N).magnitude, z, f_ywd, self._cot_theta)
        * mm,  # Required area based on shear force, as an area per unit length
        self._A_v_min,  # Minimum required area
    )
    self._V_Rd_s = (
        shear_eq.shear_reinforcement_resistance(self._A_v.to(mm).magnitude, z, f_ywd, self._cot_theta) * N
    ).to(kN)
    self._V_s_req = self._V_Rd_s
    # Maximum shear capacity is the same as the steel capacity
    self._V_Rd = self._V_Rd_s


def _check_shear_EN_1992_2004(self: "RectangularBeam", force: Forces, *, report: bool = True) -> DataFrame | None:
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
            n_legs_actual = self._stirrup_n * 2  # Ensure legs are even
            self._stirrup_s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs_actual - 1)
            (
                self._stirrup_s_max_l,
                self._stirrup_s_max_w,
            ) = max_stirrup_spacing_EN_1992_2004(self, self._alpha)

        self._DCRv = abs((self._V_Ed_2.to("kN").magnitude / self._V_Rd.to("kN").magnitude))

        if not report:
            return None

        # Design results
        results = {
            "Label": self.label,  # Beam label
            "Comb.": force.label,
            "Av,min": round(self._A_v_min.to("cm ** 2 / m").magnitude, 2),  # Minimum shear reinforcement area
            "Av,req": round(self._A_v_req.to("cm ** 2 / m").magnitude, 2),  # Required shear reinforcing area
            "Av": round(self._A_v.to("cm ** 2 / m").magnitude, 2),  # Provided stirrup reinforcement per unit length
            "NEd": self._N_Ed.to("kN").magnitude,
            "VEd,1": self._V_Ed_1.to("kN").magnitude,  # Max Vu for the design at the support
            "VEd,2": self._V_Ed_2.to("kN").magnitude,  # Max Vu for the design at d from the support
            "VRd,c": round(self._V_Rd_c.to("kN").magnitude, 2),  # Concrete contribution to shear capacity
            "VRd,s": round(self._V_Rd_s.to("kN").magnitude, 2),  # Reinforcement contribution to shear capacity
            "VRd": round(self._V_Rd.to("kN").magnitude, 2),  # Total shear capacity
            "VRd,max": round(self._V_Rd_max.to("kN").magnitude, 2),  # Maximum shear capacity
            "VEd,1≤VRd,max": self._max_shear_ok,  # Check if applied shear is within max shear capacity
            "VEd,2≤VRd": self._V_Ed_2.to("kN").magnitude
            <= self._V_Rd.to("kN").magnitude,  # Check if applied shear is within total capacity
            "DCR": round(self._DCRv, 3),
        }
        _initialize_dicts_EN_1992_2004_shear(self)
        return pd.DataFrame([results], index=[0])


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
        rho_min = flexure_eq.min_reinforcement_ratio(
            self.concrete.f_ctm.to(MPa).magnitude,
            self.steel_bar.f_y.to(MPa).magnitude,
        )
        # Set the maximum tensile reinforcement ratio
        rho_max = flexure_eq.max_reinforcement_ratio()

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


def _compression_zone_limits_EN_1992_2004(self: "RectangularBeam", d: Quantity) -> Tuple[Quantity, Quantity]:
    """Compression zone at the ductility limit, as (neutral axis, block depth).

    Both limits EN 1992-1-1 imposes are written on the NEUTRAL AXIS ratio
    ``x_u/d``: the redistribution limit of 5.5(4), ``(delta - k_1)/k_2``
    (``k_3``/``k_4`` above C50/60), and the 0.45 ductility cap. The equivalent
    rectangular stress block of 3.1.7(3) is ``lambda`` times shallower than the
    neutral axis, so the conversion happens once, here, and every caller gets
    both depths already in the right units.

    Returns
    -------
    (x_u_lim, x_eff_lim)
        Neutral axis depth and equivalent rectangular block depth at the limit.
    """
    concrete_en = cast("Concrete_EN_1992_2004", self.concrete)
    lambda_ = concrete_en._lambda_factor()
    xi_lim = flexure_eq.neutral_axis_depth_limit_ratio(  # x_u/d
        concrete_en._f_ck.to(MPa).magnitude,
        concrete_en._delta,
        concrete_en._k_1,
        concrete_en._k_2,
        concrete_en._k_3,
        concrete_en._k_4,
    )
    x_u_lim = xi_lim * d
    return x_u_lim, lambda_ * x_u_lim


def _calculate_flexural_reinforcement_EN_1992_2004(
    self: "RectangularBeam", M_Ed: Quantity, d: Quantity, d_prima: Quantity
) -> tuple[Quantity, Quantity, Quantity, Quantity]:
    """
    Calculate the required top and bottom reinforcement areas for bending.
    """
    rho_min, rho_max = _min_max_flexural_reinforcement_ratio_EN_1992_2004(self)
    # ADR-0005 boundary: convert once, compute in floats (N, mm, MPa, N·mm),
    # re-apply units on the way out.
    b = self.width.to(mm).magnitude
    d_mm = d.to(mm).magnitude
    A_s_min = rho_min * d_mm * b
    A_s_max = rho_max * d_mm * b

    # Constants and material properties
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        eta = self.concrete._eta_factor()  # Factor for concrete strength (EN 1992-1-1)
        f_cd = self._f_cd.to(MPa).magnitude
        # Define f_yd
        f_yd = (self.steel_bar.f_y / self.concrete._gamma_s).to(MPa).magnitude

        # Compression zone at the ductility limit (EC2 5.5(4) and the 0.45 cap)
        x_u_lim_q, x_eff_lim_q = _compression_zone_limits_EN_1992_2004(self, d)
        x_u_lim = x_u_lim_q.to(mm).magnitude
        x_eff_lim = x_eff_lim_q.to(mm).magnitude
        # Limit moment for compressive reinforcement
        M = M_Ed.to(_Nmm).magnitude
        M_lim = flexure_eq.limit_moment(eta, f_cd, b, x_eff_lim, d_mm)

        # Check if compressive reinforcement is required
        if M <= M_lim:
            # No compressive reinforcement required. Inverting
            # M = eta*f_cd*b*x_eff*(d - 0.5*x_eff) already yields lambda*x, so
            # lambda must NOT be applied again to the lever arm below.
            x_eff = flexure_eq.compression_block_depth_for_moment(M, b, d_mm, eta, f_cd)

            # Area of required tensile reinforcement, at least the minimum
            A_s1 = max(
                flexure_eq.reinforcement_for_moment(M, flexure_eq.lever_arm(d_mm, x_eff), f_yd),
                A_s_min,
            )
            # No compressive reinforcement required
            A_s2 = 0.0

        else:
            # Compressive reinforcement is required
            self._doubly_reinforced = True
            d_prime = d_prima.to(mm).magnitude
            # Limit tensile reinforcement area
            A_s1_lim = flexure_eq.reinforcement_for_moment(M_lim, flexure_eq.lever_arm(d_mm, x_eff_lim), f_yd)

            # Compressive reinforcement stress, from the strain at the neutral
            # axis at the plastic limit (EC2 §5.5)
            f_sd = flexure_eq.compression_steel_stress(
                x_u_lim,
                d_prime,
                self.concrete._epsilon_cu2,
                self.steel_bar._E_s.to(MPa).magnitude,
                f_yd,
            )

            # Extra moment to take with top reinforcement, on the lever arm
            # (d - d') of the compression steel couple
            A_s2 = flexure_eq.reinforcement_for_moment(M - M_lim, d_mm - d_prime, f_sd)

            # Required tensile reinforcement area
            A_s1 = max(A_s1_lim + A_s2, A_s_min)

    return A_s_min * _mm2, A_s_max * _mm2, A_s1 * _mm2, A_s2 * _mm2


def _simple_determine_nominal_moment_EN_1992_2004(
    self: "RectangularBeam",
    A_s: Quantity,
    d: Quantity,
    A_s_prime: Quantity,
    d_prime: Quantity,
) -> Quantity:
    """
    Design bending resistance M_Rd of a rectangular section.

    The reinforcement on the compression face is only part of the resisting
    couple when the section *needs* it to be in equilibrium, i.e. when the
    compression block required to balance the tension steel exceeds the
    ductility limit ``x_eff_lim``. Below that limit the concrete alone balances
    the tension steel: any bar on the compression face is then an erection or
    detailing bar, it is ignored (safe side), and the lever arm stays the singly
    reinforced one. Only above the limit does the concrete saturate at
    ``M_lim`` and the excess tension steel pair with the compression steel,
    which is what changes the lever arm to ``d - d'``.

    This mirrors :func:`_calculate_flexural_reinforcement_EN_1992_2004`, which
    sizes the tension steel with exactly the same criterion, so that a layout
    coming out of ``design_flexure`` is verified by ``check_flexure`` under the
    same mechanical model.

    For reference see Jimenez Montoya 16 Ed. P 213.
    """
    # Constants and material properties
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        f_yd = (self.steel_bar._f_y / self.concrete._gamma_s).to(MPa).magnitude
        f_cd = (self.concrete._alpha_cc * self.concrete.f_ck / self.concrete.gamma_c).to(MPa).magnitude
        eta = self.concrete._eta_factor()  # Factor for concrete strength
        b = self.width.to(mm).magnitude
        d_mm = d.to(mm).magnitude
        A_s_mm2 = A_s.to(_mm2).magnitude

        # Depth of the equivalent rectangular block that balances the whole
        # tension steel, assuming the section is singly reinforced.
        x_eff = flexure_eq.compression_block_depth_for_steel(A_s_mm2, f_yd, eta, f_cd, b)

        # Ductility limit -- same criterion as the reinforcement sizing routine.
        _, x_eff_lim_q = _compression_zone_limits_EN_1992_2004(self, d)
        x_eff_lim = x_eff_lim_q.to(mm).magnitude

        if x_eff <= x_eff_lim:
            # Singly reinforced: A_s_prime is not required for equilibrium.
            M_Rd = flexure_eq.moment_resistance_singly_reinforced(A_s_mm2, f_yd, d_mm, x_eff) * _Nmm
        else:
            # Doubly reinforced: the concrete contribution saturates at the
            # ductility limit and the excess tension steel is balanced by the
            # compression reinforcement, with lever arm (d - d').
            M_Rd = (
                flexure_eq.moment_resistance_doubly_reinforced(
                    A_s_mm2,
                    A_s_prime.to(_mm2).magnitude,
                    f_yd,
                    eta,
                    f_cd,
                    b,
                    d_mm,
                    d_prime.to(mm).magnitude,
                    x_eff_lim,
                )
                * _Nmm
            )
    return M_Rd


def _determine_nominal_moment_EN_1992_2004(self: "RectangularBeam", force: Forces) -> None:
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


# ─────────────────────────────────────────────────────────────────────────────
# Flexural design — EN 1992-2004 hooks for the shared driver
# ─────────────────────────────────────────────────────────────────────────────
# The design *strategy* -- Picard loop on the mechanical covers, discrete rebar
# selection, face reconciliation, cycle detection and final capacity
# verification -- lives in `mento.codes.flexure_design` and is shared with
# ACI 318-19. Only the two hooks below are EN-specific.


def _flexure_capacity_EN_1992_2004(self: "RectangularBeam", face: str, M_demand: Quantity) -> Quantity:
    """M_Rd of the layout currently applied to the section, on `face`.

    Recomputed rather than read from cached state so that the centroid of the
    layout just applied is taken into account.
    """
    M_abs = abs(M_demand)
    probe_force = Forces(M_y=(M_abs if face == "bot" else -M_abs))
    _determine_nominal_moment_EN_1992_2004(self, probe_force)
    return self._M_Rd_bot if face == "bot" else self._M_Rd_top


def _required_areas_EN_1992_2004(
    self: "RectangularBeam", face: str, M: Quantity, d: Quantity, d_prime: Quantity
) -> _FaceDemand:
    """Steel required by EN 1992-2004 on `face` for the moment `M`."""
    A_s_min, A_s_max, A_s1, A_s2 = _calculate_flexural_reinforcement_EN_1992_2004(self, M, d, d_prime)
    if face == "bot":
        self._A_s_min_bot, self._A_s_max_bot = A_s_min, A_s_max
    else:
        self._A_s_min_top, self._A_s_max_top = A_s_min, A_s_max
    return _FaceDemand(A_s_min, A_s_max, A_s1, A_s2)


def _design_flexure_EN_1992_2004(self: "RectangularBeam", max_M_y_bot: Quantity, max_M_y_top: Quantity) -> None:
    """Design the longitudinal reinforcement of a beam per EN 1992-2004.

    Thin wrapper: everything that is not an EN equation lives in
    ``mento.codes.flexure_design``.
    """
    _initialize_variables_EN_1992_2004(self)

    def _required(face: str, M: Quantity, d: Quantity, d_prime: Quantity) -> _FaceDemand:
        return _required_areas_EN_1992_2004(self, face, M, d, d_prime)

    def _capacity(face: str, M: Quantity) -> Quantity:
        return _flexure_capacity_EN_1992_2004(self, face, M)

    _run_flexure_design(self, max_M_y_bot, max_M_y_top, _required, _capacity)


def _check_flexure_EN_1992_2004(self: "RectangularBeam", force: Forces, *, report: bool = True) -> pd.DataFrame | None:
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
        ) = _calculate_flexural_reinforcement_EN_1992_2004(self, self._M_Ed_bot, self._d_bot, self._c_mec_top)
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

    if not report:
        return None

    # Compile the design results into a dictionary.
    results = _compile_results_EN_1992_2004_flexure_metric(self, force)

    # Initialize any additional dictionaries required for ACI 318-19 flexural checks.
    _initialize_dicts_EN_1992_2004_flexure(self)

    # Return the results as a Pandas DataFrame.
    return pd.DataFrame([results], index=[0])


def _compile_results_EN_1992_2004_flexure_metric(self: "RectangularBeam", force: Forces) -> Dict[str, Any]:
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
            self._stirrup_d_b = 0 * mm if self.concrete.unit_system == "metric" else 0 * inch
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
            "✅" if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else "❌"
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_shear_checks_passed = all(check == "✅" for check in checks)
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
        check_max = "✅" if self._max_shear_ok else "❌"
        check_DCR = "✅" if self._DCRv < 1 else "❌"
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
        # Update longitudinal rebar attributes
        self._update_longitudinal_rebar_attributes()
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
            else:
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
                "MRd",
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
                "MRd",
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
            self._all_flexure_checks_passed and (check_DCR_bot == "✅") and (check_DCR_top == "✅")
        )
