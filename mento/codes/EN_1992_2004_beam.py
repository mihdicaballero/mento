import math
from pint import Quantity
from typing import TYPE_CHECKING, Tuple, cast

# from devtools import debug


from mento.codes.en_1992_2004.equations import flexure as flexure_eq
from mento.codes.check_state import (
    ENFlexureCheckState,
    ENShearCheckState,
    apply_en_shear_state,
    new_en_flexure_state,
    new_en_shear_state,
)
from mento.codes.en_1992_2004.equations import shear as shear_eq
from mento.codes.flexure_design import _FaceDemand, _run_flexure_design
from mento.material import Concrete_EN_1992_2004
from mento.rebar import max_stirrup_spacing_EN_1992_2004
from mento.units import MPa, mm, kNm, dimensionless, kN, N
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


def _initialize_shear_variables_EN_1992_2004(self: "RectangularBeam", st: ENShearCheckState, force: Forces) -> None:
    """
    Initialize variables for EN 1992-2004 design code.
    """
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Set the initial variables
        st.N_Ed = force.N_x
        st.V_Ed_1 = force.V_z  # Consider the same shear at the edge of support and in d
        st.V_Ed_2 = force.V_z  # Consider the same shear at the edge of support and in d

        # Minimum shear reinforcement calculation. Eq. (9.5N) gives the ratio;
        # §9.2.2(5) defines it as A_sw/(s*b_w*sin(alpha)), hence the geometry here.
        rho_min = shear_eq.min_shear_reinforcement_ratio(
            self.concrete.f_ck.to(MPa).magnitude,
            self._f_ywk.to(MPa).magnitude,
        )
        st.A_v_min = rho_min * self.width * math.sin(self._alpha)

        # Consider bottom or top tension reinforcement
        st.A_s_tension = self._A_s_bot if force._M_y >= 0 * kNm else self._A_s_top
        # Compression stress, positive
        if force._M_y >= 0 * kNm:
            st.rho_l_bot = min(
                (st.A_s_tension.to("cm**2") + self._A_p.to("cm**2")) / (self.width.to("cm") * self._d_shear.to("cm")),
                0.02 * dimensionless,
            )
        else:
            st.rho_l_top = min(
                (st.A_s_tension + self._A_p) / (self.width * self._d_shear),
                0.02 * dimensionless,
            )

        # Shear calculation for sections without rebar
        st.k_value = shear_eq.size_effect_factor(self._d_shear.to(mm).magnitude)

        # Positive of compression
        st.sigma_cp = (
            shear_eq.axial_stress(
                st.N_Ed.to(N).magnitude,
                self.A_x.to(_mm2).magnitude,
                st.f_cd_shear.to(MPa).magnitude,
            )
            * MPa
        )


def _shear_without_rebar_EN_1992_2004(self: "RectangularBeam", st: ENShearCheckState) -> Quantity:
    st.theta = 0
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Total shear capacity without rebar: Eq. (6.2.a), floored by Eq. (6.2.b).
        f_ck = self.concrete.f_ck.to(MPa).magnitude
        b_w = self.width.to(mm).magnitude
        d = self._d_shear.to(mm).magnitude
        sigma_cp = st.sigma_cp.to(MPa).magnitude
        rho_l = st.rho_l_bot if self._M_Ed >= 0 * kNm else st.rho_l_top
        V_Rd_c_min = (shear_eq.min_shear_resistance_without_reinforcement(f_ck, st.k_value, sigma_cp, b_w, d) * N).to(
            kN
        )
        V_Rd_c = (
            shear_eq.shear_resistance_without_reinforcement(
                f_ck,
                self.concrete.gamma_c,
                st.k_value,
                rho_l.to(dimensionless).magnitude,
                sigma_cp,
                b_w,
                d,
            )
            * N
        ).to(kN)
    return max(V_Rd_c_min, V_Rd_c)


def _calculate_max_shear_strength_EN_1992_2004(self: "RectangularBeam", st: ENShearCheckState) -> None:
    """
    Calculate the maximum shear strength (V_Rd_max) and the corresponding strut angle (θ).
    """
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        alpha_cw = 1.0  # Non-prestressed members or members subject to tensile stress due to axial force
        nu_1 = shear_eq.strut_strength_reduction_factor(self.concrete.f_ck.to(MPa).magnitude)
        b_w = self.width.to(mm).magnitude
        f_cd = st.f_cd_shear.to(MPa).magnitude
        z = shear_eq.lever_arm(self._d_shear.to(mm).magnitude)
        st.z = z * mm  # Lever arm

        # The θ angle is limited between 21.8° ≤ θ ≤ 45° (1 ≤ cot(θ) ≤ 2.5)
        # Check the minimum strut angle θ = 21.8° (cot(θ) = 2.5)
        theta_min: float = math.radians(21.8)
        cot_theta_min: float = 1 / math.tan(theta_min)
        V_Rd_max_min_angle = (shear_eq.max_shear_resistance(alpha_cw, b_w, z, nu_1, f_cd, theta_min) * N).to(kN)

        if st.V_Ed_1 <= V_Rd_max_min_angle:
            # If within the minimum angle
            st.theta = theta_min
            st.cot_theta = cot_theta_min
            st.V_Rd_max = V_Rd_max_min_angle
            st.max_shear_ok = True
        else:
            # Check the maximum strut angle θ = 45° (cot(θ) = 1.0)
            theta_max: float = math.radians(45)
            V_Rd_max_max_angle: Quantity = (
                shear_eq.max_shear_resistance(alpha_cw, b_w, z, nu_1, f_cd, theta_max) * N
            ).to(kN)

            if st.V_Ed_1 > V_Rd_max_max_angle:
                st.theta = theta_max
                st.cot_theta = 1 / math.tan(st.theta)
                st.V_Rd_max = V_Rd_max_max_angle
                st.max_shear_ok = False
            else:
                st.max_shear_ok = True
                # Determine the angle θ of the strut based on the shear force
                st.theta = shear_eq.strut_angle(st.V_Ed_1.to(kN).magnitude, V_Rd_max_max_angle.to(kN).magnitude)
                st.cot_theta = 1 / math.tan(st.theta)
                st.V_Rd_max = (shear_eq.max_shear_resistance(alpha_cw, b_w, z, nu_1, f_cd, st.theta) * N).to(kN)


def _calculate_required_shear_reinforcement_EN_1992_2004(self: "RectangularBeam", st: ENShearCheckState) -> None:
    """
    Calculate the required shear reinforcement area (A_v_req).
    """
    z = st.z.to(mm).magnitude
    f_ywd = st.f_ywd.to(MPa).magnitude
    # Calculate the required shear reinforcement area
    st.A_v_req = max(
        shear_eq.required_shear_reinforcement(st.V_Ed_2.to(N).magnitude, z, f_ywd, st.cot_theta)
        * mm,  # Required area based on shear force, as an area per unit length
        st.A_v_min,  # Minimum required area
    )
    st.V_Rd_s = (shear_eq.shear_reinforcement_resistance(self._A_v.to(mm).magnitude, z, f_ywd, st.cot_theta) * N).to(kN)
    st.V_s_req = st.V_Rd_s
    # Maximum shear capacity is the same as the steel capacity
    st.V_Rd = st.V_Rd_s


def _check_shear_EN_1992_2004(self: "RectangularBeam", force: Forces) -> ENShearCheckState:
    """Run the EN shear check for one combination and return what it found.

    Nothing is written to the section; see the ACI counterpart for the split
    between this and the reporting path.
    """
    st = new_en_shear_state(self)
    concrete = cast("Concrete_EN_1992_2004", self.concrete)

    # The material values the shear check needs. _f_cd belongs to flexure as
    # well, so it is carried rather than written.
    st.f_ywd = self._f_ywk / concrete._gamma_s
    alpha_cc_shear = 1  # Take this as 1.00 for shear design and not 0.85, as in Eurocode Applied.
    st.f_cd_shear = alpha_cc_shear * concrete.f_ck / concrete.gamma_c
    st.f_cd = concrete._alpha_cc * concrete.f_ck / concrete.gamma_c

    _initialize_shear_variables_EN_1992_2004(self, st, force)

    if self._stirrup_n == 0:
        # The assumed stirrup diameter stays: see the note in the ACI check.
        st.V_Rd_c = _shear_without_rebar_EN_1992_2004(self, st)
        # According to EN1992-1-1 §6.2.1(4) minimum shear reinforcement should nevertheless be provided
        # according to EN1992-1-1 §9.2.2. The minimum shear reinforcement may be omitted in members where
        # transverse redistribution of loads is possible (such as slabs) and members of minor importance
        # which do not contribute significantly to the overall resistance and stability of the structure.
        st.A_v_req = st.A_v_min
        # Maximum shear capacity is the same as the concrete capacity
        st.V_Rd = st.V_Rd_c
        st.V_Rd_max = st.V_Rd
        st.max_shear_ok = st.V_Ed_1 <= st.V_Rd_max
    else:
        # The provided stirrup area per unit length is what the section already
        # carries; the setter computed it from the same geometry.
        _calculate_max_shear_strength_EN_1992_2004(self, st)
        _calculate_required_shear_reinforcement_EN_1992_2004(self, st)

        # Rebar spacing checks
        n_legs_actual = self._stirrup_n * 2  # Ensure legs are even
        st.stirrup_s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs_actual - 1)
        (
            st.stirrup_s_max_l,
            st.stirrup_s_max_w,
        ) = max_stirrup_spacing_EN_1992_2004(self, self._alpha)

    st.DCR = abs((st.V_Ed_2.to("kN").magnitude / st.V_Rd.to("kN").magnitude))
    return st


def _design_shear_EN_1992_2004(self: "RectangularBeam", force: Forces) -> None:
    """Size the shear reinforcement. Designing is meant to change the section,
    so this runs the same helpers over a state and then applies it."""
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Initialize all the code related variables
        _initialize_variables_EN_1992_2004(self)
        # Split bottom and top moments. Designing changes the section, so the
        # moment split is applied rather than kept local.
        flexure_st = new_en_flexure_state(self)
        _split_top_bot_moment(self, flexure_st, force)
        self._M_Ed, self._M_Ed_bot, self._M_Ed_top = (
            flexure_st.M_Ed,
            flexure_st.M_Ed_bot,
            flexure_st.M_Ed_top,
        )
        st = new_en_shear_state(self)
        st.f_ywd = self._f_ywd
        st.f_cd_shear = self._f_cd_shear
        st.f_cd = self._f_cd
        _initialize_shear_variables_EN_1992_2004(self, st, force)
        # Calculate maximum shear strength
        _calculate_max_shear_strength_EN_1992_2004(self, st)
        # Calculate required shear reinforcement area
        _calculate_required_shear_reinforcement_EN_1992_2004(self, st)
        apply_en_shear_state(self, st)

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
        # Derived from the concrete rather than read off the beam: it is a
        # material property, so a check has no reason to have stored it first.
        f_cd = (self.concrete._alpha_cc * self.concrete.f_ck / self.concrete.gamma_c).to(MPa).magnitude
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


def _determine_nominal_moment_EN_1992_2004(self: "RectangularBeam", st: ENFlexureCheckState, force: Forces) -> None:
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
    st.A_s_min_bot = rho_min_bot * self._d_bot * self.width
    st.A_s_max_bot = rho_max * self._d_bot * self.width
    # Determine the nominal moment for positive moments
    st.M_Rd_bot = _simple_determine_nominal_moment_EN_1992_2004(
        self, self._A_s_bot, self._d_bot, self._A_s_top, self._c_mec_top
    )
    # Determine capacity for negative moment (tension at the top)
    st.A_s_min_top = rho_min_top * self._d_top * self.width
    st.A_s_max_top = rho_max * self._d_top * self.width
    st.M_Rd_top = _simple_determine_nominal_moment_EN_1992_2004(
        self, self._A_s_top, self._d_top, self._A_s_bot, self._c_mec_bot
    )
    return None


def _split_top_bot_moment(self: "RectangularBeam", st: ENFlexureCheckState, force: Forces) -> None:
    st.M_Ed = force._M_y
    if st.M_Ed > 0 * kNm:
        st.M_Ed_bot = st.M_Ed
        st.M_Ed_top = 0 * kNm
    else:
        st.M_Ed_bot = 0 * kNm
        st.M_Ed_top = st.M_Ed


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
    st = new_en_flexure_state(self)
    _determine_nominal_moment_EN_1992_2004(self, st, probe_force)
    # The design wants the capacity on the beam too; a check never comes here.
    self._M_Rd_bot, self._M_Rd_top = st.M_Rd_bot, st.M_Rd_top
    self._A_s_min_bot, self._A_s_max_bot = st.A_s_min_bot, st.A_s_max_bot
    return st.M_Rd_bot if face == "bot" else st.M_Rd_top


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


def _check_flexure_EN_1992_2004(self: "RectangularBeam", force: Forces) -> ENFlexureCheckState:
    """Check the flexural capacity per EN 1992-1-1 and return what it found.

    Nothing is written to the section; only the reporting path copies the
    result back. See the ACI counterpart for the same split.
    """
    st = new_en_flexure_state(self)
    concrete = cast("Concrete_EN_1992_2004", self.concrete)

    # The material values this check needs, carried rather than written.
    st.f_ywd = self._f_ywk / concrete._gamma_s
    st.f_cd_shear = concrete.f_ck / concrete.gamma_c
    st.f_cd = concrete._alpha_cc * concrete.f_ck / concrete.gamma_c

    # Split bottom and top moments
    _split_top_bot_moment(self, st, force)

    # Calculate the nominal moments for both top and bottom reinforcement.
    _determine_nominal_moment_EN_1992_2004(self, st, force)
    if st.M_Ed >= 0:
        # For positive moments, calculate the reinforcement requirements for the bottom tension side.
        (
            st.A_s_min_bot,
            st.A_s_max_bot,
            st.A_s_req_bot,
            st.A_s_req_top,
        ) = _calculate_flexural_reinforcement_EN_1992_2004(self, st.M_Ed_bot, self._d_bot, self._c_mec_top)
        st.c_d_top = 0
        # Calculate the design capacity ratio for the bottom side.
        st.DCR_bot = round(
            st.M_Ed_bot.to("kN*m").magnitude / st.M_Rd_bot.to("kN*m").magnitude,
            3,
        )
        st.DCR_top = 0
    else:
        # For negative moments, calculate the reinforcement requirements for the top tension side.
        (
            st.A_s_min_top,
            st.A_s_max_top,
            st.A_s_req_top,
            st.A_s_req_bot,
        ) = _calculate_flexural_reinforcement_EN_1992_2004(
            self, abs(st.M_Ed_top / kNm) * kNm, self._d_top, self._c_mec_bot
        )
        st.c_d_bot = 0
        # Calculate the design capacity ratio for the top side.
        st.DCR_top = round(
            -st.M_Ed_top.to("kN*m").magnitude / st.M_Rd_top.to("kN*m").magnitude,
            3,
        )
        st.DCR_bot = 0

    # Determine the maximum detailing cover dimensions for top and bottom.
    st.d_b_max_top = max(self._d_b1_t, self._d_b2_t, self._d_b3_t, self._d_b4_t)
    st.d_b_max_bot = max(self._d_b1_b, self._d_b2_b, self._d_b3_b, self._d_b4_b)

    # Calculate the longitudinal reinforcement ratios for both sides.
    st.rho_l_bot = self._A_s_bot / (self._d_bot * self.width)
    st.rho_l_top = self._A_s_bot / (self._d_top * self.width)

    return st


##########################################################
# RESULTS
##########################################################
