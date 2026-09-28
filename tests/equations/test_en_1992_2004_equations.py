"""EN 1992-2004 equations, checked against values read off the code itself.

These are the Phase 1 equation tests: each case is a hand-computable evaluation
of the printed clause, so a reader with EN 1992-1-1 open can verify a row
without running mento. Nothing here imports pint or touches a beam.

EN is written in a single coherent SI system, so — unlike the ACI equations —
none of these functions takes an ``is_imperial`` flag. Everything is N, mm,
mm², MPa and N·mm.
"""

import math

import pytest

from mento.codes.en_1992_2004.equations import flexure as flex
from mento.codes.en_1992_2004.equations import shear as eq


# ===========================================================================
# SHEAR — EN 1992-1-1 §6.2, §9.2.2
# ===========================================================================


# ---------------------------------------------------------------------------
# k — §6.2.2(1)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "d_mm, expected",
    [
        # d = 200 mm is the transition: 1 + sqrt(200/200) = 2, exactly the cap.
        (200.0, 2.0),
        (50.0, 2.0),  # 1 + sqrt(4) = 3 -> capped, not > 2
        (800.0, 1.5),  # 1 + sqrt(0.25)
        (2000.0, 1 + math.sqrt(0.1)),
    ],
)
def test_size_effect_factor(d_mm, expected):
    assert eq.size_effect_factor(d_mm) == pytest.approx(expected, rel=1e-9)


def test_size_effect_factor_never_exceeds_two():
    # The cap is the point of the clause: a very shallow member gets 2, not more.
    assert all(eq.size_effect_factor(d) <= 2.0 for d in range(10, 3000, 10))


# ---------------------------------------------------------------------------
# sigma_cp — §6.2.2(1)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "N_Ed, A_c, f_cd, expected",
    [
        (0.0, 100_000.0, 16.0, 0.0),  # no axial load
        (200_000.0, 100_000.0, 16.0, 2.0),  # 200e3/100e3 = 2 MPa, under 0.2*16 = 3.2
        (1_000_000.0, 100_000.0, 16.0, 3.2),  # would be 10 MPa; capped at 0.2*f_cd
    ],
)
def test_axial_stress(N_Ed, A_c, f_cd, expected):
    assert eq.axial_stress(N_Ed, A_c, f_cd) == pytest.approx(expected, rel=1e-9)


# ---------------------------------------------------------------------------
# rho_w,min — Eq. (9.5N)
# ---------------------------------------------------------------------------


def test_min_shear_reinforcement_ratio():
    # C25/30 with B500S: 0.08*sqrt(25)/500 = 0.08*5/500 = 8.0e-4.
    # On a 200 mm web that is 8.0e-4*200 = 0.16 mm²/mm = 1.6 cm²/m.
    rho = eq.min_shear_reinforcement_ratio(25.0, 500.0)
    assert rho == pytest.approx(8.0e-4, rel=1e-12)
    assert rho * 200.0 == pytest.approx(0.16, rel=1e-12)


def test_min_shear_reinforcement_ratio_grows_with_concrete_strength():
    # sqrt(f_ck) in the numerator: stronger concrete cracks at a higher load,
    # so more stirrups are needed to hold the crack once it forms.
    assert eq.min_shear_reinforcement_ratio(50.0, 500.0) == pytest.approx(
        eq.min_shear_reinforcement_ratio(25.0, 500.0) * math.sqrt(2), rel=1e-12
    )


# ---------------------------------------------------------------------------
# V_Rd,c — Eqs. (6.2.a), (6.2.b), (6.3N)
# ---------------------------------------------------------------------------


def test_shear_resistance_without_reinforcement():
    # f_ck = 25, gamma_c = 1.5 -> C_Rd,c = 0.18/1.5 = 0.12.
    # k = 2, rho_l = 0.01, sigma_cp = 0, b_w = 200, d = 500:
    #   (100*0.01*25)**(1/3) = 25**(1/3) = 2.924017...
    #   0.12*2*2.924017 = 0.7017642 MPa
    #   * 200 * 500 = 70176.4 N
    expected = 0.12 * 2 * 25 ** (1 / 3) * 200 * 500
    got = eq.shear_resistance_without_reinforcement(25.0, 1.5, 2.0, 0.01, 0.0, 200.0, 500.0)
    assert got == pytest.approx(expected, rel=1e-9)
    assert got == pytest.approx(70_176.4, rel=1e-6)


def test_shear_resistance_without_reinforcement_adds_axial_term():
    # k_1 = 0.15 is the recommended value; sigma_cp enters as a flat stress.
    base = eq.shear_resistance_without_reinforcement(25.0, 1.5, 2.0, 0.01, 0.0, 200.0, 500.0)
    with_axial = eq.shear_resistance_without_reinforcement(25.0, 1.5, 2.0, 0.01, 2.0, 200.0, 500.0)
    assert with_axial - base == pytest.approx(eq.K_1 * 2.0 * 200.0 * 500.0, rel=1e-9)


def test_shear_resistance_without_reinforcement_scales_with_cube_root_of_rho():
    # Eq. (6.2.a) is rho_l**(1/3): eight times the steel doubles V_Rd,c.
    small = eq.shear_resistance_without_reinforcement(25.0, 1.5, 2.0, 0.0025, 0.0, 200.0, 500.0)
    large = eq.shear_resistance_without_reinforcement(25.0, 1.5, 2.0, 0.02, 0.0, 200.0, 500.0)
    assert large == pytest.approx(2 * small, rel=1e-9)


def test_min_shear_resistance_without_reinforcement():
    # v_min = 0.035*k**1.5*sqrt(f_ck) = 0.035*2**1.5*5 = 0.49497475 MPa,
    # times b_w*d = 200*500 = 100000 mm² -> 49497.5 N.
    expected = 0.035 * 2**1.5 * 5.0 * 200 * 500
    got = eq.min_shear_resistance_without_reinforcement(25.0, 2.0, 0.0, 200.0, 500.0)
    assert got == pytest.approx(expected, rel=1e-9)
    assert got == pytest.approx(49_497.47, rel=1e-6)


def test_min_shear_resistance_is_the_floor_for_light_reinforcement():
    # The point of Eq. (6.2.b): at low rho_l, Eq. (6.2.a) drops below it.
    floor = eq.min_shear_resistance_without_reinforcement(25.0, 2.0, 0.0, 200.0, 500.0)
    light = eq.shear_resistance_without_reinforcement(25.0, 1.5, 2.0, 0.0005, 0.0, 200.0, 500.0)
    heavy = eq.shear_resistance_without_reinforcement(25.0, 1.5, 2.0, 0.02, 0.0, 200.0, 500.0)
    assert light < floor < heavy


# ---------------------------------------------------------------------------
# nu_1 and z — Eq. (6.6N), §6.2.3(1)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "f_ck, expected",
    [
        (25.0, 0.54),  # 0.6*(1 - 25/250) = 0.6*0.9
        (50.0, 0.48),  # 0.6*(1 - 0.2)
        (60.0, 0.6 * 0.76),
    ],
)
def test_strut_strength_reduction_factor(f_ck, expected):
    assert eq.strut_strength_reduction_factor(f_ck) == pytest.approx(expected, rel=1e-9)


def test_shear_lever_arm():
    assert eq.lever_arm(600.0) == pytest.approx(540.0, rel=1e-12)


# ---------------------------------------------------------------------------
# V_Rd,max and theta — Eq. (6.9)
# ---------------------------------------------------------------------------


def test_max_shear_resistance_at_45_degrees():
    # alpha_cw*b_w*z*nu_1*f_cd = 1*300*500*0.5*20 = 1.5e6 N.
    # At theta = 45 deg, cot + tan = 2, so V_Rd,max = 750 kN.
    got = eq.max_shear_resistance(1.0, 300.0, 500.0, 0.5, 20.0, math.radians(45))
    assert got == pytest.approx(750_000.0, rel=1e-9)


def test_max_shear_resistance_at_the_shallowest_strut():
    # theta = 21.8 deg is the flattest strut §6.2.3(2) allows (cot ~ 2.5):
    #   cot + tan = 2.49984 + 0.40003 = 2.89987
    #   1.5e6 / 2.89987 = 517264 N
    theta = math.radians(21.8)
    expected = 1.5e6 / (1 / math.tan(theta) + math.tan(theta))
    got = eq.max_shear_resistance(1.0, 300.0, 500.0, 0.5, 20.0, theta)
    assert got == pytest.approx(expected, rel=1e-9)
    # A flatter strut crushes earlier than a 45 deg one.
    assert got < eq.max_shear_resistance(1.0, 300.0, 500.0, 0.5, 20.0, math.radians(45))


@pytest.mark.parametrize(
    "V_Ed, V_Rd_max_45, expected_deg",
    [
        (750_000.0, 750_000.0, 45.0),  # asin(1)/2 = 45 deg
        (375_000.0, 750_000.0, 15.0),  # asin(0.5)/2 = 15 deg
    ],
)
def test_strut_angle(V_Ed, V_Rd_max_45, expected_deg):
    assert math.degrees(eq.strut_angle(V_Ed, V_Rd_max_45)) == pytest.approx(expected_deg, rel=1e-9)


def test_strut_angle_inverts_eq_6_9():
    # The angle chosen for V_Ed must be the angle at which V_Rd,max equals V_Ed.
    V_Rd_max_45 = eq.max_shear_resistance(1.0, 300.0, 500.0, 0.5, 20.0, math.radians(45))
    V_Ed = 0.6 * V_Rd_max_45
    theta = eq.strut_angle(V_Ed, V_Rd_max_45)
    assert eq.max_shear_resistance(1.0, 300.0, 500.0, 0.5, 20.0, theta) == pytest.approx(V_Ed, rel=1e-9)


# ---------------------------------------------------------------------------
# A_sw/s and V_Rd,s — Eq. (6.8)
# ---------------------------------------------------------------------------


def test_shear_reinforcement_resistance():
    # A_sw/s = 0.2 mm²/mm (= 2 cm²/m), z = 500 mm, f_ywd = 400 MPa, cot = 2.5:
    #   0.2*500*400*2.5 = 100 kN
    assert eq.shear_reinforcement_resistance(0.2, 500.0, 400.0, 2.5) == pytest.approx(100_000.0, rel=1e-9)


def test_required_shear_reinforcement_inverts_the_resistance():
    assert eq.required_shear_reinforcement(100_000.0, 500.0, 400.0, 2.5) == pytest.approx(0.2, rel=1e-9)


def test_required_shear_reinforcement_falls_with_a_flatter_strut():
    # cot(theta) multiplies the number of stirrup legs the crack crosses, so
    # the flattest strut the code allows needs the least reinforcement.
    steep = eq.required_shear_reinforcement(100_000.0, 500.0, 400.0, 1.0)
    flat = eq.required_shear_reinforcement(100_000.0, 500.0, 400.0, 2.5)
    assert flat == pytest.approx(steep / 2.5, rel=1e-9)


# ===========================================================================
# FLEXURE — EN 1992-1-1 §3.1.7, §5.5, §6.1, §9.2.1
# ===========================================================================


# ---------------------------------------------------------------------------
# rho_min, rho_max — §9.2.1.1
# ---------------------------------------------------------------------------


def test_min_reinforcement_ratio_formula_governs():
    # C25/30: f_ctm = 0.3*25**(2/3) = 2.5650 MPa, f_yk = 500 MPa.
    #   0.26*2.5650/500 = 1.3338e-3, above the 0.0013 floor.
    f_ctm = 0.3 * 25 ** (2 / 3)
    assert flex.min_reinforcement_ratio(f_ctm, 500.0) == pytest.approx(0.26 * f_ctm / 500.0, rel=1e-9)
    assert flex.min_reinforcement_ratio(f_ctm, 500.0) > 0.0013


def test_min_reinforcement_ratio_floor_governs():
    # Weak concrete with strong steel: 0.26*1.6/500 = 8.32e-4, below the floor.
    assert flex.min_reinforcement_ratio(1.6, 500.0) == pytest.approx(0.0013, rel=1e-12)


def test_max_reinforcement_ratio():
    assert flex.max_reinforcement_ratio() == 0.04


# ---------------------------------------------------------------------------
# x_u/d at the ductility limit — §5.5(4)
# ---------------------------------------------------------------------------


def test_neutral_axis_depth_limit_ratio_below_c50():
    # mento's defaults for C25/30: delta = 0.85, k_1 = 0.44,
    # k_2 = 1.25*(0.6 + 0.0014/0.0035) = 1.25*1.0 = 1.25.
    #   (0.85 - 0.44)/1.25 = 0.328, below the 0.45 cap.
    k_2 = 1.25 * (0.6 + 0.0014 / 0.0035)
    got = flex.neutral_axis_depth_limit_ratio(25.0, 0.85, 0.44, k_2, 0.54, k_2)
    assert got == pytest.approx(0.328, rel=1e-9)


def test_neutral_axis_depth_limit_ratio_above_c50_uses_k3_k4():
    # C60/75: eps_cu2 = (2.6 + 35*((90-60)/100)**4)*1e-3 = 2.8835e-3, so
    # k_4 = 1.25*(0.6 + 0.0014/0.0028835) = 1.3568935.
    #   (0.85 - 0.54)/1.3568935 = 0.228461
    epsilon_cu2 = (2.6 + 35 * 0.3**4) * 1e-3
    k_4 = 1.25 * (0.6 + 0.0014 / epsilon_cu2)
    got = flex.neutral_axis_depth_limit_ratio(60.0, 0.85, 0.44, 1.25, 0.54, k_4)
    assert got == pytest.approx((0.85 - 0.54) / k_4, rel=1e-9)
    assert got == pytest.approx(0.228461, rel=1e-5)


def test_neutral_axis_depth_limit_ratio_is_capped_at_0_45():
    # No redistribution (delta = 1) would allow 0.56; the ductility cap wins.
    assert flex.neutral_axis_depth_limit_ratio(25.0, 1.0, 0.44, 1.0, 0.54, 1.0) == pytest.approx(0.45, rel=1e-12)


# ---------------------------------------------------------------------------
# M_lim and its inversion — §3.1.7(3)
# ---------------------------------------------------------------------------


def test_limit_moment():
    # eta = 1, f_cd = 20 MPa, b = 300 mm, x_eff_lim = 200 mm, d = 500 mm:
    #   C = 1*20*300*200 = 1.2e6 N on a lever arm of 500 - 100 = 400 mm
    #   -> 480e6 N·mm = 480 kNm
    assert flex.limit_moment(1.0, 20.0, 300.0, 200.0, 500.0) == pytest.approx(480e6, rel=1e-9)


def test_compression_block_depth_for_moment_inverts_the_block():
    # K = 480e6/(300*500²*1*20) = 0.32; 1 - 2K = 0.36; sqrt = 0.6
    #   x_eff = 500*(1 - 0.6) = 200 mm, which is exactly the block that
    #   produced 480 kNm in test_limit_moment.
    assert flex.compression_block_depth_for_moment(480e6, 300.0, 500.0, 1.0, 20.0) == pytest.approx(200.0, rel=1e-9)


def test_compression_block_depth_for_moment_is_zero_for_zero_moment():
    assert flex.compression_block_depth_for_moment(0.0, 300.0, 500.0, 1.0, 20.0) == pytest.approx(0.0, abs=1e-12)


def test_flexure_lever_arm():
    assert flex.lever_arm(500.0, 200.0) == pytest.approx(400.0, rel=1e-12)


def test_reinforcement_for_moment():
    # 480e6 N·mm on a 400 mm arm at 400 MPa -> 3000 mm² = 30 cm²
    assert flex.reinforcement_for_moment(480e6, 400.0, 400.0) == pytest.approx(3000.0, rel=1e-9)


# ---------------------------------------------------------------------------
# Compression steel — §6.1, §3.2.7(2)
# ---------------------------------------------------------------------------


def test_compression_steel_stress_yields():
    # x_u = 200, d' = 40: eps_s' = (160/200)*0.0035 = 0.0028,
    # times 200 GPa = 560 MPa, above f_yd = 434.78 -> capped at yield.
    got = flex.compression_steel_stress(200.0, 40.0, 0.0035, 200_000.0, 434.78)
    assert got == pytest.approx(434.78, rel=1e-12)


def test_compression_steel_stress_stays_elastic():
    # x_u = 100, d' = 50: eps_s' = 0.5*0.0035 = 0.00175 -> 350 MPa < f_yd.
    got = flex.compression_steel_stress(100.0, 50.0, 0.0035, 200_000.0, 434.78)
    assert got == pytest.approx(350.0, rel=1e-9)


# ---------------------------------------------------------------------------
# M_Rd — §6.1
# ---------------------------------------------------------------------------


def test_compression_block_depth_for_steel():
    # 3000 mm² at 400 MPa balanced by 1*20 MPa over a 300 mm width:
    #   x_eff = 3000*400/(20*300) = 200 mm
    assert flex.compression_block_depth_for_steel(3000.0, 400.0, 1.0, 20.0, 300.0) == pytest.approx(200.0, rel=1e-9)


def test_moment_resistance_singly_reinforced():
    # The mirror image of test_reinforcement_for_moment: 3000 mm² at 400 MPa
    # on a 400 mm arm gives back 480 kNm.
    assert flex.moment_resistance_singly_reinforced(3000.0, 400.0, 500.0, 200.0) == pytest.approx(480e6, rel=1e-9)


def test_moment_resistance_doubly_reinforced():
    # A_s,lim = 1*20*300*200/400 = 3000 mm², M_lim = 480e6 N·mm.
    # A_s = 4000 leaves 1000 mm² to pair with the 1000 mm² on top, on the
    # lever arm d - d' = 450 mm: 1000*400*450 = 180e6 N·mm.
    got = flex.moment_resistance_doubly_reinforced(4000.0, 1000.0, 400.0, 1.0, 20.0, 300.0, 500.0, 50.0, 200.0)
    assert got == pytest.approx(660e6, rel=1e-9)


def test_moment_resistance_doubly_reinforced_is_capped_by_the_top_steel():
    # 5000 mm² of tension steel has 2000 mm² beyond the limit, but only
    # 1000 mm² sits on the compression face, so the couple cannot grow.
    with_4000 = flex.moment_resistance_doubly_reinforced(4000.0, 1000.0, 400.0, 1.0, 20.0, 300.0, 500.0, 50.0, 200.0)
    with_5000 = flex.moment_resistance_doubly_reinforced(5000.0, 1000.0, 400.0, 1.0, 20.0, 300.0, 500.0, 50.0, 200.0)
    assert with_5000 == pytest.approx(with_4000, rel=1e-12)


def test_moment_resistance_is_continuous_at_the_ductility_limit():
    # At A_s = A_s,lim the compression couple vanishes and the doubly
    # reinforced branch must reproduce the singly reinforced value exactly.
    A_s_lim = 1.0 * 20.0 * 300.0 * 200.0 / 400.0
    singly = flex.moment_resistance_singly_reinforced(A_s_lim, 400.0, 500.0, 200.0)
    doubly = flex.moment_resistance_doubly_reinforced(A_s_lim, 1000.0, 400.0, 1.0, 20.0, 300.0, 500.0, 50.0, 200.0)
    assert doubly == pytest.approx(singly, rel=1e-12)
    assert doubly == pytest.approx(flex.limit_moment(1.0, 20.0, 300.0, 200.0, 500.0), rel=1e-12)


# ---------------------------------------------------------------------------
# Stirrup spacing limits — §9.2.2(6) and (8)
# ---------------------------------------------------------------------------


def test_max_stirrup_spacing_for_vertical_stirrups():
    # alpha = 90 deg: cot(alpha) = 0, so s_max_l = 0.75*d and s_max_w = 0.75*d.
    # d = 450 mm -> 337.5 mm, both under their caps (400 and 600).
    s_l, s_w = eq.max_stirrup_spacing(450.0, math.pi / 2)
    assert s_l == pytest.approx(337.5, rel=1e-9)
    assert s_w == pytest.approx(337.5, rel=1e-9)


def test_max_stirrup_spacing_is_capped_for_a_deep_member():
    # d = 1200 mm: 0.75*d = 900, above both caps.
    s_l, s_w = eq.max_stirrup_spacing(1200.0, math.pi / 2)
    assert s_l == pytest.approx(400.0)
    assert s_w == pytest.approx(600.0)


def test_max_stirrup_spacing_grows_with_inclined_stirrups():
    """Inclined stirrups cross the crack over a longer run, so the clause allows
    them further apart along the member — until the 400 mm cap bites."""
    vertical, _ = eq.max_stirrup_spacing(300.0, math.pi / 2)
    inclined, _ = eq.max_stirrup_spacing(300.0, math.radians(45))
    assert inclined > vertical
    # 0.75*300*(1+1) = 450 -> capped at 400
    assert inclined == pytest.approx(400.0)
