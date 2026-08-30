"""ACI 318-19 shear equations, checked against values read off the code itself.

These are the Phase 1 equation tests: each case is a hand-computable evaluation
of the printed clause, so a reader with the code open can verify a row without
running mento. Nothing here imports pint or touches a beam.
"""

import math

import pytest

from mento.codes.aci_318_19.equations import shear as eq


# ---------------------------------------------------------------------------
# lambda_s — Eq. 22.5.5.1.3
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "d_mm, expected",
    [
        # d = 250 mm is the transition: 2/(1+0.004*250) = 1, so lambda_s = 1.
        (250.0, 1.0),
        (100.0, 1.0),  # shallower than the transition -> capped, not > 1
        (500.0, math.sqrt(2 / 3)),  # 2/(1+2.0)
        (1000.0, math.sqrt(0.4)),  # 2/(1+4.0)
    ],
)
def test_size_effect_factor_si(d_mm, expected):
    assert eq.size_effect_factor(d_mm) == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize(
    "d_in, expected",
    [
        (10.0, 1.0),  # transition, 2/(1+1)
        (4.0, 1.0),  # capped
        (20.0, math.sqrt(2 / 3)),
        (40.0, math.sqrt(0.4)),
    ],
)
def test_size_effect_factor_us(d_in, expected):
    assert eq.size_effect_factor(d_in, is_imperial=True) == pytest.approx(expected, rel=1e-9)


def test_size_effect_factor_never_exceeds_one():
    # The cap is the point of the clause: it may only reduce V_c.
    assert all(eq.size_effect_factor(d) <= 1.0 for d in range(10, 2000, 10))


# ---------------------------------------------------------------------------
# sigma_Nu — §22.5.5.1
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "N_u, A_g, f_c, expected",
    [
        (0.0, 100_000.0, 25.0, 0.0),  # no axial load
        (600_000.0, 100_000.0, 25.0, 1.0),  # 600e3/(6*100e3) = 1.0 MPa, under the cap
        (6_000_000.0, 100_000.0, 25.0, 1.25),  # would be 10 MPa; capped at 0.05*f_c
    ],
)
def test_axial_stress_influence(N_u, A_g, f_c, expected):
    assert eq.axial_stress_influence(N_u, A_g, f_c) == pytest.approx(expected, rel=1e-9)


# ---------------------------------------------------------------------------
# v_c — Table 22.5.5.1
# ---------------------------------------------------------------------------


def test_concrete_shear_stress_with_min_rebar_si():
    # f_c = 25 MPa, lambda = 1, rho_w = 0.01, sigma_Nu = 0.
    #   row (a): 0.17*sqrt(25)          = 0.85
    #   row (b): 0.66*0.01**(1/3)*sqrt(25) = 0.7109...
    # The member may take the larger.
    got = eq.concrete_shear_stress(25.0, 1.0, 0.01, 0.0, 1.0, has_min_rebar=True)
    assert got == pytest.approx(0.85, rel=1e-9)


def test_concrete_shear_stress_without_min_rebar_si():
    # Row (c): 0.66*lambda_s*lambda*rho_w**(1/3)*sqrt(f_c)
    expected = 0.66 * 1.0 * 1.0 * 0.01 ** (1 / 3) * 5.0
    got = eq.concrete_shear_stress(25.0, 1.0, 0.01, 0.0, 1.0, has_min_rebar=False)
    assert got == pytest.approx(expected, rel=1e-9)


def test_concrete_shear_stress_row_c_scales_with_size_effect():
    # Without stirrups the size effect factor is what reduces v_c in a deep
    # member; with stirrups it does not appear at all.
    without = eq.concrete_shear_stress(25.0, 1.0, 0.01, 0.0, 0.8, has_min_rebar=False)
    full = eq.concrete_shear_stress(25.0, 1.0, 0.01, 0.0, 1.0, has_min_rebar=False)
    assert without == pytest.approx(0.8 * full, rel=1e-9)

    with_rebar_small = eq.concrete_shear_stress(25.0, 1.0, 0.01, 0.0, 0.5, has_min_rebar=True)
    with_rebar_full = eq.concrete_shear_stress(25.0, 1.0, 0.01, 0.0, 1.0, has_min_rebar=True)
    assert with_rebar_small == with_rebar_full


def test_concrete_shear_stress_us():
    # f_c = 4000 psi, rho_w = 0.01: row (a) 2*sqrt(4000) = 126.49 psi governs
    # over row (b) 8*0.01**(1/3)*sqrt(4000) = 109.0 psi.
    got = eq.concrete_shear_stress(4000.0, 1.0, 0.01, 0.0, 1.0, has_min_rebar=True, is_imperial=True)
    assert got == pytest.approx(2 * math.sqrt(4000.0), rel=1e-9)


def test_concrete_shear_stress_adds_axial_term():
    base = eq.concrete_shear_stress(25.0, 1.0, 0.01, 0.0, 1.0, has_min_rebar=True)
    with_axial = eq.concrete_shear_stress(25.0, 1.0, 0.01, 0.4, 1.0, has_min_rebar=True)
    assert with_axial == pytest.approx(base + 0.4, rel=1e-9)


# ---------------------------------------------------------------------------
# Stress limits — §22.5.5.1.1, §22.5.1.2, §9.6.3.1
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fn, f_c, is_imperial, expected",
    [
        (eq.max_concrete_shear_stress, 25.0, False, 0.42 * 5.0),
        (eq.max_concrete_shear_stress, 4000.0, True, 5 * math.sqrt(4000.0)),
        (eq.shear_stress_capacity_increment, 25.0, False, 0.66 * 5.0),
        (eq.shear_stress_capacity_increment, 4000.0, True, 8 * math.sqrt(4000.0)),
        (eq.min_shear_reinforcement_threshold_stress, 25.0, False, 0.083 * 5.0),
        (eq.min_shear_reinforcement_threshold_stress, 4000.0, True, math.sqrt(4000.0)),
    ],
)
def test_stress_limits(fn, f_c, is_imperial, expected):
    assert fn(f_c, 1.0, is_imperial=is_imperial) == pytest.approx(expected, rel=1e-9)


def test_stress_limits_scale_with_lambda():
    # Lightweight concrete reduces every sqrt(f_c) term linearly.
    assert eq.max_concrete_shear_stress(25.0, 0.75) == pytest.approx(0.75 * eq.max_concrete_shear_stress(25.0, 1.0))


# ---------------------------------------------------------------------------
# f_yt — §20.2.2.4
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("f_y, expected", [(280.0, 280.0), (420.0, 420.0), (500.0, 420.0), (600.0, 420.0)])
def test_max_yield_strength_for_shear_si(f_y, expected):
    assert eq.max_yield_strength_for_shear(f_y) == pytest.approx(expected)


@pytest.mark.parametrize("f_y, expected", [(40_000.0, 40_000.0), (60_000.0, 60_000.0), (75_000.0, 60_000.0)])
def test_max_yield_strength_for_shear_us(f_y, expected):
    assert eq.max_yield_strength_for_shear(f_y, is_imperial=True) == pytest.approx(expected)


# ---------------------------------------------------------------------------
# A_v,min — Table 9.6.3.4
# ---------------------------------------------------------------------------


def test_min_shear_reinforcement_ratio_floor_governs():
    # f_c = 25 MPa, f_yt = 420 MPa, b_w = 300 mm.
    #   0.062*sqrt(25)/420 = 7.381e-4   <-- sqrt term
    #   0.35/420           = 8.333e-4   <-- floor governs at low f_c
    # 8.333e-4 * 300 = 0.25 mm²/mm = 2.5 cm²/m, the familiar minimum.
    assert eq.min_shear_reinforcement_ratio(25.0, 420.0, 300.0) == pytest.approx(0.25, rel=1e-6)


def test_min_shear_reinforcement_ratio_sqrt_term_governs_at_high_fc():
    # At f_c = 64 MPa the sqrt term overtakes the floor: 0.062*8/420 = 1.181e-3.
    expected = 0.062 * 8.0 / 420.0 * 300.0
    assert eq.min_shear_reinforcement_ratio(64.0, 420.0, 300.0) == pytest.approx(expected, rel=1e-9)
    assert eq.min_shear_reinforcement_ratio(64.0, 420.0, 300.0) > eq.min_shear_reinforcement_ratio(25.0, 420.0, 300.0)


def test_min_shear_reinforcement_ratio_us():
    # f_c = 4000 psi, f_yt = 60 ksi, b_w = 12 in:
    #   0.75*sqrt(4000)/60000 = 7.906e-4 ; 50/60000 = 8.333e-4 -> floor governs
    assert eq.min_shear_reinforcement_ratio(4000.0, 60_000.0, 12.0, is_imperial=True) == pytest.approx(0.01, rel=1e-6)


def test_min_shear_reinforcement_ratio_scales_with_width():
    single = eq.min_shear_reinforcement_ratio(25.0, 420.0, 300.0)
    double = eq.min_shear_reinforcement_ratio(25.0, 420.0, 600.0)
    assert double == pytest.approx(2 * single, rel=1e-12)


# ---------------------------------------------------------------------------
# V_s — §22.5.8.5.3
# ---------------------------------------------------------------------------


def test_shear_strength_of_reinforcement():
    # A_v/s = 0.25 mm²/mm, f_yt = 420 MPa, d = 450 mm -> 47.25 kN
    assert eq.shear_strength_of_reinforcement(0.25, 420.0, 450.0) == pytest.approx(47_250.0, rel=1e-9)


def test_shear_strength_of_reinforcement_is_linear():
    base = eq.shear_strength_of_reinforcement(0.25, 420.0, 450.0)
    assert eq.shear_strength_of_reinforcement(0.5, 420.0, 450.0) == pytest.approx(2 * base)
    assert eq.shear_strength_of_reinforcement(0.25, 420.0, 900.0) == pytest.approx(2 * base)
