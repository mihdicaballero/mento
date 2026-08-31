"""ACI 318-19 structural wall equations, checked against the printed clauses."""

import math

import pytest

from mento.codes.aci_318_19.equations import wall as eq


# ---------------------------------------------------------------------------
# alpha_c — §11.5.4.6
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "hw_lw, expected",
    [
        (0.5, 0.25),  # squat wall, flat below the transition
        (1.5, 0.25),  # start of the transition
        (1.75, 0.21),  # midway: 0.25 + 0.5*(0.17-0.25)
        (2.0, 0.17),  # end of the transition
        (4.0, 0.17),  # slender wall, flat above
    ],
)
def test_alpha_c_si(hw_lw, expected):
    assert eq.alpha_c(hw_lw) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize(
    "hw_lw, expected",
    [(1.0, 3.0), (1.5, 3.0), (1.75, 2.5), (2.0, 2.0), (3.0, 2.0)],
)
def test_alpha_c_us(hw_lw, expected):
    assert eq.alpha_c(hw_lw, is_imperial=True) == pytest.approx(expected, rel=1e-12)


def test_alpha_c_never_increases_with_slenderness():
    ratios = [i / 10 for i in range(5, 41)]
    values = [eq.alpha_c(r) for r in ratios]
    assert all(a >= b for a, b in zip(values, values[1:]))


# ---------------------------------------------------------------------------
# Shear stresses — §11.5.4.3 and §11.5.4.6
# ---------------------------------------------------------------------------


def test_concrete_shear_stress():
    # alpha_c = 0.17, lambda = 1, f_c = 25 -> 0.17*sqrt(25) = 0.85 MPa
    assert eq.concrete_shear_stress(25.0, 0.17, 1.0) == pytest.approx(0.85, rel=1e-12)


def test_concrete_shear_stress_uses_alpha_c_for_the_unit_system():
    """The system-dependent constant lives in alpha_c, so this expression is
    shared: feeding it the US alpha_c and psi must give the US answer."""
    got = eq.concrete_shear_stress(4000.0, 2.0, 1.0)
    assert got == pytest.approx(2.0 * math.sqrt(4000.0), rel=1e-12)


@pytest.mark.parametrize(
    "f_c, is_imperial, expected",
    [(25.0, False, 0.66 * 5.0), (4000.0, True, 8.0 * math.sqrt(4000.0))],
)
def test_max_shear_stress(f_c, is_imperial, expected):
    assert eq.max_shear_stress(f_c, 1.0, is_imperial=is_imperial) == pytest.approx(expected, rel=1e-12)


def test_reinforcement_shear_stress():
    # rho_t = 0.0025 (the code minimum), f_yt = 420 -> 1.05 MPa
    assert eq.reinforcement_shear_stress(0.0025, 420.0) == pytest.approx(1.05, rel=1e-12)


# ---------------------------------------------------------------------------
# rho_l,min — Eq. (11.6.2)
# ---------------------------------------------------------------------------


def test_min_vertical_ratio_equals_horizontal_for_a_squat_wall():
    # hw/lw <= 0.5: 0.0025 + 0.5*2.0*(rho_t - 0.0025) = rho_t
    assert eq.min_vertical_reinforcement_ratio(0.5, 0.005) == pytest.approx(0.005, rel=1e-12)


def test_min_vertical_ratio_clamps_below_half():
    assert eq.min_vertical_reinforcement_ratio(0.1, 0.005) == eq.min_vertical_reinforcement_ratio(0.5, 0.005)


def test_min_vertical_ratio_drops_to_the_floor_for_a_slender_wall():
    # hw/lw >= 2.5 leaves only the 0.0025 minimum.
    assert eq.min_vertical_reinforcement_ratio(2.5, 0.005) == pytest.approx(0.0025, rel=1e-12)
    assert eq.min_vertical_reinforcement_ratio(6.0, 0.005) == pytest.approx(0.0025, rel=1e-12)


def test_min_vertical_ratio_interpolates_in_between():
    # hw/lw = 1.5: 0.0025 + 0.5*1.0*(0.005-0.0025) = 0.00375
    assert eq.min_vertical_reinforcement_ratio(1.5, 0.005) == pytest.approx(0.00375, rel=1e-12)


def test_min_vertical_ratio_never_falls_below_the_floor():
    # A rho_t below the minimum would drive the equation under 0.0025.
    assert eq.min_vertical_reinforcement_ratio(0.5, 0.001) == pytest.approx(eq.MIN_REINFORCEMENT_RATIO, rel=1e-12)


# ---------------------------------------------------------------------------
# Spacing limits — §11.7.3
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "l_w, t, expected_h, expected_v",
    [
        # lw/5 = 600, 3t = 600, cap 450 -> the absolute cap governs both
        (3000.0, 200.0, 450.0, 450.0),
        # lw/5 = 200 governs horizontally; 3t = 300 governs vertically
        (1000.0, 100.0, 200.0, 300.0),
    ],
)
def test_spacing_limits_si(l_w, t, expected_h, expected_v):
    assert eq.max_horizontal_spacing(l_w, t) == pytest.approx(expected_h, rel=1e-12)
    assert eq.max_vertical_spacing(l_w, t) == pytest.approx(expected_v, rel=1e-12)


def test_spacing_limits_us():
    # lw = 120 in, t = 8 in: lw/5 = 24, 3t = 24, cap 18 -> cap governs
    assert eq.max_horizontal_spacing(120.0, 8.0, is_imperial=True) == pytest.approx(18.0)
    assert eq.max_vertical_spacing(120.0, 8.0, is_imperial=True) == pytest.approx(18.0)


def test_vertical_spacing_is_never_tighter_than_horizontal():
    """The horizontal bars resist in-plane shear, so their limit is the
    stricter of the two for every geometry."""
    for l_w in (500.0, 1000.0, 2000.0, 5000.0):
        for t in (100.0, 150.0, 300.0):
            assert eq.max_horizontal_spacing(l_w, t) <= eq.max_vertical_spacing(l_w, t)
