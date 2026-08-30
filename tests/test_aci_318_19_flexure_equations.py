"""ACI 318-19 flexure equations, checked against hand-computable values.

Same contract as the shear equation tests: every case can be verified with the
printed code and a calculator. Units are N, mm, MPa (or lb, in, psi).
"""

import math

import pytest

from mento.codes.aci_318_19.equations import flexure as eq

# A recurring section: f_c 25 MPa, f_y 420 MPa, 300x500 with d = 450 mm.
F_C = 25.0
F_Y = 420.0
E_S = 200_000.0
EPS_Y = F_Y / E_S  # 0.0021
EPS_C = 0.003
BETA_1 = 0.85
B = 300.0
D = 450.0


# ---------------------------------------------------------------------------
# rho_max — §21.2.2
# ---------------------------------------------------------------------------


def test_max_reinforcement_ratio():
    # 0.85*0.85*25/420 * (0.003/(0.0021+0.006))
    expected = 0.85 * BETA_1 * F_C / F_Y * (EPS_C / (EPS_Y + 2 * EPS_C))
    assert eq.max_reinforcement_ratio(F_C, F_Y, BETA_1, EPS_C, EPS_Y) == pytest.approx(expected, rel=1e-12)
    # Sanity: 0.85*0.85*25/420 = 0.043006, times 0.003/0.0081 = 0.370370.
    assert eq.max_reinforcement_ratio(F_C, F_Y, BETA_1, EPS_C, EPS_Y) == pytest.approx(0.0159281, rel=1e-5)


def test_max_reinforcement_ratio_falls_with_stronger_steel():
    # Higher f_y means the steel yields at a larger strain, so the section turns
    # brittle sooner and less steel is allowed.
    weak = eq.max_reinforcement_ratio(F_C, 280.0, BETA_1, EPS_C, 280.0 / E_S)
    strong = eq.max_reinforcement_ratio(F_C, 500.0, BETA_1, EPS_C, 500.0 / E_S)
    assert strong < weak


# ---------------------------------------------------------------------------
# rho_min — §9.6.1.2
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "f_c, f_y, expected, governing",
    [
        # 0.25*sqrt(25)/420 = 2.976e-3 vs 1.4/420 = 3.333e-3 -> floor
        (25.0, 420.0, 1.4 / 420.0, "floor"),
        # 0.25*sqrt(64)/420 = 4.762e-3 vs 3.333e-3 -> sqrt term
        (64.0, 420.0, 0.25 * 8.0 / 420.0, "sqrt"),
    ],
)
def test_min_reinforcement_ratio_si(f_c, f_y, expected, governing):
    assert eq.min_reinforcement_ratio(f_c, f_y) == pytest.approx(expected, rel=1e-12)


def test_min_reinforcement_ratio_us_matches_the_validated_case():
    # Test_Etabs_05 in the beam suite: f_c = 6000 psi, f_y = 60 ksi.
    # 3*sqrt(6000)/60000 = 0.003873, which governs over 200/60000 = 0.003333.
    got = eq.min_reinforcement_ratio(6000.0, 60_000.0, is_imperial=True)
    assert got == pytest.approx(0.003873, rel=1e-4)
    assert got == pytest.approx(3 * math.sqrt(6000.0) / 60_000.0, rel=1e-12)


def test_min_reinforcement_ratio_us_floor_governs_at_low_fc():
    # f_c = 3000 psi: 3*sqrt(3000)/60000 = 0.002739 < 200/60000 = 0.003333
    assert eq.min_reinforcement_ratio(3000.0, 60_000.0, is_imperial=True) == pytest.approx(200 / 60_000.0, rel=1e-12)


# ---------------------------------------------------------------------------
# Ductility limit — §21.2.2 and §22.2
# ---------------------------------------------------------------------------


def test_neutral_axis_at_ductility_limit():
    # 0.003*450/(0.0021+0.006) = 1.35/0.0081 = 166.667 mm
    assert eq.neutral_axis_at_ductility_limit(D, EPS_Y) == pytest.approx(166.6667, rel=1e-6)


def test_neutral_axis_at_ductility_limit_is_proportional_to_d():
    assert eq.neutral_axis_at_ductility_limit(900.0, EPS_Y) == pytest.approx(
        2 * eq.neutral_axis_at_ductility_limit(450.0, EPS_Y)
    )


def test_compression_steel_net_stress_when_the_bar_yields():
    # c_t = 166.667, d' = 50: eps_s' = 0.003*(1-0.3) = 0.0021 -> f_s' = 420 (capped
    # at f_y). Net of the displaced concrete: 420 - 0.85*25 = 398.75 MPa.
    c_t = eq.neutral_axis_at_ductility_limit(D, EPS_Y)
    assert eq.compression_steel_net_stress(50.0, c_t, E_S, F_Y, F_C) == pytest.approx(398.75, rel=1e-6)


def test_compression_steel_net_stress_when_the_bar_stays_elastic():
    # d' = 120 of c_t = 166.667: eps_s' = 0.003*(1-0.72) = 8.4e-4
    # f_s' = 8.4e-4*200000 = 168 MPa, below f_y, so no cap. 168 - 21.25 = 146.75
    c_t = eq.neutral_axis_at_ductility_limit(D, EPS_Y)
    assert eq.compression_steel_net_stress(120.0, c_t, E_S, F_Y, F_C) == pytest.approx(146.75, rel=1e-6)


def test_compression_steel_net_stress_is_capped_at_yield():
    c_t = eq.neutral_axis_at_ductility_limit(D, EPS_Y)
    deep = eq.compression_steel_net_stress(10.0, c_t, E_S, F_Y, F_C)
    assert deep == pytest.approx(F_Y - 0.85 * F_C, rel=1e-12)


# ---------------------------------------------------------------------------
# Required steel — §22.2/§22.3
# ---------------------------------------------------------------------------


def test_flexural_resistance_factor():
    # M_u = 120 kNm = 120e6 N*mm, phi = 0.9: 120e6/(0.9*300*450**2) = 2.1948 MPa
    assert eq.flexural_resistance_factor(120e6, 0.9, B, D) == pytest.approx(2.194787, rel=1e-6)


def test_singly_reinforced_discriminant():
    R_n = eq.flexural_resistance_factor(120e6, 0.9, B, D)
    # 1 - 2*2.1948/(0.85*25) = 1 - 4.3896/21.25 = 0.79343
    assert eq.singly_reinforced_discriminant(R_n, F_C) == pytest.approx(0.793433, rel=1e-5)


def test_singly_reinforced_discriminant_goes_negative_when_overloaded():
    # A moment the section cannot carry with tension steel alone.
    R_n = eq.flexural_resistance_factor(900e6, 0.9, B, D)
    assert eq.singly_reinforced_discriminant(R_n, F_C) < 0


def test_tension_steel_for_moment():
    # 0.85*25*300*450/420 * (1 - sqrt(0.793433)) = 6830.36 * 0.109252 = 746.2 mm²
    R_n = eq.flexural_resistance_factor(120e6, 0.9, B, D)
    assert eq.tension_steel_for_moment(R_n, F_C, F_Y, B, D) == pytest.approx(746.2, rel=1e-3)


def test_tension_steel_round_trips_through_the_nominal_moment():
    """The strongest check available: solve for A_s, then put it back.

    phi*M_n computed from the resulting steel must return the demand.
    """
    M_u = 120e6
    R_n = eq.flexural_resistance_factor(M_u, 0.9, B, D)
    A_s = eq.tension_steel_for_moment(R_n, F_C, F_Y, B, D)
    M_n = eq.nominal_moment_singly_reinforced(A_s, F_Y, F_C, B, D)
    assert 0.9 * M_n == pytest.approx(M_u, rel=1e-9)


def test_neutral_axis_depth():
    # c = 746.2*420/(0.85*25*300*0.85) = 313404/5418.75 = 57.84 mm
    assert eq.neutral_axis_depth(746.2, F_Y, F_C, B, BETA_1) == pytest.approx(57.84, rel=1e-3)


def test_neutral_axis_depth_grows_with_steel():
    assert eq.neutral_axis_depth(1500.0, F_Y, F_C, B, BETA_1) == pytest.approx(
        2 * eq.neutral_axis_depth(750.0, F_Y, F_C, B, BETA_1)
    )


# ---------------------------------------------------------------------------
# Nominal moment — §22.3
# ---------------------------------------------------------------------------


def test_nominal_moment_singly_reinforced():
    # a = 1000*420/(0.85*25*300) = 65.882 mm
    # M_n = 1000*420*(450 - 32.941) = 175.16e6 N*mm
    assert eq.nominal_moment_singly_reinforced(1000.0, F_Y, F_C, B, D) == pytest.approx(175.165e6, rel=1e-4)


def test_doubly_reinforced_approaches_singly_when_compression_steel_vanishes():
    """With A_s' -> 0 the two formulas must agree: no compression bar, no
    difference in how the section resolves."""
    singly = eq.nominal_moment_singly_reinforced(1000.0, F_Y, F_C, B, D)
    doubly = eq.nominal_moment_doubly_reinforced(1000.0, 1e-9, F_Y, F_C, B, D, 50.0, BETA_1, EPS_C, EPS_Y, E_S)
    assert doubly == pytest.approx(singly, rel=1e-6)


def test_doubly_reinforced_yielding_branch():
    """A_s' shallow enough to yield: hand-check the closed form."""
    A_s, A_s_prime, d_prime = 4000.0, 1000.0, 50.0
    c = (A_s * F_Y - A_s_prime * (F_Y - 0.85 * F_C)) / (0.85 * F_C * B * BETA_1)
    eps_s = (c - d_prime) / c * EPS_C
    assert eps_s >= EPS_Y, "this case is meant to exercise the yielding branch"

    a = c * BETA_1
    expected = 0.85 * F_C * a * B * (D - a / 2) + A_s_prime * (F_Y - 0.85 * F_C) * (D - d_prime)
    got = eq.nominal_moment_doubly_reinforced(A_s, A_s_prime, F_Y, F_C, B, D, d_prime, BETA_1, EPS_C, EPS_Y, E_S)
    assert got == pytest.approx(expected, rel=1e-12)


def test_doubly_reinforced_elastic_branch_solves_the_quadratic():
    """A_s' deep enough that it cannot yield: the quadratic branch must run and
    the neutral axis it finds must satisfy the equilibrium it was derived from."""
    A_s, A_s_prime, d_prime = 2000.0, 1500.0, 150.0
    c_assumed = (A_s * F_Y - A_s_prime * (F_Y - 0.85 * F_C)) / (0.85 * F_C * B * BETA_1)
    eps_s = (c_assumed - d_prime) / c_assumed * EPS_C
    assert eps_s < EPS_Y, "this case is meant to exercise the elastic branch"

    A = 0.85 * F_C * B * BETA_1
    Bq = A_s_prime * (EPS_C * E_S - 0.85 * F_C) - A_s * F_Y
    C = -d_prime * A_s_prime * EPS_C * E_S
    c = (-Bq + math.sqrt(Bq**2 - 4 * A * C)) / (2 * A)
    assert A * c**2 + Bq * c + C == pytest.approx(0.0, abs=1e-3)

    f_s_net = (c - d_prime) / c * EPS_C * E_S - 0.85 * F_C
    a = c * BETA_1
    expected = (A_s - A_s_prime * f_s_net / F_Y) * F_Y * (D - a / 2) + A_s_prime * f_s_net * (D - d_prime)
    got = eq.nominal_moment_doubly_reinforced(A_s, A_s_prime, F_Y, F_C, B, D, d_prime, BETA_1, EPS_C, EPS_Y, E_S)
    assert got == pytest.approx(expected, rel=1e-12)


def test_doubly_reinforced_handles_equal_top_and_bottom_steel():
    """A_s == A_s' drives the assumed neutral axis to zero; the guard in the
    equation must send it to the elastic branch instead of dividing by zero."""
    got = eq.nominal_moment_doubly_reinforced(2000.0, 2000.0, F_Y, F_C, B, D, 50.0, BETA_1, EPS_C, EPS_Y, E_S)
    assert got > 0
