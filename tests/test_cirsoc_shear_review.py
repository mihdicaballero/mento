"""Regression cases for the ACI/CIRSOC shear review.

The expected values come from the cited clauses or invariance under a change
of load sign/order, rather than from copying the implementation's outputs.
"""

import math
import warnings

import pytest

from mento import (
    Concrete_ACI_318_19,
    Concrete_CIRSOC_201_25,
    Concrete_EN_1992_2004,
    Footing,
    Forces,
    OneWaySlab,
    RectangularBeam,
    ShearWall,
    SteelBar,
)
from mento.codes.registry import design_code
from mento.codes.ACI_318_19_beam import _calculate_max_shear_capacity_aci
from mento.codes.aci_318_19.equations.flexure import min_reinforcement_ratio
from mento.units import MPa, cm, inch, kN, kNm, mm, psi


ACI_CODES = [Concrete_ACI_318_19, Concrete_CIRSOC_201_25]


def beam(concrete_cls=Concrete_ACI_318_19, fc=25, *, height=60, stirrups=True):
    section = RectangularBeam(
        concrete=concrete_cls(name="Concrete", f_c=fc * MPa),
        steel_bar=SteelBar(name="Steel", f_y=420 * MPa),
        width=30 * cm,
        height=height * cm,
        c_c=25 * mm,
    )
    section.set_longitudinal_rebar_bot(n1=3, d_b1=20 * mm)
    section.set_longitudinal_rebar_top(n1=2, d_b1=12 * mm)
    if stirrups:
        section.set_transverse_rebar(n_stirrups=1, d_b=10 * mm, s_l=20 * cm)
    else:
        section.set_transverse_rebar(n_stirrups=0, d_b=0 * mm, s_l=0 * mm)
    return section


def state(section, vu, *, nu=0, moment=0):
    code = design_code(section.concrete)
    check = code.check_shear_wall if isinstance(section, ShearWall) else code.check_shear
    return check(section, Forces(V_z=vu * kN, N_x=nu * kN, M_y=moment * kNm))


@pytest.mark.parametrize("concrete_cls", ACI_CODES)
@pytest.mark.parametrize("fc", [25, 80])
def test_provided_minimum_controls_vc_independently_of_demand(concrete_cls, fc):
    section = beam(concrete_cls, fc)
    checks = [state(section, vu) for vu in (0, 20, 250)]
    rho = checks[0].rho_w
    # Table 22.5.5.1(a)/(b); §22.5.3.2 permits the full root in this beam.
    expected = max(0.17, 0.66 * rho ** (1 / 3)) * math.sqrt(fc)
    assert [s.k_c_min for s in checks] == pytest.approx([expected] * 3)
    assert checks[0].A_v_min == 0  # waived demand minimum, not missing stirrups
    assert checks[-1].A_v_min > 0
    assert checks[0].V_s_req == 0  # required nominal strength cannot be negative


def slab(concrete_cls=Concrete_ACI_318_19, fc=25):
    section = OneWaySlab(
        concrete=concrete_cls(name="Concrete", f_c=fc * MPa),
        steel_bar=SteelBar(name="Steel", f_y=420 * MPa),
        width=100 * cm,
        height=25 * cm,
        c_c=25 * mm,
    )
    section.set_slab_longitudinal_rebar_bot(d_b1=16 * mm, s_b1=15 * cm)
    return section


@pytest.mark.parametrize("concrete_cls", ACI_CODES)
def test_slab_requires_table_minimum_just_above_concrete_capacity(concrete_cls):
    section = slab(concrete_cls)
    bare = state(section, 0)
    # §7.6.3.1/§7.6.3.3: no minimum at equality, table minimum just above.
    assert state(section, bare.phi_V_c / 1000).A_v_min == 0
    loaded = state(section, bare.phi_V_c / 1000 + 0.1)
    expected = 0.35 * 1000 / 420  # mm²/mm; Table 9.6.3.4, H25
    assert loaded.A_v_min == pytest.approx(expected)
    assert loaded.A_v_req == pytest.approx(expected)
    section.design_shear([Forces(V_z=(bare.phi_V_c / 1000 + 0.1) * kN)])
    assert section.shear_design.A_v.to("mm**2/mm").magnitude >= expected


@pytest.mark.parametrize("concrete_cls", ACI_CODES)
def test_reinforced_slab_uses_table_rows_a_b_but_keeps_root_cap(concrete_cls):
    section = slab(concrete_cls, fc=80)
    section.set_transverse_rebar(n_stirrups=4, d_b=10 * mm, s_l=10 * cm)
    result = state(section, 20)
    # §22.5.3.2 exempts beams/joists, not slabs. Table row still depends on Av.
    assert result.k_c_min == pytest.approx(max(0.17, 0.66 * result.rho_w ** (1 / 3)) * 8.3)


def test_beam_minimum_threshold_is_specific_to_each_code():
    aci, cirsoc = (beam(cls, stirrups=False) for cls in ACI_CODES)
    area = state(aci, 0).A_cv
    demand = 0.75 * 0.084 * math.sqrt(25) * area / 1000
    assert state(aci, demand).A_v_min > 0
    assert state(cirsoc, demand).A_v_min == 0
    for section, coefficient in ((aci, 0.083), (cirsoc, 0.085)):
        threshold = 0.75 * coefficient * math.sqrt(25) * area / 1000
        assert state(section, threshold * (1 - 1e-9)).A_v_min == 0
        assert state(section, threshold * (1 + 1e-9)).A_v_min > 0


@pytest.mark.parametrize("concrete_cls", ACI_CODES)
def test_shallow_beam_minimum_exemption_does_not_remove_design_cage(concrete_cls):
    section = beam(concrete_cls, height=20, stirrups=False)
    bare = state(section, 0)
    demand = 0.9 * bare.phi_V_c / 1000
    assert state(section, demand).A_v_min == 0
    assert state(section, bare.phi_V_c / 1000 + 0.1).A_v_min > 0
    section.design_shear([Forces(V_z=demand * kN)])
    assert section.shear_design.A_v.to("mm**2/mm").magnitude >= 0.35 * 300 / 420


@pytest.mark.parametrize("concrete_cls", [*ACI_CODES, Concrete_EN_1992_2004])
def test_reversing_shear_preserves_check_and_design(concrete_cls):
    positive, negative = beam(concrete_cls), beam(concrete_cls)
    checks = [state(positive, 300), state(negative, -300)]
    assert checks[0].A_v_req == pytest.approx(checks[1].A_v_req)
    assert checks[0].DCR == pytest.approx(checks[1].DCR)
    assert checks[0].max_shear_ok == checks[1].max_shear_ok
    for section, sign in ((positive, 1), (negative, -1)):
        section.design_shear([Forces(V_z=sign * 300 * kN)])
    assert positive.reinforcement.transverse == negative.reinforcement.transverse


def wall(concrete_cls=Concrete_ACI_318_19):
    return ShearWall(
        concrete=concrete_cls(name="Concrete", f_c=25 * MPa),
        steel_bar=SteelBar(name="Steel", f_y=420 * MPa),
        thickness=25 * cm,
        length=400 * cm,
        height=350 * cm,
        c_c=25 * mm,
    )


@pytest.mark.parametrize("concrete_cls, divisor", [(Concrete_ACI_318_19, 3.45), (Concrete_CIRSOC_201_25, 3.5)])
@pytest.mark.parametrize("fraction", [0.5, 1.0, 1.5])
def test_wall_net_axial_tension_reduces_concrete_shear(concrete_cls, divisor, fraction):
    section = wall(concrete_cls)
    # A_g = 1 m². N_u = -fraction * divisor MPa * A_g.
    result = state(section, 1200, nu=-fraction * divisor * 1000)
    alpha = max(0.17 * (1 - fraction), 0)
    assert result.alpha_c == pytest.approx(alpha)
    assert result.V_c_wall.to("kN").magnitude == pytest.approx(alpha * 5 * 1000)
    assert result.rho_t_req.magnitude == pytest.approx(max((1200 / 0.75 / 1000 - alpha * 5) / 420, 0.0025))


def test_wall_tension_uses_imperial_clause():
    section = ShearWall(
        concrete=Concrete_ACI_318_19(name="Concrete", f_c=4000 * psi),
        steel_bar=SteelBar(name="Steel", f_y=60000 * psi),
        thickness=10 * inch,
        length=160 * inch,
        height=140 * inch,
        c_c=1 * inch,
    )
    force = Forces(V_z=1200 * kN, N_x=-250 * psi * (1600 * inch**2))
    result = design_code(section.concrete).check_shear_wall(section, force)
    assert result.alpha_c == pytest.approx(1.0)  # 2*(1-250/500)


@pytest.mark.parametrize("concrete_cls", ACI_CODES)
def test_wall_design_envelope_includes_tension_and_both_shear_signs(concrete_cls):
    section = wall(concrete_cls)
    compression = state(section, 1200, nu=500)
    tension = state(section, -1200, nu=-1750)
    assert compression.alpha_c == pytest.approx(0.25)
    assert tension.rho_t_req > compression.rho_t_req
    assert tension.DCR >= 0
    design_code(section.concrete).design_shear_wall(
        section, [Forces(V_z=-1200 * kN, N_x=-1750 * kN), Forces(V_z=1200 * kN, N_x=500 * kN)]
    )
    assert section._rho_t >= tension.rho_t_req


def test_en_shear_uses_current_tension_face_without_flexure_history():
    section = beam(Concrete_EN_1992_2004, stirrups=False)
    section.set_longitudinal_rebar_top(n1=4, d_b1=25 * mm)
    first = state(section, 60, moment=-20)
    d = section._d_shear.to(mm).magnitude
    rho = min((4 * math.pi * 25**2 / 4) / (300 * d), 0.02)
    k = min(1 + math.sqrt(200 / d), 2)
    expected = max(0.12 * k * (100 * rho * 25) ** (1 / 3), 0.035 * k**1.5 * math.sqrt(25)) * 300 * d
    assert first.V_Rd_c == pytest.approx(expected)
    for moment in (20, -20):
        section.check_flexure([Forces(M_y=moment * kNm)])
        assert state(section, 60, moment=-20).V_Rd_c == pytest.approx(expected)
    flexural_ratio = section._rho_l_top
    section.check_shear([Forces(V_z=60 * kN, M_y=-20 * kNm)])
    assert section._rho_l_top == flexural_ratio
    assert section._shear_concrete["Value"][0] == round(rho, 4)


def test_en_shear_capacity_and_dcr_include_concrete_strut_limit():
    section = beam(Concrete_EN_1992_2004)
    section.set_transverse_rebar(n_stirrups=2, d_b=16 * mm, s_l=5 * cm)
    result = state(section, 1000)
    assert result.V_Rd_s > result.V_Rd_max
    assert not result.max_shear_ok
    assert result.V_Rd == result.V_Rd_max
    assert result.DCR == pytest.approx(1_000_000 / result.V_Rd_max)
    assert result.DCR > 1


def test_en_zero_concrete_capacity_reports_failure_without_negative_capacity():
    section = beam(Concrete_EN_1992_2004, stirrups=False)
    result = state(section, 50, nu=-5000)
    assert result.V_Rd_c == 0
    assert result.V_Rd == 0
    assert math.isinf(result.DCR)


@pytest.mark.parametrize("concrete_cls", ACI_CODES)
def test_maximum_shear_capacity_boundary_is_inclusive(concrete_cls):
    section = beam(concrete_cls)
    result = state(section, 0)
    # Exercise exact equality without rounding N -> kN -> N in Forces.
    result.V_u = result.phi_V_max
    _calculate_max_shear_capacity_aci(section, result)
    assert result.max_shear_ok
    result.V_u += 0.01
    _calculate_max_shear_capacity_aci(section, result)
    assert not result.max_shear_ok


@pytest.mark.parametrize("imperial, fc, fy, cap", [(False, 25, 600, 550), (True, 4000, 100000, 80000)])
def test_minimum_flexural_equation_keeps_two_argument_interface(imperial, fc, fy, cap):
    assert min_reinforcement_ratio(fc, fy, is_imperial=imperial) == pytest.approx(
        min_reinforcement_ratio(fc, fy, cap, is_imperial=imperial)
    )


@pytest.mark.parametrize("concrete_cls", ACI_CODES)
@pytest.mark.parametrize("check", ["check_flexure", "check_shear"])
def test_footing_warns_when_actual_bars_reduce_depth_below_minimum(concrete_cls, check):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        section = Footing(
            concrete=concrete_cls(name="Concrete", f_c=25 * MPa),
            steel_bar=SteelBar(name="Steel", f_y=420 * MPa),
            width=100 * cm,
            height=24 * cm,
            c_c=75 * mm,
        )
    assert not caught  # thinnest Ø10 mat leaves 160 mm
    section.set_slab_longitudinal_rebar_bot(d_b1=32 * mm, s_b1=20 * cm)
    assert section._d_bot.to(mm).magnitude == pytest.approx(149)
    with pytest.warns(UserWarning, match="effective depth.*150"):
        getattr(section, check)([Forces(M_y=20 * kNm)])


@pytest.mark.parametrize("concrete_cls", ACI_CODES)
def test_footing_design_checks_depth_with_selected_bars(concrete_cls):
    section = Footing(
        concrete=concrete_cls(name="Concrete", f_c=25 * MPa),
        steel_bar=SteelBar(name="Steel", f_y=420 * MPa),
        width=100 * cm,
        height=232 * mm,
        c_c=75 * mm,
    )
    with pytest.warns(UserWarning, match="effective depth.*150"):
        section.design_flexure([Forces(M_y=60 * kNm)])
    assert section._d_bot < 150 * mm
