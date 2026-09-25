"""Properties every flexure design has to have, swept over sections and moments.

The examples elsewhere pin numbers; these pin behaviour across a grid wide enough
to reach the cases a single example misses -- a section past the tension-controlled
limit, a catalogue with no bars between A_s,req and A_s,max, compression steel that
only works in one layer. Each property was broken before it was written down:

* A design passes its own check. If it found a layout, that layout carries the
  moment and keeps within the code's limits; if it did not, it says so. In a sweep
  of 960 ACI 318-19 designs (b 12-25 cm, h 25-60 cm, f'c 20-30 MPa, up to twice
  the singly reinforced limit) the code before this left 311 with DCR > 1 and no
  warning, and 195 that were not tension-controlled, also without one.
* The capacity a check reports is never more than the section has. It is compared
  here with strain compatibility written out again, independently of mento: the
  check used to cut the tension steel back to A_s,max + A_s'*f_s'/f_y and keep
  phi = 0.90, which in the same sweep passed 32 designs whose real DCR was up to
  1.07.
* The steel asked for never drops as the moment grows. It fell back to A_s,max once
  the moment passed what tension steel alone can carry.
"""

import math
from typing import Any, List, Tuple

import pytest

from mento import (
    Concrete_ACI_318_19,
    Concrete_CIRSOC_201_25,
    Concrete_EN_1992_2004,
    Forces,
    Node,
    RectangularBeam,
    SteelBar,
)
from mento.units import MPa, cm, kN, kNm, mm

_F_Y = 420.0
_E_S = 200000.0
_CONCRETES = {
    "ACI 318-19": Concrete_ACI_318_19,
    "CIRSOC 201-25": Concrete_CIRSOC_201_25,
    "EN 1992-2004": Concrete_EN_1992_2004,
}
#: Flags a design raises when it could not find a layout that works.
_CANNOT = {"As_below_required", "bars_do_not_fit"}


def _independent_capacity(A_s: float, A_sp: float, d: float, dp: float, b: float, f_c: float) -> Tuple[float, float]:
    """(phi*Mn in kN·m, eps_t) by strain compatibility, ACI 318-19 §22.2 and Table 21.2.2.

    Written apart from mento on purpose: its own bisection, every bar at the stress
    its strain gives it, the compression bar net of the concrete it displaces, and
    phi from the strain the tension steel reaches. mm, MPa.
    """
    beta_1 = min(0.85, max(0.65, 0.85 - 0.05 * (f_c - 28) / 7))

    def unbalance(c: float) -> float:
        f_sp = max(-_F_Y, min(_F_Y, _E_S * 0.003 * (c - dp) / c))
        f_s = max(-_F_Y, min(_F_Y, _E_S * 0.003 * (d - c) / c))
        return 0.85 * f_c * beta_1 * c * b + A_sp * (f_sp - 0.85 * f_c) - A_s * f_s

    lo, hi = 1e-9, d
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if unbalance(mid) > 0 else (mid, hi)
    c = (lo + hi) / 2
    a = beta_1 * c
    f_sp = max(-_F_Y, min(_F_Y, _E_S * 0.003 * (c - dp) / c))
    M_n = 0.85 * f_c * a * b * (d - a / 2) + A_sp * (f_sp - 0.85 * f_c) * (d - dp)
    eps_t = 0.003 * (d - c) / c
    phi = max(0.65, min(0.90, 0.65 + 0.25 * (eps_t - _F_Y / _E_S) / 0.003))
    return phi * M_n / 1e6, eps_t


def _moment_at_singly_limit(b_cm: float, h_cm: float, f_c: float) -> float:
    """phi*Mn of a b x d section at rho_max, kN·m: where compression steel starts."""
    d = h_cm * 10 - 45
    beta_1 = min(0.85, max(0.65, 0.85 - 0.05 * (f_c - 28) / 7))
    a = beta_1 * 0.003 * d / (_F_Y / _E_S + 0.006)
    return 0.9 * 0.85 * f_c * a * b_cm * 10 * (d - a / 2) / 1e6


def _grid() -> List[Tuple[float, float, float, float]]:
    """Narrow and wide webs, shallow and deep, singly reinforced to well past it."""
    cases = []
    for b_cm in (12, 15, 20):
        for h_cm in (25, 30, 50):
            for f_c in (20, 30):
                M_lim = _moment_at_singly_limit(b_cm, h_cm, f_c)
                for k in (0.9, 1.15, 1.6):
                    for sign in (1, -1):
                        cases.append((b_cm, h_cm, f_c, round(sign * k * M_lim, 2)))
    return cases


def _designed(code: str, b_cm: float, h_cm: float, f_c: float, M: float) -> Tuple[RectangularBeam, Node]:
    beam = RectangularBeam(
        label="P",
        concrete=_CONCRETES[code](name="C", f_c=f_c * MPa),
        steel_bar=SteelBar(name="S", f_y=_F_Y * MPa),
        width=b_cm * cm,
        height=h_cm * cm,
        c_c=25 * mm,
    )
    node = Node(section=beam, forces=[Forces(label="M", M_y=M * kNm)])
    node.design_flexure()
    node.check_flexure()
    return beam, node


@pytest.mark.parametrize("code", ["ACI 318-19", "CIRSOC 201-25", "EN 1992-2004"])
def test_a_design_passes_its_own_check(code: str) -> None:
    """Either the layout works -- DCR <= 1, within the maximum -- or the design says it found none.

    Under ACI 318-19 and CIRSOC 201-25 the capacity it is judged by is also checked
    against strain compatibility written out apart from mento, with every bar and
    phi from eps_t: never above it, and a section that is not tension-controlled
    is reported as such. One pass over the grid does both, the designs being the
    costly part.
    """
    independent = code != "EN 1992-2004"
    silent: List[Any] = []
    overstated: List[Any] = []
    for b_cm, h_cm, f_c, M in _grid():
        beam, node = _designed(code, b_cm, h_cm, f_c, M)
        check = beam.flexure_checks[0]
        face = check.bottom if M > 0 else check.top
        codes = {warning.code for warning in node.warnings}
        assert face.M_capacity is not None and face.M_capacity.magnitude > 0, (b_cm, h_cm, f_c, M)
        if not codes & _CANNOT and (face.DCR > 1 + 1e-9 or "As_above_max" in codes):
            silent.append((b_cm, h_cm, f_c, M, round(face.DCR, 3), sorted(codes)))
        if not independent:
            continue
        if M > 0:
            A_s, A_sp = beam._A_s_bot.to("mm**2").magnitude, beam._A_s_top.to("mm**2").magnitude
            d, dp = beam._d_bot.to("mm").magnitude, beam._c_mec_top.to("mm").magnitude
        else:
            A_s, A_sp = beam._A_s_top.to("mm**2").magnitude, beam._A_s_bot.to("mm**2").magnitude
            d, dp = beam._d_top.to("mm").magnitude, beam._c_mec_bot.to("mm").magnitude
        phi_M_n, eps_t = _independent_capacity(A_s, A_sp, d, dp, b_cm * 10, f_c)
        # Up to A_s_max mento keeps the closed form, which leaves the compression
        # steel out: never more than the section has, and equal to it past the limit.
        reported = face.M_capacity.to("kN*m").magnitude
        if reported > phi_M_n * (1 + 1e-6):
            overstated.append((b_cm, h_cm, f_c, M, round(reported, 2), round(phi_M_n, 2)))
        if eps_t < _F_Y / _E_S + 0.003 - 1e-9:
            assert "As_above_max" in codes, (b_cm, h_cm, f_c, M, eps_t)
    assert not silent, f"{code}: designs that fail their own check without saying so: {silent}"
    assert not overstated, f"{code}: capacity above strain compatibility: {overstated}"


def test_a_full_design_passes_its_own_check_with_the_stirrups_it_ends_with() -> None:
    """ACI 318-19 20x50, f'c 25 MPa, ADN 420, c_c 25 mm, Mu = 150 kN·m, Vu = 250 kN.

    The flexure design runs at the depth of the 8 mm starter stirrup, where
    2Ø25 = 981.7 mm² carry the moment: d = 500 − 25 − 8 − 12.5 = 454.5 mm,
    a = 981.7·420 / (0.85·25·200) = 97.0 mm, φMn = 0.9·981.7·420·(454.5 − 48.5)
    = 150.66 kN·m, DCR 0.996. The shear design then settles on 1eØ10/11 for
    250 kN, the bars sink 2 mm to d = 452.5 mm and φMn = 0.9·412314·404.0 =
    149.92 kN·m: DCR 1.0005, and nothing warned, since the search never saw
    a shortfall (A_s,req = 9.823 cm² at that depth against 9.817 placed).
    ``design()`` now designs the flexure again with the stirrup it ended
    with: 2Ø25 + 1Ø20 at DCR 0.787, the same bars every time it runs. Fails
    before that with DCR 1.0005 and no warning.
    """
    beam = RectangularBeam(
        label="V",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=20 * cm,
        height=50 * cm,
        c_c=25 * mm,
    )
    node = Node(section=beam, forces=[Forces(label="U", M_y=150 * kNm, V_z=250 * kN)])
    node.design()
    node.check()
    first = str(beam.reinforcement.bottom)

    assert beam._stirrup_d_b.to("mm").magnitude == pytest.approx(10.0)
    assert beam.flexure_checks[0].bottom.DCR <= 1.0
    assert beam.shear_checks[0].DCR <= 1.0
    assert node.warnings == ()
    assert beam.flexure_design.bottom.A_s >= beam.flexure_design.bottom.A_s_req

    node.design()
    assert str(beam.reinforcement.bottom) == first


def test_a_lighter_stirrup_does_not_lift_the_minimum_past_the_bars_a_design_placed() -> None:
    """EN 1992-1-1 C25/30, B500S, 30x80, c_c 25 mm, M_Ed = 30 kN·m, V_Ed = 80 kN.

    A_s,min = 0.26·f_ctm/f_yk·b_t·d ≥ 0.0013·b_t·d (§9.2.1.1(1)), with
    f_ctm = 0.30·25^(2/3) = 2.565 MPa: 0.001334·b_t·d, which grows with d.
    At the Ø8 starter depth the design places 2Ø12 + 1Ø10 = 304.7 mm² on
    d = 761.3 mm against A_s,min = 0.001334·300·761.3 = 304.6 mm², and
    passes. The shear design then settles on 1eØ6/23, the bars rise 2 mm to
    d = 763.3 mm, A_s,min = 305.4 mm², and the design warned ``As_below_min``
    on its own bars. It now redoes the flexure with the Ø6 on the section and
    places 4Ø10 = 314.2 mm² (A_s,min 305.7 at the depth they sit). Fails
    before that with ``As_below_min``.
    """
    beam = RectangularBeam(
        label="V",
        concrete=Concrete_EN_1992_2004(name="C25", f_c=25 * MPa),
        steel_bar=SteelBar(name="B500S", f_y=500 * MPa),
        width=30 * cm,
        height=80 * cm,
        c_c=25 * mm,
    )
    node = Node(section=beam, forces=[Forces(label="U", M_y=30 * kNm, V_z=80 * kN)])
    node.design()
    node.check()

    assert beam._stirrup_d_b.to("mm").magnitude == pytest.approx(6.0)
    assert node.warnings == ()
    bottom = beam.flexure_design.bottom
    assert bottom.A_s >= bottom.A_s_min
    assert bottom.A_s.to("mm**2").magnitude == pytest.approx(314.16, abs=0.05)


@pytest.mark.parametrize("code", ["ACI 318-19", "CIRSOC 201-25"])
def test_the_steel_asked_for_grows_with_the_moment(code: str) -> None:
    """A_s,req on the tension face, and the compression it asks for, never drop as M grows.

    On a fixed section, from singly reinforced to three times the singly
    reinforced limit: the tension requirement used to fall back to A_s_max once
    the moment passed what tension steel alone carries, with no compression
    steel at all.
    """
    beam = RectangularBeam(
        label="M",
        concrete=_CONCRETES[code](name="C", f_c=25 * MPa),
        steel_bar=SteelBar(name="S", f_y=_F_Y * MPa),
        width=15 * cm,
        height=30 * cm,
        c_c=25 * mm,
    )
    beam.set_longitudinal_rebar_bot(2, 20 * mm, 0, None, 2, 20 * mm)
    beam.set_longitudinal_rebar_top(2, 20 * mm, 0, None, 2, 20 * mm)
    M_lim = _moment_at_singly_limit(15, 30, 25)
    moments = [M_lim * k / 10 for k in range(2, 31)]
    for sign in (1, -1):
        checks = beam.flexure_check_results([Forces(label=str(M), M_y=sign * M * kNm) for M in moments])
        tension = [(c.bottom if sign > 0 else c.top).A_s_req.to("cm**2").magnitude for c in checks]
        compression = [(c.top if sign > 0 else c.bottom).A_s_req.to("cm**2").magnitude for c in checks]
        for series in (tension, compression):
            assert all(b >= a - 1e-9 for a, b in zip(series, series[1:])), series
        # Past the singly reinforced limit the section asks for compression steel.
        assert compression[-1] > 0
        assert math.isfinite(tension[-1])
