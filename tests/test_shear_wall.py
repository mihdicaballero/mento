"""
Tests for ShearWall — ACI 318-19 Section 11 shear check and design.

Reference calculation (metric):
    lw=4.0 m, t=0.25 m, hw=3.5 m, f'c=25 MPa, fy=420 MPa
    hw/lw = 0.875 ≤ 1.5  →  α_c = 0.25
    Acv = 1.0 m²

    Ø12@200 mm:  ρt = 113.097/(250×200) = 0.002262
    Vc  = 0.25 × 5 × 1.0 MPa·m² = 1 250 kN
    Vs  = 0.002262 × 420 × 1.0 MPa·m² ≈ 950 kN
    Vn  = 2 200 kN,  Vn,max = 3 300 kN
    φVn = 1 650 kN,  φVn,max = 2 475 kN
    Vu=1 200 kN → DCR = 1 200/1 650 ≈ 0.727

    Ø12@150 mm:  ρt = 113.097/(250×150) = 0.003016  (≥ 0.0025 → passes)
    s_h,max = min(4000/5, 3×250, 450) = 450 mm
    s_v,max = min(4000/3, 3×250, 450) = 450 mm
"""

import math

import pytest

from mento.forces import Forces
from mento.material import Concrete_ACI_318_19, SteelBar
from mento.shear_wall import ShearWall
from mento.units import MPa, cm, kN, m, mm


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def wall_metric() -> ShearWall:
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)
    return ShearWall(
        label="W1",
        concrete=concrete,
        steel_bar=steel,
        thickness=25 * cm,
        length=4.0 * m,
        height=3.5 * m,
        c_c=20 * mm,
    )


@pytest.fixture
def wall_high_hw() -> ShearWall:
    """hw/lw = 9/4 = 2.25 ≥ 2.0  →  α_c = 0.17."""
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)
    return ShearWall(
        label="W_high",
        concrete=concrete,
        steel_bar=steel,
        thickness=25 * cm,
        length=4.0 * m,
        height=9.0 * m,
        c_c=20 * mm,
    )


@pytest.fixture
def wall_interp() -> ShearWall:
    """hw/lw = 7/4 = 1.75 (between 1.5 and 2.0)  →  α_c = 0.21."""
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN420", f_y=420 * MPa)
    return ShearWall(
        label="W_interp",
        concrete=concrete,
        steel_bar=steel,
        thickness=25 * cm,
        length=4.0 * m,
        height=7.0 * m,
        c_c=20 * mm,
    )


# ---------------------------------------------------------------------------
# Initialisation
# ---------------------------------------------------------------------------


class TestShearWallInit:
    def test_mode_is_shear_wall(self, wall_metric: ShearWall) -> None:
        assert wall_metric.mode == "shear_wall"

    def test_rho_t_starts_zero(self, wall_metric: ShearWall) -> None:
        assert wall_metric._rho_t.to("").magnitude == pytest.approx(0.0)

    def test_set_horizontal_rebar_updates_rho_t(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        expected = (math.pi / 4 * 12**2) / (250 * 200)  # mm²/(mm·mm)
        assert wall_metric._rho_t.to("").magnitude == pytest.approx(expected, rel=1e-4)

    def test_set_vertical_rebar_updates_rho_l(self, wall_metric: ShearWall) -> None:
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=150 * mm)
        expected = (math.pi / 4 * 12**2) / (250 * 150)
        assert wall_metric._rho_l.to("").magnitude == pytest.approx(expected, rel=1e-4)


# ---------------------------------------------------------------------------
# α_c
# ---------------------------------------------------------------------------


class TestAlphaC:
    def test_alpha_c_low_hw_lw(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=500 * kN)])
        assert wall_metric._hw_lw == pytest.approx(0.875, rel=1e-4)
        assert wall_metric._alpha_c == pytest.approx(0.25, rel=1e-4)

    def test_alpha_c_high_hw_lw(self, wall_high_hw: ShearWall) -> None:
        wall_high_hw.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_high_hw.check_shear([Forces(V_z=500 * kN)])
        assert wall_high_hw._alpha_c == pytest.approx(0.17, rel=1e-4)

    def test_alpha_c_interpolated(self, wall_interp: ShearWall) -> None:
        wall_interp.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_interp.check_shear([Forces(V_z=500 * kN)])
        # hw/lw=1.75 → t=0.5  → α_c = 0.25 + 0.5*(0.17-0.25) = 0.21
        assert wall_interp._alpha_c == pytest.approx(0.21, rel=1e-4)


# ---------------------------------------------------------------------------
# Shear strength
# ---------------------------------------------------------------------------


class TestShearStrength:
    def test_Vc(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=1200 * kN)])
        # Vc = 0.25 × 1.0 × √25 MPa × 1.0 m² = 1 250 kN
        assert wall_metric._V_c_wall.to("kN").magnitude == pytest.approx(1250.0, rel=1e-3)

    def test_phi_Vc(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=1200 * kN)])
        assert (wall_metric.concrete.phi_v * wall_metric._V_c_wall).to("kN").magnitude == pytest.approx(937.5, rel=1e-3)

    def test_phi_Vn(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=1200 * kN)])
        assert wall_metric._phi_V_n_wall.to("kN").magnitude == pytest.approx(1650.0, rel=1e-2)

    def test_phi_Vn_max(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=1200 * kN)])
        assert wall_metric._phi_V_n_max_wall.to("kN").magnitude == pytest.approx(2475.0, rel=1e-2)

    def test_DCR(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        result = wall_metric.check_shear([Forces(V_z=1200 * kN)])
        # data row is index 1 (index 0 is units row)
        dcr = float(result.iloc[1]["DCR"])
        assert dcr == pytest.approx(1200.0 / 1650.0, rel=1e-3)

    def test_Vu_le_phi_Vn(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        result = wall_metric.check_shear([Forces(V_z=1200 * kN)])
        assert result.iloc[1]["Vu≤ØVn"] is True

    def test_Vu_le_phi_Vn_max(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        result = wall_metric.check_shear([Forces(V_z=1200 * kN)])
        assert result.iloc[1]["Vu≤ØVn,max"] is True


# ---------------------------------------------------------------------------
# Minimum reinforcement
# ---------------------------------------------------------------------------


class TestMinReinforcement:
    def test_rho_t_min_is_0025(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=500 * kN)])
        assert wall_metric._rho_t_min.to("").magnitude == pytest.approx(0.0025, rel=1e-4)

    def test_rho_l_min_is_0025_when_hw_lw_gt2(self, wall_high_hw: ShearWall) -> None:
        # hw/lw=2.25 > 2.0 → ρl,min = 0.0025 regardless of ρt
        wall_high_hw.set_horizontal_rebar(d_b=16 * mm, s=100 * mm)
        wall_high_hw.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall_high_hw.check_shear([Forces(V_z=500 * kN)])
        assert wall_high_hw._rho_l_min.to("").magnitude == pytest.approx(0.0025, rel=1e-4)

    def test_rho_l_min_equals_rho_t_when_hw_lw_le2(self, wall_metric: ShearWall) -> None:
        # hw/lw=0.875 ≤ 2.0 → ρl,min = max(0.0025, ρt)
        # Ø16@100: ρt = π/4×256/(250×100) = 201.06/25000 ≈ 0.00804
        wall_metric.set_horizontal_rebar(d_b=16 * mm, s=100 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=500 * kN)])
        expected_rho_t = (math.pi / 4 * 16**2) / (250 * 100)
        assert wall_metric._rho_l_min.to("").magnitude == pytest.approx(expected_rho_t, rel=1e-4)

    def test_rho_t_below_min_flagged(self, wall_metric: ShearWall) -> None:
        # Ø12@200 → ρt=0.002262 < 0.0025 → ❌
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=500 * kN)])
        assert wall_metric._data_min_max_wall["Ok?"][0] == "❌"

    def test_rho_t_above_min_passes(self, wall_metric: ShearWall) -> None:
        # Ø12@150 → ρt = 113.097/(250×150) ≈ 0.003016 > 0.0025 → ✅
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=150 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=150 * mm)
        wall_metric.check_shear([Forces(V_z=500 * kN)])
        assert wall_metric._data_min_max_wall["Ok?"][0] == "✅"


# ---------------------------------------------------------------------------
# Spacing limits
# ---------------------------------------------------------------------------


class TestSpacingLimits:
    def test_s_h_max_metric(self, wall_metric: ShearWall) -> None:
        # min(4000/5, 3×250, 450) = min(800, 750, 450) = 450 mm
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=500 * kN)])
        assert wall_metric._s_h_max.to("mm").magnitude == pytest.approx(450.0, rel=1e-4)

    def test_s_v_max_metric(self, wall_metric: ShearWall) -> None:
        # min(4000/3, 3×250, 450) = min(1333, 750, 450) = 450 mm
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=500 * kN)])
        assert wall_metric._s_v_max.to("mm").magnitude == pytest.approx(450.0, rel=1e-4)

    def test_spacing_ok_flag_when_within_limit(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=500 * kN)])
        # s_h=200 ≤ 450 → ✅
        assert wall_metric._data_min_max_wall["Ok?"][2] == "✅"
        assert wall_metric._data_min_max_wall["Ok?"][3] == "✅"

    def test_spacing_fail_flag_when_exceeds_limit(self) -> None:
        concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
        steel = SteelBar(name="ADN420", f_y=420 * MPa)
        # Small wall: lw=500 mm, t=150 mm → s_h,max = min(100, 450, 450) = 100 mm
        wall = ShearWall(
            label="W_small",
            concrete=concrete,
            steel_bar=steel,
            thickness=150 * mm,
            length=500 * mm,
            height=3000 * mm,
            c_c=20 * mm,
        )
        wall.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)  # 200 > s_h,max
        wall.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall.check_shear([Forces(V_z=100 * kN)])
        assert wall._data_min_max_wall["Ok?"][2] == "❌"


# ---------------------------------------------------------------------------
# DataFrame output
# ---------------------------------------------------------------------------


class TestDataFrameOutput:
    def test_result_columns(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        result = wall_metric.check_shear([Forces(V_z=1200 * kN)])
        expected_cols = [
            "Label",
            "Comb.",
            "ρt,min",
            "ρt,req",
            "ρt",
            "ρl,min",
            "ρl",
            "Vu",
            "ØVc",
            "ØVs",
            "ØVn",
            "ØVn,max",
            "Vu≤ØVn,max",
            "Vu≤ØVn",
            "DCR",
        ]
        assert list(result.columns) == expected_cols

    def test_result_two_rows_for_one_force(self, wall_metric: ShearWall) -> None:
        # Row 0 = units header, row 1 = data
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        result = wall_metric.check_shear([Forces(V_z=1200 * kN)])
        assert len(result) == 2

    def test_units_row_has_kN(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        result = wall_metric.check_shear([Forces(V_z=500 * kN)])
        assert result.iloc[0]["Vu"] == "kN"

    def test_shear_checked_flag_set(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=500 * kN)])
        assert wall_metric._shear_checked is True
        assert wall_metric._shear_wall_checked is True


# ---------------------------------------------------------------------------
# Multiple forces
# ---------------------------------------------------------------------------


class TestMultipleForces:
    def test_three_rows_for_two_forces(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=150 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=150 * mm)
        f1 = Forces(label="C1", V_z=500 * kN)
        f2 = Forces(label="C2", V_z=1500 * kN)
        result = wall_metric.check_shear([f1, f2])
        assert len(result) == 3  # units row + 2 data rows

    def test_limiting_case_has_highest_dcr(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=150 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=150 * mm)
        f1 = Forces(label="C1", V_z=500 * kN)
        f2 = Forces(label="C2", V_z=1500 * kN)
        result = wall_metric.check_shear([f1, f2])
        data = result.iloc[1:]  # skip units row
        data_dcr = data["DCR"].astype(float)
        assert data_dcr.idxmax() == data.index[-1]  # C2 is the worst case


# ---------------------------------------------------------------------------
# Design
# ---------------------------------------------------------------------------


class TestDesignShear:
    def test_design_returns_dataframe(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        result = wall_metric.design_shear([Forces(V_z=500 * kN)])
        import pandas as pd

        assert isinstance(result, pd.DataFrame)

    def test_rho_t_req_at_least_rho_t_min(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.design_shear([Forces(V_z=200 * kN)])
        assert wall_metric._rho_t_req.to("").magnitude >= 0.0025 - 1e-9
