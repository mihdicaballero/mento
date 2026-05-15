"""
Tests for ShearWall — ACI 318-19 Section 11 shear check and design.

Reference calculation (metric):
    lw=4.0 m, t=0.25 m, hw=3.5 m, f'c=25 MPa, fy=420 MPa
    hw/lw = 0.875 ≤ 1.5  →  α_c = 0.25
    Acv = 1.0 m²

    Distributed mesh is on both faces (E.F.): ρ = 2·Ab/(t·s).

    Ø12@200 mm E.F.:  ρt = 2×113.097/(250×200) = 0.004524
    Vc  = 0.25 × 5 × 1.0 MPa·m² = 1 250 kN
    Vs  = 0.004524 × 420 × 1.0 MPa·m² ≈ 1 900 kN
    Vn  = 3 150 kN,  Vn,max = 3 300 kN
    φVn = 2 362.5 kN,  φVn,max = 2 475 kN
    Vu=1 200 kN → DCR = 1 200/2 362.5 ≈ 0.508

    Ø12@150 mm E.F.:  ρt = 2×113.097/(250×150) = 0.006032  (≥ 0.0025 → passes)
    s_h,max = min(4000/5, 3×250, 450) = 450 mm
    s_v,max = min(4000/3, 3×250, 450) = 450 mm
"""

import math

import matplotlib
import pytest

matplotlib.use("Agg")  # headless backend — no display needed for plot tests

from mento.forces import Forces
from mento.material import Concrete_ACI_318_19, Concrete_EN_1992_2004, SteelBar
from mento.shear_wall import ShearWall
from mento.units import MPa, cm, inch, kip, kN, ksi, m, mm, psi


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
        # Mesh on both faces (E.F.): ρt = 2·Ab/(t·s)
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        expected = 2 * (math.pi / 4 * 12**2) / (250 * 200)  # mm²/(mm·mm)
        assert wall_metric._rho_t.to("").magnitude == pytest.approx(expected, rel=1e-4)

    def test_set_vertical_rebar_updates_rho_l(self, wall_metric: ShearWall) -> None:
        # Mesh on both faces (E.F.): ρl = 2·Ab/(t·s)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=150 * mm)
        expected = 2 * (math.pi / 4 * 12**2) / (250 * 150)
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
        # Ø12@200 E.F.: Vs ≈ 1900 kN → Vn ≈ 3150 kN → φVn ≈ 2362.5 kN
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=1200 * kN)])
        assert wall_metric._phi_V_n_wall.to("kN").magnitude == pytest.approx(2362.5, rel=1e-2)

    def test_phi_Vn_max(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=1200 * kN)])
        assert wall_metric._phi_V_n_max_wall.to("kN").magnitude == pytest.approx(2475.0, rel=1e-2)

    def test_DCR(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        result = wall_metric.check_shear([Forces(V_z=1200 * kN)])
        # data row is index 1 (index 0 is units row)
        dcr = float(result.iloc[1]["DCR"])
        assert dcr == pytest.approx(1200.0 / 2362.5, rel=1e-2)

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

    def test_rho_l_min_interpolated_per_11_6_2(self, wall_metric: ShearWall) -> None:
        # ACI 318-19 §11.6.2 Eq.(11.6.2):
        #   ρl,min = max(0.0025, 0.0025 + 0.5·(2.5 − r_hw)·(ρt,req − 0.0025))
        # High Vu so ρt,req exceeds the 0.0025 minimum and the interpolation bites.
        wall_metric.set_horizontal_rebar(d_b=16 * mm, s=100 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=2000 * kN)])
        rho_t_req = wall_metric._rho_t_req.to("").magnitude
        assert rho_t_req > 0.0025  # interpolation is non-trivial
        r_hw = max(0.5, min(wall_metric._hw_lw, 2.5))
        expected = max(0.0025, 0.0025 + 0.5 * (2.5 - r_hw) * (rho_t_req - 0.0025))
        assert wall_metric._rho_l_min.to("").magnitude == pytest.approx(expected, rel=1e-4)

    def test_rho_t_below_min_flagged(self, wall_metric: ShearWall) -> None:
        # Ø12@400 E.F. → ρt = 2×113.097/(250×400) = 0.002262 < 0.0025 → ❌
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=400 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=400 * mm)
        wall_metric.check_shear([Forces(V_z=500 * kN)])
        assert wall_metric._data_min_max_wall["Ok?"][0] == "❌"

    def test_rho_t_above_min_passes(self, wall_metric: ShearWall) -> None:
        # Ø12@150 E.F. → ρt = 2×113.097/(250×150) ≈ 0.006032 > 0.0025 → ✅
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
        wall_metric.design_shear([Forces(V_z=200 * kN)])
        assert wall_metric._rho_t_req.to("").magnitude >= 0.0025 - 1e-9

    def test_design_assigns_mesh(self, wall_metric: ShearWall) -> None:
        """design_shear must select and apply a valid mesh in both directions."""
        wall_metric.design_shear([Forces(V_z=1200 * kN)])
        assert wall_metric._s_h.magnitude > 0
        assert wall_metric._s_v.magnitude > 0
        assert wall_metric._rho_t.to("").magnitude >= wall_metric._rho_t_req.to("").magnitude - 1e-9
        assert wall_metric._rho_l.to("").magnitude >= wall_metric._rho_l_min.to("").magnitude - 1e-9
        assert wall_metric._s_h <= wall_metric._s_h_max
        assert wall_metric._s_v <= wall_metric._s_v_max

    def test_design_respects_bar_cap(self, wall_metric: ShearWall) -> None:
        """A low-demand wall should stay at or below the Ø12 mm crack-control cap."""
        wall_metric.design_shear([Forces(V_z=600 * kN)])
        assert wall_metric._d_b_h <= 12 * mm
        assert wall_metric._d_b_v <= 12 * mm

    def test_design_spacing_on_grid(self, wall_metric: ShearWall) -> None:
        """Selected spacings are integer centimetres on the 2.5 cm-derived grid."""
        wall_metric.design_shear([Forces(V_z=900 * kN)])
        for s in (wall_metric._s_h, wall_metric._s_v):
            s_cm = s.to("cm").magnitude
            assert abs(s_cm - round(s_cm)) < 1e-6  # whole centimetres
            assert s_cm >= 5.0  # practical floor

    def test_design_worst_case_across_forces(self, wall_metric: ShearWall) -> None:
        """The designed mesh must satisfy the most demanding combination."""
        result = wall_metric.design_shear([Forces(label="c1", V_z=400 * kN), Forces(label="c2", V_z=1300 * kN)])
        assert result["DCR"].iloc[1:].astype(float).max() <= 1.0

    def test_design_cirsoc_allows_6mm_transverse(self) -> None:
        """CIRSOC vertical mesh stays ≥ Ø10 mm; transverse may use Ø6 mm."""
        from mento.material import Concrete_CIRSOC_201_25

        conc = Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa)
        steel = SteelBar(name="ADN420", f_y=420 * MPa)
        wall = ShearWall(
            label="WC",
            concrete=conc,
            steel_bar=steel,
            thickness=20 * cm,
            length=3.0 * m,
            height=2.5 * m,
            c_c=20 * mm,
        )
        wall.design_shear([Forces(V_z=150 * kN)])
        assert wall._d_b_v >= 10 * mm  # CIRSOC vertical minimum
        assert wall._d_b_h >= 6 * mm  # CIRSOC transverse minimum

    def test_design_empty_forces_raises(self, wall_metric: ShearWall) -> None:
        with pytest.raises(ValueError):
            wall_metric.design_shear([])

    def test_design_cirsoc_reference_wall(self) -> None:
        """CIRSOC 201-25 reference wall — full design + check.

        Shear Wall W1: lw=400 cm, t=20 cm, hw=350 cm, c_c=2 cm,
        Concrete C25, Rebar ADN 420.
        Expected design (Vu = 800 kN):
            Horizontal rebar:        Ø8/20 cm E.F.  → ρt = 0.00251
            Minimum vertical rebar:  Ø10/30 cm E.F. → ρl = 0.00262
            φVn = 1383.35 kN  →  DCR = 0.578
        """
        from mento.material import Concrete_CIRSOC_201_25

        concrete = Concrete_CIRSOC_201_25(name="C25", f_c=25 * MPa)
        steel = SteelBar(name="ADN 420", f_y=420 * MPa)
        wall = ShearWall(
            label="W1",
            concrete=concrete,
            steel_bar=steel,
            thickness=20 * cm,
            length=400 * cm,
            height=350 * cm,
            c_c=2 * cm,
        )
        wall.design_shear([Forces(V_z=800 * kN)])

        # Horizontal (transverse) mesh
        assert wall._d_b_h.to("mm").magnitude == pytest.approx(8.0)
        assert wall._s_h.to("cm").magnitude == pytest.approx(20.0)
        assert wall._rho_t.to("").magnitude == pytest.approx(0.00251, abs=1e-5)

        # Minimum vertical mesh
        assert wall._d_b_v.to("mm").magnitude == pytest.approx(10.0)
        assert wall._s_v.to("cm").magnitude == pytest.approx(30.0)
        assert wall._rho_l.to("").magnitude == pytest.approx(0.00262, abs=1e-5)

        # Capacity and demand-capacity ratio
        assert wall._phi_V_n_wall.to("kN").magnitude == pytest.approx(1383.35, rel=1e-3)
        assert wall._DCRv_wall == pytest.approx(0.578, abs=1e-3)


# ---------------------------------------------------------------------------
# Imperial unit system
# ---------------------------------------------------------------------------


@pytest.fixture
def wall_imperial() -> ShearWall:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steel = SteelBar(name="G60", f_y=60 * ksi)
    return ShearWall(
        label="WI",
        concrete=concrete,
        steel_bar=steel,
        thickness=10 * inch,
        length=160 * inch,
        height=120 * inch,
        c_c=1 * inch,
    )


class TestImperialWall:
    def test_imperial_check_units_row(self, wall_imperial: ShearWall) -> None:
        wall_imperial.set_horizontal_rebar(d_b=0.5 * inch, s=12 * inch)
        result = wall_imperial.check_shear([Forces(V_z=100 * kip)])
        assert result.iloc[0]["Vu"] == "kip"

    def test_imperial_alpha_c(self, wall_imperial: ShearWall) -> None:
        # hw/lw = 120/160 = 0.75 ≤ 1.5 → α_c = 3.0 (imperial)
        wall_imperial.set_horizontal_rebar(d_b=0.5 * inch, s=12 * inch)
        wall_imperial.check_shear([Forces(V_z=100 * kip)])
        assert wall_imperial._alpha_c == pytest.approx(3.0, rel=1e-4)

    def test_imperial_design(self, wall_imperial: ShearWall) -> None:
        wall_imperial.design_shear([Forces(V_z=120 * kip)])
        assert wall_imperial._s_h.magnitude > 0
        assert wall_imperial._s_v.magnitude > 0
        assert wall_imperial._rho_t.to("").magnitude >= wall_imperial._rho_t_req.to("").magnitude - 1e-9


# ---------------------------------------------------------------------------
# Reporting, display, plot, and error paths (coverage)
# ---------------------------------------------------------------------------


class TestShearWallReporting:
    # --- unsupported design code ----------------------------------------
    def _en_wall(self) -> ShearWall:
        concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
        steel = SteelBar(name="ADN420", f_y=420 * MPa)
        return ShearWall(
            label="WE",
            concrete=concrete,
            steel_bar=steel,
            thickness=20 * cm,
            length=4.0 * m,
            height=3.5 * m,
            c_c=20 * mm,
        )

    def test_check_shear_unsupported_code_raises(self) -> None:
        with pytest.raises(NotImplementedError):
            self._en_wall().check_shear([Forces(V_z=100 * kN)])

    def test_design_shear_unsupported_code_raises(self) -> None:
        with pytest.raises(NotImplementedError):
            self._en_wall().design_shear([Forces(V_z=100 * kN)])

    # --- check()/design() wrappers --------------------------------------
    def test_check_wrapper(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        assert wall_metric.check([Forces(V_z=500 * kN)]) is None
        assert wall_metric._shear_wall_checked is True

    def test_design_wrapper(self, wall_metric: ShearWall) -> None:
        assert wall_metric.design([Forces(V_z=500 * kN)]) is None
        assert wall_metric._s_h.magnitude > 0

    # --- flexure stubs ---------------------------------------------------
    def test_flexure_methods_raise(self, wall_metric: ShearWall) -> None:
        f = Forces(V_z=1 * kN)
        with pytest.raises(NotImplementedError):
            wall_metric.check_flexure([f])
        with pytest.raises(NotImplementedError):
            wall_metric.design_flexure([f])
        with pytest.raises(NotImplementedError):
            wall_metric.flexure_results_detailed()
        with pytest.raises(NotImplementedError):
            wall_metric.flexure_results_detailed_doc()

    # --- data / shear_results / results ---------------------------------
    def test_data_property(self, wall_metric: ShearWall) -> None:
        assert wall_metric.data is None
        assert "Shear Wall" in wall_metric._md_data

    def test_shear_results_not_checked(self, wall_metric: ShearWall) -> None:
        assert wall_metric.shear_results is None
        assert wall_metric._md_shear_results == "Shear results are not available."

    def test_shear_results_after_check(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=800 * kN)])
        assert wall_metric.shear_results is None
        assert "Horizontal rebar" in wall_metric._md_shear_results
        assert "E.F." in wall_metric._md_shear_results

    def test_shear_results_no_rebar_assigned(self, wall_metric: ShearWall) -> None:
        wall_metric.check_shear([Forces(V_z=200 * kN)])
        wall_metric.shear_results
        assert "not assigned" in wall_metric._md_shear_results

    def test_shear_results_no_capacity_branch(self, wall_metric: ShearWall) -> None:
        # White-box: checked flag set but no limiting-case details available.
        wall_metric._shear_wall_checked = True
        wall_metric._limiting_case_shear_details = None
        assert wall_metric.shear_results is None
        assert wall_metric._md_shear_results == "No shear to check."

    def test_results_property(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=800 * kN)])
        assert wall_metric.results is None

    # --- detailed results ------------------------------------------------
    def test_shear_results_detailed_not_checked(self, wall_metric: ShearWall) -> None:
        assert wall_metric.shear_results_detailed() is None

    def test_shear_results_detailed_limiting(self, wall_metric: ShearWall, capsys) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=800 * kN)])
        wall_metric.shear_results_detailed()
        assert "SHEAR WALL DETAILED RESULTS" in capsys.readouterr().out

    def test_shear_results_detailed_specific_force(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        f = Forces(V_z=800 * kN)
        wall_metric.check_shear([f])
        wall_metric.shear_results_detailed(f)  # explicit force

    def test_shear_results_detailed_bad_force_raises(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=800 * kN)])
        with pytest.raises(ValueError):
            wall_metric.shear_results_detailed(Forces(V_z=1 * kN))

    def test_shear_results_detailed_doc(self, wall_metric: ShearWall, tmp_path, monkeypatch) -> None:
        monkeypatch.chdir(tmp_path)
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=800 * kN)])
        wall_metric.shear_results_detailed_doc()
        assert any(p.suffix == ".docx" for p in tmp_path.iterdir())

    def test_shear_results_detailed_doc_not_checked(self, wall_metric: ShearWall) -> None:
        assert wall_metric.shear_results_detailed_doc() is None

    def test_shear_results_detailed_doc_specific_force(self, wall_metric: ShearWall, tmp_path, monkeypatch) -> None:
        monkeypatch.chdir(tmp_path)
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=12 * mm, s=200 * mm)
        f = Forces(V_z=800 * kN)
        wall_metric.check_shear([f])
        wall_metric.shear_results_detailed_doc(f)
        assert any(p.suffix == ".docx" for p in tmp_path.iterdir())

    def test_shear_results_detailed_doc_bad_force_raises(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.check_shear([Forces(V_z=800 * kN)])
        with pytest.raises(ValueError):
            wall_metric.shear_results_detailed_doc(Forces(V_z=1 * kN))

    # --- plot ------------------------------------------------------------
    def test_plot_with_rebar(self, wall_metric: ShearWall) -> None:
        wall_metric.set_horizontal_rebar(d_b=12 * mm, s=200 * mm)
        wall_metric.set_vertical_rebar(d_b=10 * mm, s=300 * mm)
        fig = wall_metric.plot()
        assert fig is not None

    def test_plot_without_rebar(self, wall_metric: ShearWall) -> None:
        import warnings

        # show=True also exercises the plt.show() branch; under the Agg
        # backend it emits a harmless "non-interactive" warning we silence.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fig = wall_metric.plot(show=True)
        assert fig is not None


# ---------------------------------------------------------------------------
# Code-function guards (white-box)
# ---------------------------------------------------------------------------


class TestWallCodeGuards:
    def test_check_function_rejects_non_aci_concrete(self) -> None:
        from mento.codes.ACI_318_19_wall import _check_shear_ACI_318_19_wall

        concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
        steel = SteelBar(name="ADN420", f_y=420 * MPa)
        wall = ShearWall(
            label="WE",
            concrete=concrete,
            steel_bar=steel,
            thickness=20 * cm,
            length=4.0 * m,
            height=3.5 * m,
            c_c=20 * mm,
        )
        with pytest.raises(TypeError):
            _check_shear_ACI_318_19_wall(wall, Forces(V_z=100 * kN))

    def test_design_core_empty_forces_raises(self, wall_metric: ShearWall) -> None:
        from mento.codes.ACI_318_19_wall import _ACI_WALL_BARS_METRIC, _design_shear_wall_core

        with pytest.raises(ValueError):
            _design_shear_wall_core(wall_metric, [], _ACI_WALL_BARS_METRIC, _ACI_WALL_BARS_METRIC)

    def test_design_no_bar_fits_raises(self) -> None:
        # Very short wall → s_h,max < 5 cm practical floor → empty grid → ValueError.
        concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
        steel = SteelBar(name="ADN420", f_y=420 * MPa)
        wall = ShearWall(
            label="WT",
            concrete=concrete,
            steel_bar=steel,
            thickness=10 * cm,
            length=20 * cm,
            height=200 * cm,
            c_c=20 * mm,
        )
        with pytest.raises(ValueError):
            wall.design_shear([Forces(V_z=50 * kN)])
