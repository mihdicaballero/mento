"""Tests for ShearWallSummary class."""

import pytest
import pandas as pd
import os

from mento import (
    Concrete_ACI_318_19,
    SteelBar,
    ShearWallSummary,
    MPa,
    psi,
    ksi,
)


# ------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------


@pytest.fixture
def concrete():
    return Concrete_ACI_318_19(name="H25", f_c=25 * MPa)


@pytest.fixture
def steel():
    return SteelBar(name="ADN 420", f_y=420 * MPa)


@pytest.fixture
def sample_df():
    """DataFrame with 4 walls: (Level 1, M1), (Level 2, M1), (Level 1, M2), (Level 2, M2)."""
    data = {
        "Level": [
            "",
            "Level 1",
            "Level 1",
            "Level 1",
            "Level 1",
            "Level 2",
            "Level 2",
            "Level 2",
            "Level 2",
            "Level 1",
            "Level 1",
            "Level 1",
            "Level 1",
            "Level 2",
            "Level 2",
            "Level 2",
            "Level 2",
        ],
        "Label": ["", "M1", "M1", "M1", "M1", "M1", "M1", "M1", "M1", "M2", "M2", "M2", "M2", "M2", "M2", "M2", "M2"],
        "Comb.": [
            "",
            "ELU 1",
            "ELU 2",
            "ELU 3",
            "ELU 4",
            "ELU 1",
            "ELU 2",
            "ELU 3",
            "ELU 4",
            "ELU 1",
            "ELU 2",
            "ELU 3",
            "ELU 4",
            "ELU 1",
            "ELU 2",
            "ELU 3",
            "ELU 4",
        ],
        "t": ["cm", 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20],
        "lw": ["m", 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
        "hw": ["m", 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0],
        "cc": ["mm", 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25],
        "Nx": ["kN", 0, 0, 0, -301, -150, 55.5, 282, -4.5, -240, -163, -17, 332, -150, 55.5, -163, 55.5],
        "Vz": ["kN", 264, 138, 123, 152, 32.3, 163, 19, 88.15, 61.2, 29, 47, 21, 32.3, 163, 29, 163],
        "My": ["kNm", -172, -90, -81, -234, 143, -278, 159, -97, -38, 60, -39, 46.13, 143, -278, 60, -278],
        "dbh": ["mm", 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
        "sh": ["cm", 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20],
        "dbv": ["mm", 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
        "sv": ["cm", 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15],
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_df_no_rebar():
    """Same as sample_df but with no rebar assigned."""
    data = {
        "Level": ["", "Level 1", "Level 1", "Level 2", "Level 2"],
        "Label": ["", "M1", "M1", "M1", "M1"],
        "Comb.": ["", "ELU 1", "ELU 2", "ELU 1", "ELU 2"],
        "t": ["cm", 20, 20, 20, 20],
        "lw": ["m", 3.0, 3.0, 3.0, 3.0],
        "hw": ["m", 3.0, 3.0, 3.0, 3.0],
        "cc": ["mm", 25, 25, 25, 25],
        "Nx": ["kN", 0, -301, -150, 55.5],
        "Vz": ["kN", 264, 152, 32.3, 163],
        "My": ["kNm", -172, -234, 143, -278],
        "dbh": ["mm", 0, 0, 0, 0],
        "sh": ["cm", 0, 0, 0, 0],
        "dbv": ["mm", 0, 0, 0, 0],
        "sv": ["cm", 0, 0, 0, 0],
    }
    return pd.DataFrame(data)


@pytest.fixture
def summary(concrete, steel, sample_df):
    return ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=sample_df)


# ------------------------------------------------------------------
# Init / Grouping
# ------------------------------------------------------------------


class TestShearWallSummaryInit:
    def test_correct_number_of_nodes(self, summary):
        assert len(summary.nodes) == 4

    def test_wall_keys(self, summary):
        assert summary.wall_keys == [
            ("Level 1", "M1"),
            ("Level 2", "M1"),
            ("Level 1", "M2"),
            ("Level 2", "M2"),
        ]

    def test_forces_per_node(self, summary):
        for node in summary.nodes:
            assert len(node.forces) == 4

    def test_wall_level_attribute(self, summary):
        wall_0 = summary.nodes[0].section
        assert wall_0.level == "Level 1"
        wall_1 = summary.nodes[1].section
        assert wall_1.level == "Level 2"

    def test_wall_label_attribute(self, summary):
        wall_0 = summary.nodes[0].section
        assert wall_0.label == "M1"
        wall_2 = summary.nodes[2].section
        assert wall_2.label == "M2"

    def test_wall_geometry(self, summary):
        wall = summary.nodes[0].section
        assert wall.thickness.to("cm").magnitude == pytest.approx(20)
        assert wall.length.to("m").magnitude == pytest.approx(3.0)
        assert wall.height.to("m").magnitude == pytest.approx(3.0)

    def test_rebar_set(self, summary):
        wall = summary.nodes[0].section
        assert wall._d_b_h.to("mm").magnitude == pytest.approx(8)
        assert wall._s_h.to("cm").magnitude == pytest.approx(20)
        assert wall._d_b_v.to("mm").magnitude == pytest.approx(12)
        assert wall._s_v.to("cm").magnitude == pytest.approx(15)

    def test_no_rebar_walls(self, concrete, steel, sample_df_no_rebar):
        s = ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=sample_df_no_rebar)
        wall = s.nodes[0].section
        assert wall._d_b_h.magnitude == 0  # type: ignore
        assert wall._d_b_v.magnitude == 0  # type: ignore


# ------------------------------------------------------------------
# Input validation
# ------------------------------------------------------------------


class TestShearWallSummaryValidation:
    def test_invalid_units(self, concrete: Concrete_ACI_318_19, steel: SteelBar):
        data = {
            "Level": [""],
            "Label": [""],
            "Comb.": [""],
            "t": ["parsecs"],
            "lw": ["m"],
            "hw": ["m"],
            "cc": ["mm"],
            "Nx": ["kN"],
            "Vz": ["kN"],
            "My": ["kNm"],
            "dbh": ["mm"],
            "sh": ["cm"],
            "dbv": ["mm"],
            "sv": ["cm"],
        }
        df = pd.DataFrame(data)
        with pytest.raises(ValueError, match="Invalid unit"):
            ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=df)

    def test_geometry_mismatch_raises(self, concrete: Concrete_ACI_318_19, steel: SteelBar):
        data = {
            "Level": ["", "Level 1", "Level 1"],
            "Label": ["", "M1", "M1"],
            "Comb.": ["", "ELU 1", "ELU 2"],
            "t": ["cm", 20, 25],  # mismatch!
            "lw": ["m", 3.0, 3.0],
            "hw": ["m", 3.0, 3.0],
            "cc": ["mm", 25, 25],
            "Nx": ["kN", 0, 0],
            "Vz": ["kN", 264, 138],
            "My": ["kNm", -172, -90],
            "dbh": ["mm", 8, 8],
            "sh": ["cm", 20, 20],
            "dbv": ["mm", 12, 12],
            "sv": ["cm", 15, 15],
        }
        df = pd.DataFrame(data)
        with pytest.raises(ValueError, match="Geometry mismatch"):
            ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=df)


# ------------------------------------------------------------------
# Check
# ------------------------------------------------------------------


class TestShearWallSummaryCheck:
    def test_check_returns_correct_rows(self, summary):
        df = summary.check()
        # 1 units row + 4 wall rows
        assert len(df) == 5

    def test_check_has_dcr_column(self, summary):
        df = summary.check()
        assert "DCR" in df.columns

    def test_check_has_status_column(self, summary):
        df = summary.check()
        assert "Status" in df.columns
        # All should pass for this input
        statuses = df["Status"].iloc[1:].tolist()
        assert all(s == "✅" for s in statuses)

    def test_check_dcr_values(self, summary):
        df = summary.check()
        dcr_values = df["DCR"].iloc[1:].tolist()
        for dcr in dcr_values:
            assert 0 < dcr < 1

    def test_check_raises_without_rebar(self, concrete, steel, sample_df_no_rebar):
        s = ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=sample_df_no_rebar)
        with pytest.raises(ValueError, match="no horizontal rebar"):
            s.check()


# ------------------------------------------------------------------
# Design
# ------------------------------------------------------------------


class TestShearWallSummaryDesign:
    def test_design_returns_dataframe(self, concrete, steel, sample_df_no_rebar):
        s = ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=sample_df_no_rebar)
        result = s.design()
        assert isinstance(result, pd.DataFrame)

    def test_design_fills_rebar(self, concrete, steel, sample_df_no_rebar):
        s = ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=sample_df_no_rebar)
        s.design()
        wall = s.nodes[0].section
        assert wall._d_b_h.to("mm").magnitude > 0  # type: ignore
        assert wall._s_h.to("cm").magnitude > 0  # type: ignore
        assert wall._d_b_v.to("mm").magnitude > 0  # type: ignore
        assert wall._s_v.to("cm").magnitude > 0  # type: ignore

    def test_design_sets_design_data(self, concrete, steel, sample_df_no_rebar):
        s = ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=sample_df_no_rebar)
        s.design()
        assert hasattr(s, "design_data")

    def test_check_after_design_passes(self, concrete, steel, sample_df_no_rebar):
        s = ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=sample_df_no_rebar)
        s.design()
        df = s.check()
        statuses = df["Status"].iloc[1:].tolist()
        assert all(s == "✅" for s in statuses)


# ------------------------------------------------------------------
# Shear results
# ------------------------------------------------------------------


class TestShearWallSummaryShearResults:
    def test_shear_results_all(self, summary):
        df = summary.shear_results()
        # 1 units row + 16 data rows (4 walls × 4 combos)
        assert len(df) == 17

    def test_shear_results_single(self, summary):
        df = summary.shear_results(index=1)
        # 1 units row + 4 combos for wall 1
        assert len(df) == 5

    def test_shear_results_index_out_of_range(self, summary):
        with pytest.raises(IndexError):
            summary.shear_results(index=99)

    def test_shear_results_index_zero(self, summary):
        with pytest.raises(IndexError):
            summary.shear_results(index=0)

    def test_shear_results_has_dcr(self, summary):
        df = summary.shear_results(index=1)
        assert "DCR" in df.columns


# ------------------------------------------------------------------
# Export / Import
# ------------------------------------------------------------------


class TestShearWallSummaryExportImport:
    def test_export_without_design_raises(self, summary):
        with pytest.raises(AttributeError, match="No design data"):
            summary.export_design("test.xlsx")

    def test_export_import_roundtrip(self, concrete, steel, sample_df_no_rebar, tmp_path):
        s = ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=sample_df_no_rebar)
        s.design()

        path = str(tmp_path / "walls.xlsx")
        s.export_design(path)
        assert os.path.exists(path)

        s.import_design(path)
        assert len(s.nodes) == 2  # 2 walls: (Level 1, M1) and (Level 2, M1)

    def test_reimport_check_matches(self, concrete, steel, sample_df_no_rebar, tmp_path):
        s = ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=sample_df_no_rebar)
        s.design()
        check_before = s.check()

        path = str(tmp_path / "walls.xlsx")
        s.export_design(path)
        s.import_design(path)
        check_after = s.check()

        dcr_before = check_before["DCR"].iloc[1:].tolist()
        dcr_after = check_after["DCR"].iloc[1:].tolist()
        for a, b in zip(dcr_before, dcr_after):
            assert abs(a - b) < 0.01


# ------------------------------------------------------------------
# Word export
# ------------------------------------------------------------------


class TestShearWallSummaryDoc:
    def test_results_detailed_doc(self, summary, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        summary.results_detailed_doc(index=1)
        expected_file = tmp_path / f"Shear_Wall_Summary_{summary.concrete.design_code}.docx"
        assert expected_file.exists()

    def test_results_detailed_doc_index_out_of_range(self, summary):
        with pytest.raises(IndexError):
            summary.results_detailed_doc(index=99)


# ------------------------------------------------------------------
# Level attribute
# ------------------------------------------------------------------


class TestShearWallLevel:
    def test_level_stored_on_wall(self, summary):
        for i, node in enumerate(summary.nodes):
            wall = node.section
            expected_level = summary.wall_keys[i][0]
            assert wall.level == expected_level


# ------------------------------------------------------------------
# Coverage: unrecognized unit (line 70)
# ------------------------------------------------------------------


class TestShearWallSummaryGetUnitVariable:
    def test_unrecognized_unit_raises(self, concrete, steel):
        data = {
            "Level": ["", "Level 1"],
            "Label": ["", "M1"],
            "Comb.": ["", "ELU 1"],
            "t": ["furlongs", 20],
            "lw": ["m", 3.0],
            "hw": ["m", 3.0],
            "cc": ["mm", 25],
            "Nx": ["kN", 0],
            "Vz": ["kN", 264],
            "My": ["kNm", -172],
            "dbh": ["mm", 8],
            "sh": ["cm", 20],
            "dbv": ["mm", 12],
            "sv": ["cm", 15],
        }
        df = pd.DataFrame(data)
        with pytest.raises(ValueError, match="Invalid unit"):
            ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=df)

    def test_get_unit_variable_unrecognized(self, summary):
        with pytest.raises(ValueError, match="not recognized"):
            summary.get_unit_variable("parsecs")


# ------------------------------------------------------------------
# Coverage: missing vertical rebar only (line 162)
# ------------------------------------------------------------------


class TestShearWallSummaryCheckVerticalMissing:
    def test_check_raises_missing_vertical_rebar(self, concrete, steel):
        data = {
            "Level": ["", "Level 1"],
            "Label": ["", "M1"],
            "Comb.": ["", "ELU 1"],
            "t": ["cm", 20],
            "lw": ["m", 3.0],
            "hw": ["m", 3.0],
            "cc": ["mm", 25],
            "Nx": ["kN", 0],
            "Vz": ["kN", 264],
            "My": ["kNm", -172],
            "dbh": ["mm", 8],
            "sh": ["cm", 20],
            "dbv": ["mm", 0],
            "sv": ["cm", 0],
        }
        df = pd.DataFrame(data)
        s = ShearWallSummary(concrete=concrete, steel_bar=steel, wall_list=df)
        with pytest.raises(ValueError, match="no vertical rebar"):
            s.check()


# ------------------------------------------------------------------
# Coverage: imperial unit branch (line 205)
# ------------------------------------------------------------------


class TestShearWallSummaryImperial:
    def test_check_imperial_units(self):
        concrete_imp = Concrete_ACI_318_19(name="C4000", f_c=4000 * psi)
        steel_imp = SteelBar(name="G60", f_y=60 * ksi)
        data = {
            "Level": ["", "Level 1"],
            "Label": ["", "W1"],
            "Comb.": ["", "ELU 1"],
            "t": ["inch", 10],
            "lw": ["ft", 12],
            "hw": ["ft", 10],
            "cc": ["inch", 1.5],
            "Nx": ["kN", 0],
            "Vz": ["kN", 200],
            "My": ["kNm", 0],
            "dbh": ["inch", 0.5],
            "sh": ["inch", 8],
            "dbv": ["inch", 0.5],
            "sv": ["inch", 8],
        }
        df = pd.DataFrame(data)
        s = ShearWallSummary(concrete=concrete_imp, steel_bar=steel_imp, wall_list=df)
        check_df = s.check()
        # Units row should show "kip" for imperial
        assert check_df["Vu,max"].iloc[0] == "kip"
        assert check_df["ØVn"].iloc[0] == "kip"
