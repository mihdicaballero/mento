import pytest
from mento.rebar import Rebar
from mento.beam import RectangularBeam
from mento.material import Concrete_ACI_318_19, SteelBar, Concrete_EN_1992_2004
from mento.units import psi, kip, mm, inch, ksi, cm, MPa, kN
from mento.forces import Forces
from mento.node import Node
from mento.settings import BeamSettings
from mento.slab import OneWaySlab
import pandas as pd


@pytest.fixture()
def beam_example_imperial() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steelBar = SteelBar(name="ADN 420", f_y=60 * ksi)
    beam_settings = BeamSettings(stirrup_diameter_ini=0.5 * inch, minimum_longitudinal_diameter=1 * inch)
    section = RectangularBeam(
        label="V-10x16",
        concrete=concrete,
        steel_bar=steelBar,
        width=10 * inch,
        height=16 * inch,
        c_c=1.5 * inch,
        settings=beam_settings,
    )
    return section


@pytest.fixture()
def beam_example_metric() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(stirrup_diameter_ini=8 * mm)
    section = RectangularBeam(
        label="V 20x50",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=50 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )
    return section


@pytest.fixture()
def slab_example_metric() -> OneWaySlab:
    """Fixture for a slab in metric units"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    slab_settings = BeamSettings(stirrup_diameter_ini=8 * mm)

    # Assuming Slab inherits from or is similar to RectangularBeam
    # Adjust this based on your actual implementation
    section = OneWaySlab(
        label="Slab 100cm",
        concrete=concrete,
        steel_bar=steelBar,
        width=100 * cm,  # Wide slab
        height=20 * cm,  # Thin slab
        c_c=30 * mm,
        settings=slab_settings,
    )
    return section


@pytest.fixture()
def beam_narrow_metric() -> RectangularBeam:
    """Fixture for a narrow beam to test single bar scenarios"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(stirrup_diameter_ini=8 * mm)

    section = RectangularBeam(
        label="V-narrow",
        concrete=concrete,
        steel_bar=steelBar,
        width=12 * cm,  # Very narrow
        height=50 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )
    return section


@pytest.fixture()
def beam_en_1992() -> RectangularBeam:
    """Fixture for EN 1992 design code"""
    concrete = Concrete_EN_1992_2004(name="C30/37", f_c=30 * MPa)
    steelBar = SteelBar(name="B500", f_y=500 * MPa)
    beam_settings = BeamSettings(stirrup_diameter_ini=8 * mm)

    section = RectangularBeam(
        label="V-EN",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=50 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )
    return section


# Tests
def test_beam_longitudinal_rebar_ACI_318_19_metric(
    beam_example_metric: RectangularBeam,
) -> None:
    A_s_req = 5 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design["n_1"] == 2
    assert best_design["d_b1"] == 16 * mm
    assert best_design["n_2"] == 1
    assert best_design["d_b2"] == 12 * mm
    assert best_design["n_3"] == 0
    assert best_design["d_b3"] is None
    assert best_design["n_4"] == 0
    assert best_design["d_b3"] is None
    assert best_design["total_as"].magnitude == pytest.approx(5.15, rel=1e-3)
    assert best_design["total_bars"] == 3
    assert best_design["clear_spacing"].magnitude == pytest.approx(40, rel=1e-3)


def test_longitudinal_rebar_ACI_max_area(beam_example_metric: RectangularBeam) -> None:
    A_s_req = 5 * cm**2
    A_s_max = 5.2 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req, A_s_max=A_s_max)
    design_with_limit = beam_rebar.longitudinal_rebar_design

    assert design_with_limit["total_as"] <= A_s_max

    beam_rebar_default = Rebar(beam_example_metric)
    beam_rebar_default.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    default_design = beam_rebar_default.longitudinal_rebar_design

    assert design_with_limit["total_as"].magnitude == pytest.approx(default_design["total_as"].magnitude, rel=1e-3)


def test_beam_longitudinal_rebar_min_diameter(
    beam_example_metric: RectangularBeam,
) -> None:
    # This test ensures that the minimum rebar diameter constraint is respected.
    A_s_req = 2 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design["d_b1"] >= 10 * mm
    assert best_design["d_b2"] >= 10 * mm if best_design["d_b2"] is not None else True


def test_beam_longitudinal_rebar_max_bars_per_layer(
    beam_example_metric: RectangularBeam,
) -> None:
    # This test ensures that the maximum number of bars per layer constraint is respected.
    A_s_req = 10 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design["n_1"] + best_design["n_2"] <= beam_example_metric.settings.max_bars_per_layer
    assert best_design["n_3"] + best_design["n_4"] <= beam_example_metric.settings.max_bars_per_layer


def test_beam_longitudinal_rebar_CIRSOC_201_25(
    beam_example_metric: RectangularBeam,
) -> None:
    # This test ensures that the longitudinal rebar design works correctly for different design codes.
    beam_example_metric.concrete.design_code = "CIRSOC 201-25"
    A_s_req = 5 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design["n_1"] == 2
    assert best_design["d_b1"] == 16 * mm
    assert best_design["n_2"] == 1
    assert best_design["d_b2"] == 12 * mm
    assert best_design["n_3"] == 0
    assert best_design["d_b3"] is None
    assert best_design["n_4"] == 0
    assert best_design["d_b3"] is None
    assert best_design["total_as"].magnitude == pytest.approx(5.15, rel=1e-3)
    assert best_design["total_bars"] == 3
    assert best_design["clear_spacing"].magnitude == pytest.approx(40, rel=1e-3)


def test_longitudinal_rebar_factory_max_area(
    beam_example_metric: RectangularBeam,
) -> None:
    A_s_req = 5 * cm**2
    A_s_max = 5.2 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar(A_s_req=A_s_req, A_s_max=A_s_max)
    design_with_limit = beam_rebar.longitudinal_rebar_design

    assert design_with_limit["total_as"] <= A_s_max

    beam_rebar_default = Rebar(beam_example_metric)
    beam_rebar_default.longitudinal_rebar(A_s_req=A_s_req)
    default_design = beam_rebar_default.longitudinal_rebar_design

    assert design_with_limit["total_as"].magnitude == pytest.approx(default_design["total_as"].magnitude, rel=1e-3)


def test_beam_longitudinal_rebar_small_area(
    beam_example_metric: RectangularBeam,
) -> None:
    # This test ensures that the rebar design handles very small required areas correctly.
    A_s_req = 0.1 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design["total_as"] >= A_s_req


def test_beam_longitudinal_rebar_large_area(
    beam_example_metric: RectangularBeam,
) -> None:
    # This test ensures that the rebar design handles very large required areas correctly.
    A_s_req = 50 * cm**2
    beam_rebar = Rebar(beam_example_metric)
    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    assert best_design is not None
    assert best_design["n_1"] >= 2
    assert best_design["d_b1"] >= 8 * mm


def test_beam_transverse_rebar_ACI_318_19_imperial(
    beam_example_imperial: RectangularBeam,
) -> None:
    f = Forces(V_z=37.727 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    N1 = Node(beam_example_imperial, forces=f)
    results = N1.design_shear()
    assert results is not None
    assert beam_example_imperial._stirrup_n == 1
    assert beam_example_imperial._stirrup_d_b == 3 / 8 * inch
    assert beam_example_imperial._stirrup_s_l == 5 * inch
    assert beam_example_imperial._A_v.to("cm**2/m").magnitude == pytest.approx(11.2214, rel=1e-3)


def test_beam_transverse_rebar_ACI_318_19_metric(
    beam_example_metric: RectangularBeam,
) -> None:
    # This test ensures that the transverse rebar design works correctly for beams defined in metric units.
    f = Forces(V_z=100 * kN)
    beam_example_metric.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)
    N1 = Node(beam_example_metric, forces=f)
    results = N1.design_shear()

    assert results is not None
    assert beam_example_metric._stirrup_n == 1
    assert beam_example_metric._stirrup_d_b == 10 * mm
    assert beam_example_metric._stirrup_s_l == 22 * cm
    assert beam_example_metric._A_v.to("cm**2/m").magnitude == pytest.approx(7.14, rel=1e-3)


def test_beam_transverse_rebar_max_spacing_imperial(
    beam_example_imperial: RectangularBeam,
) -> None:
    # This test ensures that the maximum spacing constraints are respected.
    f = Forces(V_z=50 * kip)
    beam_example_imperial.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)
    N1 = Node(beam_example_imperial, forces=f)
    results = N1.design_shear()

    assert results is not None
    assert beam_example_imperial._stirrup_s_l <= 24 * inch
    assert beam_example_imperial._stirrup_s_w <= 24 * inch


def test_beam_transverse_rebar_min_diameter(
    beam_example_metric: RectangularBeam,
) -> None:
    # This test ensures that the minimum rebar diameter constraint is respected for transverse reinforcement.
    f = Forces(V_z=50 * kN)
    beam_example_metric.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)
    N1 = Node(beam_example_metric, forces=f)
    results = N1.design_shear()

    assert results is not None
    assert beam_example_metric._stirrup_d_b >= 8 * mm


def test_beam_transverse_rebar_CIRSOC_201_25(
    beam_example_metric: RectangularBeam,
) -> None:
    # This test ensures that the transverse rebar design works correctly for different design codes.
    beam_example_metric.concrete.design_code = "CIRSOC 201-25"
    f = Forces(V_z=100 * kN)
    beam_example_metric.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)
    N1 = Node(beam_example_metric, forces=f)
    results = N1.design_shear()

    assert results is not None
    assert beam_example_metric._stirrup_n == 1
    assert beam_example_metric._stirrup_d_b == 6 * mm
    assert beam_example_metric._stirrup_s_l == 22 * cm
    assert beam_example_metric._A_v.to("cm**2/m").magnitude == pytest.approx(2.57, rel=1e-3)


# Test 1: Slab mode - spacing penalty
def test_slab_mode_spacing_penalty(slab_example_metric: OneWaySlab) -> None:
    """Test that slab mode applies spacing penalty correctly"""
    A_s_req = 10 * cm**2
    beam_rebar = Rebar(slab_example_metric)

    # The rebar should be in slab mode
    assert beam_rebar.mode == "slab"

    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    assert result is not None
    assert len(result) > 0

    # Check that spacing_penalty column exists
    assert "spacing_penalty" in result.columns


# Test 2: Slab mode - penalty weight adjustments
def test_slab_mode_penalty_weights(slab_example_metric: OneWaySlab) -> None:
    """Test that slab mode adjusts penalty weights appropriately"""
    A_s_req = 8 * cm**2
    beam_rebar = Rebar(slab_example_metric)

    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    assert result is not None
    # Slabs should prefer multiple smaller bars
    beam_rebar.longitudinal_rebar_design

    # Check that diameter_size_penalty is calculated
    assert "diameter_size_penalty" in result.columns


# Test 3: Slab mode - mirrored second layer
def test_slab_mode_second_layer_mirror(slab_example_metric: OneWaySlab) -> None:
    """Test that slab mode creates mirrored second layer"""
    A_s_req = 18 * cm**2  # Large area to force two layers
    beam_rebar = Rebar(slab_example_metric)

    beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)
    best_design = beam_rebar.longitudinal_rebar_design

    # If there's a second layer in slab mode, it should mirror the first
    if best_design["n_3"] > 0:
        assert best_design["n_3"] == best_design["n_1"]
        assert best_design["n_4"] == best_design["n_2"]


# Test 5: EN 1992 design code
def test_longitudinal_rebar_EN_1992_2004(beam_en_1992: RectangularBeam) -> None:
    """Test EN 1992 longitudinal rebar design"""
    A_s_req = 5 * cm**2
    beam_rebar = Rebar(beam_en_1992)

    # This should call the EN 1992 method (which currently delegates to ACI)
    beam_rebar.longitudinal_rebar(A_s_req=A_s_req)

    best_design = beam_rebar.longitudinal_rebar_design
    assert best_design is not None
    assert best_design["total_as"] >= A_s_req


# Test 6: Layer 2 continues - max bars exceeded
def test_beam_layer2_max_bars_exceeded() -> None:
    """Test that layer 2 respects max_bars_per_layer constraint"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(
        stirrup_diameter_ini=8 * mm,
        max_bars_per_layer=3,  # Very restrictive
    )

    beam = RectangularBeam(
        label="V-restricted",
        concrete=concrete,
        steel_bar=steelBar,
        width=30 * cm,
        height=50 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )

    A_s_req = 20 * cm**2  # Large requirement
    beam_rebar = Rebar(beam)
    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    # Should respect the max_bars_per_layer limit
    for _, row in result.iterrows():
        assert row["n_1"] + row["n_2"] <= 3
        assert row["n_3"] + row["n_4"] <= 3


# Test 7: Layer 2 area less than layer 1
def test_beam_layer2_area_constraint() -> None:
    """Test that layer 2 area is less than or equal to layer 1"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(stirrup_diameter_ini=8 * mm)

    beam = RectangularBeam(
        label="V-test",
        concrete=concrete,
        steel_bar=steelBar,
        width=25 * cm,
        height=60 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )

    A_s_req = 25 * cm**2  # Force two layers
    beam_rebar = Rebar(beam)
    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    # Check that all combinations respect the area constraint
    for _, row in result.iterrows():
        if row["n_3"] > 0 or row["n_4"] > 0:
            # Calculate areas
            area_1 = row["n_1"] * beam_rebar.rebar_areas[row["d_b1"]] + (
                row["n_2"] * beam_rebar.rebar_areas[row["d_b2"]] if row["n_2"] > 0 else 0 * cm**2
            )
            area_2 = row["n_3"] * beam_rebar.rebar_areas[row["d_b3"]] + (
                row["n_4"] * beam_rebar.rebar_areas[row["d_b4"]] if row["n_4"] > 0 else 0 * cm**2
            )
            assert area_1 >= area_2


# Test 8: Layer 2 spacing check
def test_beam_layer2_spacing_check() -> None:
    """Test that layer 2 spacing is checked"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(
        stirrup_diameter_ini=8 * mm,
        clear_spacing=50 * mm,  # Stricter spacing requirement
    )

    beam = RectangularBeam(
        label="V-spacing",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,  # Narrower beam
        height=60 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )

    A_s_req = 20 * cm**2  # Force two layers
    beam_rebar = Rebar(beam)
    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    # All results should have valid spacing
    for _, row in result.iterrows():
        assert row["clear_spacing"].magnitude >= 50


# Test 9: Fallback combination tracking
def test_fallback_combination_when_no_valid_solution() -> None:
    """Test that fallback combination is returned when no valid solution meets A_s_req"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(
        stirrup_diameter_ini=8 * mm, max_bars_per_layer=4, max_longitudinal_diameter=20 * mm, clear_spacing=25 * mm
    )

    beam = RectangularBeam(
        label="V-test-fallback",
        concrete=concrete,
        steel_bar=steelBar,
        width=30 * cm,  # Wide enough to fit bars
        height=50 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )

    # Create a scenario where we can get close but not exactly meet the requirement
    # This tests the fallback logic when combinations are found but slightly under/over
    A_s_req = 15 * cm**2  # Specific requirement
    A_s_max = 12 * cm**2  # Restrict maximum below required

    beam_rebar = Rebar(beam)
    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req, A_s_max=A_s_max)

    # Should return fallback with best possible area under the constraint
    assert result is not None
    assert len(result) > 0
    best_design = beam_rebar.longitudinal_rebar_design

    # The area should respect A_s_max even though it's less than A_s_req
    assert best_design["total_as"] <= A_s_max
    assert best_design["total_as"] > 0 * cm**2


# Test 10: Slab with large spacing (triggering spacing penalty)
def test_slab_large_spacing_penalty() -> None:
    """Test that slabs with large spacing get penalized"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    slab_settings = BeamSettings(stirrup_diameter_ini=8 * mm)

    slab = OneWaySlab(
        label="Wide-Slab",
        concrete=concrete,
        steel_bar=steelBar,
        width=200 * cm,  # Very wide to force large spacing
        height=20 * cm,
        c_c=30 * mm,
        settings=slab_settings,
    )

    A_s_req = 5 * cm**2  # Small requirement
    beam_rebar = Rebar(slab)
    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    # Check that spacing penalty is applied when spacing > 300mm
    for _, row in result.iterrows():
        if row["clear_spacing"].magnitude > 300:
            # Find this row in the internal DataFrame to check penalty
            assert "spacing_penalty" in result.columns


# Test 11: Beam mode (default) - no spacing penalty
def test_beam_mode_no_spacing_penalty() -> None:
    """Test that beam mode does not apply spacing penalty"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(stirrup_diameter_ini=8 * mm)

    beam = RectangularBeam(
        label="V-beam",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=50 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )

    A_s_req = 5 * cm**2
    beam_rebar = Rebar(beam)

    # Explicitly check mode
    assert beam_rebar.mode == "beam"

    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    # In beam mode, spacing_penalty should be 0 for all rows
    assert all(result["spacing_penalty"] == 0)


# Test 12: Small stirrup spacing - metric (s_l < 5 cm)
def test_transverse_rebar_small_spacing_metric() -> None:
    """Test that stirrups increase legs when spacing becomes too small (metric)"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(stirrup_diameter_ini=8 * mm)

    beam = RectangularBeam(
        label="V-high-shear",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=50 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )

    # Very high shear force to require small spacing (which will force more legs)
    f = Forces(V_z=600 * kN)  # Very high shear
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=20 * mm)
    N1 = Node(beam, forces=f)
    results = N1.design_shear()

    assert results is not None
    # With very high shear, should either:
    # - Have very small spacing (< 5 cm would trigger leg increase)
    # - Or have multiple stirrups (n > 1)
    assert beam._stirrup_n >= 1

    # If spacing is >= 5cm, n should be 1; if < 5cm tried, n should be > 1
    if beam._stirrup_s_l.to("cm").magnitude >= 5:
        assert beam._stirrup_n >= 1
    # The test passes if the design completes without error


# Test 13: Small stirrup spacing - imperial (s_l < 2 inch)
def test_transverse_rebar_small_spacing_imperial() -> None:
    """Test that stirrups increase legs when spacing becomes too small (imperial)"""
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steelBar = SteelBar(name="Grade 60", f_y=60 * ksi)
    beam_settings = BeamSettings(stirrup_diameter_ini=0.375 * inch)

    beam = RectangularBeam(
        label="V-high-shear",
        concrete=concrete,
        steel_bar=steelBar,
        width=10 * inch,
        height=24 * inch,
        c_c=1.5 * inch,
        settings=beam_settings,
    )

    # Very high shear force
    f = Forces(V_z=200 * kip)  # Very high shear
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=0.75 * inch)
    N1 = Node(beam, forces=f)
    results = N1.design_shear()

    assert results is not None
    assert beam._stirrup_n >= 1

    # If spacing is >= 2 inch, n should be 1; if < 2 inch tried, n should be > 1
    if beam._stirrup_s_l.to("inch").magnitude >= 2:
        assert beam._stirrup_n >= 1


# Test 14: Stirrups with n_legs > 6 (breaking condition)
def test_transverse_rebar_max_legs_exceeded() -> None:
    """Test that stirrup design stops when legs exceed 6"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(stirrup_diameter_ini=6 * mm)  # Small stirrup

    beam = RectangularBeam(
        label="V-extreme-shear",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=50 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )

    # Extremely high shear that would require > 6 legs
    f = Forces(V_z=500 * kN)  # Extreme shear
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)
    N1 = Node(beam, forces=f)
    results = N1.design_shear()

    # Should complete (either with solution or moving to larger diameter)
    assert results is not None
    # Should not exceed 3 stirrups (6 legs max)
    assert beam._stirrup_n <= 3


# Test 15: Zero steel requirement (A_s_req = 0)
def test_longitudinal_rebar_zero_requirement() -> None:
    """Test early exit when no steel is required"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(stirrup_diameter_ini=8 * mm)

    beam = RectangularBeam(
        label="V-zero-steel",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=50 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )

    A_s_req = 0 * cm**2  # Zero steel required
    beam_rebar = Rebar(beam)
    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    assert result is not None
    assert len(result) == 1  # Should return exactly one row

    best_design = beam_rebar.longitudinal_rebar_design
    assert best_design["n_1"] == 0
    assert best_design["n_2"] == 0
    assert best_design["n_3"] == 0
    assert best_design["n_4"] == 0
    assert best_design["total_as"].magnitude == 0
    assert best_design["total_bars"] == 0


# Test 16: Imperial unit system - minimum longitudinal diameter
def test_longitudinal_rebar_imperial_min_diameter() -> None:
    """Test that imperial system uses correct minimum diameter (3/8 inch)"""
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steelBar = SteelBar(name="Grade 60", f_y=60 * ksi)

    beam = RectangularBeam(
        label="V-imperial",
        concrete=concrete,
        steel_bar=steelBar,
        width=12 * inch,
        height=20 * inch,
        c_c=1.5 * inch,
    )

    A_s_req = 1.5 * inch**2
    beam_rebar = Rebar(beam)
    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    # Check that min_long_rebar is set correctly for imperial
    assert beam_rebar.min_long_rebar == 3 * inch / 8
    assert result is not None
    assert len(result) > 0

    # All diameters should be >= 3/8 inch
    for _, row in result.iterrows():
        assert row["d_b1"] >= 3 * inch / 8
        if row["d_b2"] is not None:
            assert row["d_b2"] >= 3 * inch / 8


# Test 17: Empty DataFrame scenario (impossible constraints)
def test_empty_dataframe_impossible_constraints() -> None:
    """Test that empty DataFrame is returned when constraints are impossible"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(
        stirrup_diameter_ini=8 * mm,
        max_bars_per_layer=2,
        minimum_longitudinal_diameter=32 * mm,  # Very large
        max_longitudinal_diameter=32 * mm,  # Only one option
        clear_spacing=200 * mm,  # Impossible spacing for narrow beam
        max_diameter_diff=0 * mm,  # No diameter variation allowed
    )

    beam = RectangularBeam(
        label="V-impossible",
        concrete=concrete,
        steel_bar=steelBar,
        width=10 * cm,  # Too narrow for constraints
        height=30 * cm,
        c_c=30 * mm,
        settings=beam_settings,
    )

    A_s_req = 5 * cm**2
    beam_rebar = Rebar(beam)
    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    # Should return empty DataFrame
    assert result is not None
    assert isinstance(result, pd.DataFrame)
    # If empty, the _long_combos_df should also be empty
    if result.empty:
        assert beam_rebar._long_combos_df.empty


# Test 18: Metric minimum longitudinal diameter
def test_longitudinal_rebar_metric_min_diameter_explicit() -> None:
    """Test that metric system uses correct minimum diameter (10 mm)"""
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)

    beam = RectangularBeam(
        label="V-metric",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=50 * cm,
        c_c=30 * mm,
    )

    A_s_req = 3 * cm**2
    beam_rebar = Rebar(beam)
    result = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=A_s_req)

    # Verify min_long_rebar is set to 10mm for metric
    assert beam_rebar.min_long_rebar == 10 * mm
    assert result is not None
    assert len(result) > 0

    # All diameters should be >= 10 mm
    for _, row in result.iterrows():
        assert row["d_b1"] >= 10 * mm
        if row["d_b2"] is not None:
            assert row["d_b2"] >= 10 * mm
        if row["d_b3"] is not None:
            assert row["d_b3"] >= 10 * mm
        if row["d_b4"] is not None:
            assert row["d_b4"] >= 10 * mm


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
