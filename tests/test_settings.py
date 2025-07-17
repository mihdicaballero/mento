import pytest
from dataclasses import fields
from pint import Quantity

from mento.settings import BeamSettings, GLOBAL_BEAM_SETTINGS
from mento.units import mm, inch


# Helper to check if a quantity matches another, allowing for slight float differences
def assert_quantity_equal(q1: Quantity, q2: Quantity, rtol: float = 1e-9) -> None:
    """Asserts two Pint Quantities are numerically and unit-wise equal."""
    assert q1.units == q2.units, f"Units mismatch: {q1.units} vs {q2.units}"
    assert (
        pytest.approx(q1.magnitude, rel=rtol) == q2.magnitude
    ), f"Magnitudes mismatch: {q1.magnitude} vs {q2.magnitude}"


# --- Fixtures ---


@pytest.fixture
def metric_settings() -> BeamSettings:
    """Fixture for default metric BeamSettings."""
    return BeamSettings(unit_system="metric")


@pytest.fixture
def imperial_settings() -> BeamSettings:
    """Fixture for default imperial BeamSettings."""
    return BeamSettings(unit_system="imperial")


# --- Tests for Default Values ---


def test_metric_defaults_all_set(metric_settings: BeamSettings) -> None:
    """Test that all metric defaults are correctly applied when no custom values are given."""
    expected_defaults = (
        BeamSettings._metric_defaults
    )  # Access the class variable for comparison

    for field_info in fields(BeamSettings):
        if field_info.name in expected_defaults:  # Check only fields with defaults
            actual_value = getattr(metric_settings, field_info.name)
            expected_value = expected_defaults[field_info.name]

            if isinstance(expected_value, Quantity):
                assert_quantity_equal(actual_value, expected_value)
            else:  # For max_bars_per_layer (int)
                assert actual_value == expected_value

    assert metric_settings.unit_system == "metric"


def test_imperial_defaults_all_set(imperial_settings: BeamSettings) -> None:
    """Test that all imperial defaults are correctly applied when no custom values are given."""
    expected_defaults = BeamSettings._imperial_defaults

    for field_info in fields(BeamSettings):
        if field_info.name in expected_defaults:
            actual_value = getattr(imperial_settings, field_info.name)
            expected_value = expected_defaults[field_info.name]

            if isinstance(expected_value, Quantity):
                assert_quantity_equal(actual_value, expected_value)
            else:
                assert actual_value == expected_value

    assert imperial_settings.unit_system == "imperial"


def test_override_specific_metric_setting() -> None:
    """Test overriding a specific setting in metric system."""
    custom_clear_spacing = 30 * mm
    settings = BeamSettings(unit_system="metric", clear_spacing=custom_clear_spacing)
    assert_quantity_equal(settings.clear_spacing, custom_clear_spacing)
    # Ensure other settings remain default metric
    assert_quantity_equal(
        settings.stirrup_diameter_ini,
        BeamSettings._metric_defaults["stirrup_diameter_ini"],
    )


def test_override_specific_imperial_setting() -> None:
    """Test overriding a specific setting in imperial system."""
    custom_stirrup_diameter = 0.5 * inch
    settings = BeamSettings(
        unit_system="imperial", stirrup_diameter_ini=custom_stirrup_diameter
    )
    assert_quantity_equal(settings.stirrup_diameter_ini, custom_stirrup_diameter)
    # Ensure other settings remain default imperial
    assert_quantity_equal(
        settings.clear_spacing, BeamSettings._imperial_defaults["clear_spacing"]
    )


def test_custom_unit_system_with_no_defaults_set() -> None:
    """Test that if a unit system is specified without a corresponding default dict, it uses metric."""
    # Temporarily rename _metric_defaults to simulate a missing imperial_defaults scenario
    # (though in your code, imperial_defaults is always present).
    # This scenario would require a more complex mock/patch, or testing how it handles an unknown 'unit_system' string.
    # Given the current implementation, it will always fall back to metric if not "imperial".
    settings = BeamSettings(unit_system="unknown_system")
    assert settings.unit_system == "unknown_system"  # Unit system itself is set
    assert_quantity_equal(
        settings.clear_spacing, BeamSettings._metric_defaults["clear_spacing"]
    )
    # This implies that any system other than "imperial" defaults to metric. This is good to test.


# --- Tests for Validation ---


def test_max_bars_per_layer_validation_zero() -> None:
    """Test that max_bars_per_layer must be at least 1 (zero case)."""
    with pytest.raises(ValueError, match="max_bars_per_layer must be at least 1"):
        BeamSettings(unit_system="metric", max_bars_per_layer=0)


def test_max_bars_per_layer_validation_negative() -> None:
    """Test that max_bars_per_layer must be at least 1 (negative case)."""
    with pytest.raises(ValueError, match="max_bars_per_layer must be at least 1"):
        BeamSettings(unit_system="metric", max_bars_per_layer=-5)


def test_max_bars_per_layer_valid() -> None:
    """Test that a valid max_bars_per_layer does not raise an error."""
    try:
        settings = BeamSettings(unit_system="metric", max_bars_per_layer=1)
        assert settings.max_bars_per_layer == 1
        settings = BeamSettings(unit_system="metric", max_bars_per_layer=5)  # Default
        assert settings.max_bars_per_layer == 5
    except ValueError as e:
        pytest.fail(f"Valid max_bars_per_layer raised an error: {e}")


# --- Tests for __str__ method ---


def test_str_metric_defaults(metric_settings: BeamSettings) -> None:
    """Test the string representation for default metric settings."""
    s = str(metric_settings)
    expected_lines = [
        "clear_spacing: 25.00 mm",
        "stirrup_diameter_ini: 8.00 mm",
        "vibrator_size: 30.00 mm",
        "layers_spacing: 25.00 mm",
        "max_diameter_diff: 5.00 mm",
        "minimum_longitudinal_diameter: 8.00 mm",
        "max_bars_per_layer: 5",
    ]
    for line in expected_lines:
        assert line in s
    # Ensure no unit_system in output
    assert "unit_system" not in s


def test_str_imperial_defaults(imperial_settings: BeamSettings) -> None:
    """Test the string representation for default imperial settings."""
    s = str(imperial_settings)
    # Note: Pint's default for 3/8 inch might be 0.375 inch or similar.
    # Using 'in' for inches is common with '~' format.
    expected_lines = [
        "clear_spacing: 1.00 in",
        "stirrup_diameter_ini: 0.38 in",  # 3/8 is 0.375, .2f will round
        "vibrator_size: 1.25 in",
        "layers_spacing: 1.00 in",
        "max_diameter_diff: 0.25 in",
        "minimum_longitudinal_diameter: 0.38 in",
        "max_bars_per_layer: 5",
    ]
    for line in expected_lines:
        assert line in s
    assert "unit_system" not in s


def test_str_with_custom_values() -> None:
    """Test the string representation with some custom values."""
    settings = BeamSettings(
        unit_system="metric", clear_spacing=40 * mm, max_bars_per_layer=7
    )
    s = str(settings)
    assert "clear_spacing: 40.00 mm" in s
    assert "max_bars_per_layer: 7" in s
    # Check that other defaults are still present
    assert "stirrup_diameter_ini: 8.00 mm" in s
    assert "unit_system" not in s


# --- Tests for GLOBAL_BEAM_SETTINGS ---


def test_global_beam_settings_exists() -> None:
    """Test that GLOBAL_BEAM_SETTINGS is an instance of BeamSettings."""
    assert isinstance(GLOBAL_BEAM_SETTINGS, BeamSettings)


def test_global_beam_settings_is_metric_by_default(
    metric_settings: BeamSettings,
) -> None:
    """Test that GLOBAL_BEAM_SETTINGS uses metric defaults."""
    # Compare with a fresh default metric settings instance
    assert GLOBAL_BEAM_SETTINGS == metric_settings
    # Note: Dataclass equality works by comparing all fields.
    # If any internal _id or _last_id logic affects equality, this might need adjustment.
    # Given your BeamSettings is immutable and no _id logic, direct comparison is fine.
    assert GLOBAL_BEAM_SETTINGS.unit_system == "metric"
    assert_quantity_equal(GLOBAL_BEAM_SETTINGS.clear_spacing, 25 * mm)
