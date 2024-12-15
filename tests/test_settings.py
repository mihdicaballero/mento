import pytest
from typing import Dict, Any

from mento.units import mm, inch, psi, MPa

from mento.settings import Settings  # Replace 'your_module' with the actual module name
from mento.material import Concrete

@pytest.mark.parametrize("unit_system, expected_defaults", [
    ("metric", {
        'clear_cover': 25 * mm, 
        'clear_spacing': 25 * mm, 
        'stirrup_diameter_ini': 8 * mm,
        'vibrator_size': 30 * mm, 
        'layers_spacing': 25 * mm,
        'max_diameter_diff': 5*mm,
        'max_bars_per_layer': 5,
        'minimum_longitudinal_diameter': 10*mm,
    }),
    ("imperial", {
        'clear_cover': 1 * inch, 
        'clear_spacing': 1 * inch, 
        'stirrup_diameter_ini': 3/8*inch,
        'vibrator_size': 1.25*inch, 
        'layers_spacing': 1*inch,
        'max_diameter_diff': 2/8*inch,
        'max_bars_per_layer': 5,
        'minimum_longitudinal_diameter': 4/8*inch,
    })
])
def test_default_settings(unit_system: str, expected_defaults: Dict[str, Any]) -> None:
    # Create a concrete object with the desired unit system
    concrete = Concrete(name='conc', f_c=25*MPa if unit_system == "metric" else 3000*psi)
    settings = Settings(concrete=concrete)

    # Verify that all expected defaults match
    for key, expected_value in expected_defaults.items():
        assert settings.get_setting(key) == expected_value, f"Mismatch for {key} in {unit_system} system"

def test_get_setting_existing_key() -> None:
    settings = Settings()
    assert settings.get_setting('clear_cover') == 25 * mm

def test_custom_settings() -> None:
    custom_settings = {
        'clear_cover': 30 * mm,
        'stirrup_diameter': 8 * mm,
    }
    settings = Settings(settings=custom_settings)
    assert settings.get_setting('clear_cover') == 30 * mm
    assert settings.get_setting('stirrup_diameter') == 8 * mm
    assert settings.get_setting('clear_spacing') == 25 * mm 
    
def test_get_setting_non_existent_key() -> None:
    settings = Settings()
    with pytest.raises(KeyError):
        settings.get_setting('non_existent_setting')

def test_set_setting_existing_key() -> None:
    settings = Settings()
    settings.set_setting('clear_cover', 35 * mm)
    assert settings.get_setting('clear_cover') == 35 * mm

def test_set_setting_non_existent_key() -> None:
    settings = Settings()
    with pytest.raises(KeyError):
        settings.set_setting('non_existent_setting', 50 * mm)
