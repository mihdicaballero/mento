import pytest
from mento.units import mm

from mento.settings import Settings  # Replace 'your_module' with the actual module name

def test_default_settings() -> None:
    settings = Settings()
    assert settings.get_setting('clear_cover') == 25 * mm
    assert settings.get_setting('clear_spacing') == 25 * mm
    assert settings.get_setting('stirrup_diameter_ini') == 10 * mm
    assert settings.get_setting('longitudinal_diameter_ini') == 16 * mm
    assert settings.get_setting('vibrator_size') == 30 * mm
    assert settings.get_setting('layers_spacing') == 25 * mm

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
