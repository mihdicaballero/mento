import pytest
from mento.units import mm

from mento.settings import Settings  # Replace 'your_module' with the actual module name

def test_default_settings():
    settings = Settings()
    assert settings.get_setting('clear_cover') == 25 * mm #type: ignore
    assert settings.get_setting('clear_spacing') == 20 * mm #type: ignore
    assert settings.get_setting('stirrup_diameter') == 6 * mm #type: ignore
    assert settings.get_setting('vibrator_size') == 30 * mm #type: ignore
    assert settings.get_setting('layers_spacing') == 25 * mm #type: ignore

def test_custom_settings():
    custom_settings = {
        'clear_cover': 30 * mm, #type: ignore
        'stirrup_diameter': 8 * mm, #type: ignore
    }
    settings = Settings(settings_dict=custom_settings)
    assert settings.get_setting('clear_cover') == 30 * mm #type: ignore
    assert settings.get_setting('stirrup_diameter') == 8 * mm #type: ignore
    assert settings.get_setting('clear_spacing') == 20 * mm  #type: ignore

def test_get_setting_existing_key():
    settings = Settings()
    assert settings.get_setting('clear_cover') == 25 * mm #type: ignore

def test_get_setting_non_existent_key():
    settings = Settings()
    with pytest.raises(KeyError):
        settings.get_setting('non_existent_setting')

def test_set_setting_existing_key():
    settings = Settings()
    settings.set_setting('clear_cover', 35 * mm) #type: ignore
    assert settings.get_setting('clear_cover') == 35 * mm #type: ignore

def test_set_setting_non_existent_key():
    settings = Settings()
    with pytest.raises(KeyError):
        settings.set_setting('non_existent_setting', 50 * mm) #type: ignore
