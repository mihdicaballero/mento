import pytest
from mento.material import SteelBar, Concrete
from mento.units import MPa
from mento.section import Section     # Assuming this is the module name
from mento.settings import Settings
from typing import Dict, Any

@pytest.fixture()
def setup_basic_section() -> Section:
    concrete = Concrete('C25')
    steel_bar = SteelBar(name='B500S', f_y=500 * MPa)
    section = Section(label="V101",concrete=concrete, steel_bar=steel_bar)
    return section


def test_section_id_increment(setup_basic_section: Section) -> None:
    section1 = setup_basic_section
    section2 = Section(label="V102",concrete=section1.concrete, steel_bar=section1.steel_bar)
    assert section2.id == section1.id + 1


def test_initial_settings(setup_basic_section: Section) -> None:
    section = setup_basic_section
    assert isinstance(section.settings, Settings)


def test_update_settings(setup_basic_section: Section) -> None:
    section = setup_basic_section
    new_settings: Dict[str, Any] = {'clear_cover': 20}
    section.update_settings(new_settings)
    assert section.settings.get_setting('clear_cover') == 20


def test_get_settings(setup_basic_section: Section) -> None:
    section = setup_basic_section
    settings = section.get_settings
    assert isinstance(settings, dict)
    assert 'clear_cover' in settings
