import sys
import os

# Add the project root to PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
from structurelab.material import create_concrete
from structurelab.settings import Settings
from structurelab.concrete.beam import Beam
from structurelab.rebar import Rebar
import forallpeople

forallpeople.environment('structural', top_level=True)
cm = 1e-2 * forallpeople.m  # type: ignore

concrete = create_concrete(name="H30", f_c=30*MPa, design_code="ACI 318-19") # type: ignore
beam = Beam(
    name="V 20x50",
    concrete=concrete,
    steelBar="Barras Longitudinales",
    width=20*cm,
    depth=50*cm,
)
custom_settings=Settings({
        'clear_cover': 30 * mm, # type: ignore
        'stirrup_diameter': 8 * mm, # type: ignore
    })

# Tests
def test_rebar_initialization(section, custom_settings):
    as_nec = 15 * cm**2
    rebar = Rebar(required_as=as_nec, beam=section, settings=custom_settings)

    assert rebar.required_as == as_nec
    assert rebar.beam == section
    assert rebar.settings == custom_settings

def test_rebar_calculate_rebars(section, custom_settings):
    as_nec = 15 * cm**2
    rebar = Rebar(required_as=as_nec, beam=section, settings=custom_settings)
    best_combination = rebar.calculate_rebars()

    assert best_combination is not None
    assert 'layer_1' in best_combination
    assert 'num_bars_1' in best_combination
    assert 'diameter_1' in best_combination
    assert 'total_as' in best_combination
    assert 'DCR' in best_combination
    assert 'available_spacing_1' in best_combination

    # Check if total_as meets or exceeds the required area
    assert best_combination['total_as'] >= as_nec

def test_rebar_no_solution_possible():
    # Adjust the settings to a scenario where no solution is possible
    impossible_settings = Settings({
        'clear_cover': 30 * mm,  # type: ignore
        'stirrup_diameter': 32 * mm,  # type: ignore
        'clear_spacing': 2 * mm  # type: ignore
    })

    concrete = create_concrete(name="H30", f_c=30*MPa, design_code="ACI 318-19")  # type: ignore
    section = Beam(
        name="V 20x50",
        concrete=concrete,
        steelBar="Barras Longitudinales",
        width=20*cm,
        depth=50*cm,
    )

    rebar = Rebar(required_as=100 * cm**2, beam=section, settings=impossible_settings)

    with pytest.raises(ValueError, match="Cannot fit the required reinforcement"):
        rebar.calculate_rebars()

# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()