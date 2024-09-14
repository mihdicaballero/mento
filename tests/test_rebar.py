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
section = Beam(
    name="V 20x50",
    concrete=concrete,
    steelBar="Barras Longitudinales",
    width=20*cm,
    depth=50*cm,
)
custom_settings={
        'clear_cover': 30 * mm, # type: ignore
        'stirrup_diameter': 8 * mm, # type: ignore
    }
section.update_settings(custom_settings)

# Tests
def test_rebar_initialization():
    as_nec = 15 * cm**2
    rebar = Rebar(required_as=as_nec, beam=section)

    assert rebar.required_as == as_nec
    assert rebar.beam == section

def test_rebar_calculate_rebars():
    as_nec = 16.08 * cm**2
    rebar = Rebar(required_as=as_nec, beam=section)
    best_combination = rebar.calculate_rebars()

    assert best_combination is not None
    assert 'layer_1' in best_combination
    assert 'num_bars_1' in best_combination
    assert 'diameter_1' in best_combination
    assert 'total_as' in best_combination
    assert 'available_spacing_1' in best_combination

    # Check if total_as meets or exceeds the required area
    assert best_combination['total_as'].value == pytest.approx(as_nec.value,rel=1e-3)

def test_rebar_no_solution_possible():
    as_nec = 150 * cm**2
    rebar = Rebar(required_as=as_nec, beam=section)

    with pytest.raises(ValueError, match="Cannot fit the required reinforcement within the beam width considering clear cover and spacing."):
        rebar.calculate_rebars()
        raise ValueError('Cannot fit the required reinforcement within the beam width considering clear cover and spacing')

# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()