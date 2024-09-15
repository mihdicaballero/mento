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
def test_beam_longitudinal_ACI_318_19():
    as_nec = 5 * cm**2
    best_combination = Rebar(beam=section).beam_longitudinal_ACI_318_19(as_req=as_nec)

    assert best_combination is not None
    assert best_combination['layer_1'] == 1
    assert best_combination['num_bars_1'] == 3
    assert best_combination['diameter_1'] ==  16*mm # type: ignore
    assert best_combination['total_as'].value == pytest.approx(0.0006033,rel=1e-3)
    assert best_combination['available_spacing_1'] == 38*mm # type: ignore

def test_beam_longitudinal_ACI_318_19_fail():
    as_nec = 150 * cm**2
    # Expecting a ValueError to be raised
    with pytest.raises(ValueError, match="Cannot fit the required reinforcement within the beam width considering clear cover and spacing."):
        Rebar(beam=section).beam_longitudinal_ACI_318_19(as_req=as_nec)

def test_beam_transverse_ACI_318_19():
    av_nec = 0.00104*m #type: ignore
    best_combination = Rebar(beam=section).beam_transverse_ACI_318_19(A_v_req=av_nec,V_s_req = 24.92*kip, lambda_factor= 1, f_c=4000*psi, d=13.5*inch) #type: ignore

    assert best_combination is not None
    assert best_combination['n_stirrups'] == 1
    assert best_combination['d_b'] == 12*mm # type: ignore
    assert best_combination['s'] ==  6*inch # type: ignore
    assert best_combination['A_v'].value == pytest.approx(0.001484,rel=1e-3)

# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()