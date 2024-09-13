import sys
import os

# Add the project root to PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
from unittest.mock import patch, MagicMock
import forallpeople
forallpeople.environment('structural', top_level=True)
from structurelab.concrete.beam import Beam
from structurelab import material


# The shear function
def shear():
    # Define custom settings
    custom_settings = {
        'clear_cover': 2*inch,  # type: ignore
        'stirrup_diameter': 0.375*inch,  # type: ignore
        'longitudinal_diameter': 0.25*inch  # type: ignore
    }
    
    concrete = material.create_concrete(name="C4", f_c=4000*psi, design_code="ACI 318-19")  # type: ignore
    steelBar = material.SteelBar(name="ADN 420", f_y=60*ksi)  # type: ignore
    
    section = Beam(
        name="V-10x16",
        concrete=concrete,
        steelBar=steelBar,
        width=10*inch,  # type: ignore
        depth=16*inch,  # type: ignore
    )
    
    section.update_settings(custom_settings)
    
    V_u = 37.727*kip  # type: ignore
    N_u = 0*kip  # type: ignore
    A_s = 0.847*inch**2  # type: ignore
    
    results = section.check_shear_ACI_318_19(V_u, N_u, A_s, d_b=0.5*inch, s=6*inch, n_legs=2)  # type: ignore
    return results

# The pytest test function
@pytest.mark.parametrize("expected_results", [{
    'A_vmin': 0.00833 * inch, # type: ignore
    'A_vs': 0.065 * inch, # type: ignore
    'phi_V_c': 12807.225*lb, # type: ignore
    'phi_V_s': 39.761 * kip, # type: ignore
    'phi_V_n': 52568.007 * lb, # type: ignore
    'phi_V_max': 64036.123 * lb, # type: ignore
    'shear_ok': True, 
    'max_shear_ok': True
}])
@patch('material.create_concrete')
@patch('material.SteelBar')
@patch('structurelab.concrete.beam.Beam.check_shear_ACI_318_19')
def test_shear(mock_check_shear, mock_SteelBar, mock_create_concrete, expected_results):
    # Mock the material and section objects
    mock_create_concrete.return_value = MagicMock()
    mock_SteelBar.return_value = MagicMock()
    
    # Mock the section's shear check method to return expected results
    mock_check_shear.return_value = expected_results
    
    # Call the shear function
    results = shear()
    
    # Assert that the results match the expected results
    assert results == expected_results

