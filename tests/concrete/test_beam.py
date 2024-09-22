import sys
import os
from structurelab.concrete.beam import Beam
from structurelab import material

# Add the project root to PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
import forallpeople
forallpeople.environment('structural', top_level=True)


@pytest.fixture()
def beam_example():
    # Define custom settings
    custom_settings = {
        'clear_cover': 1.5*inch,  # type: ignore
        'stirrup_diameter': 0.5*inch,  # type: ignore
        'longitudinal_diameter': 1*inch  # type: ignore
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
    return section

def test_shear_check_ACI_318_19(beam_example):
    V_u = 37.727*kip  # type: ignore
    N_u = 0*kip  # type: ignore
    A_s = 0.847*inch**2  # type: ignore
    results = beam_example.check_shear_ACI_318_19(V_u, N_u, A_s, d_b=0.5*inch, s=6*inch, n_legs=2)  # type: ignore

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results['A_v_min'].value  == pytest.approx(0.0002116, rel=1e-3)
    assert results['A_v_req'].value  == pytest.approx(0.001042, rel=1e-3)
    assert results['A_v'].value  == pytest.approx(0.001662, rel=1e-3)
    assert results['phi_V_c'].value  == pytest.approx(56969.4, rel=1e-3)
    assert results['phi_V_s'].value  == pytest.approx(176865, rel=1e-3)
    assert results['phi_V_n'].value  == pytest.approx(233834, rel=1e-3)
    assert results['phi_V_max'].value  == pytest.approx(284847, rel=1e-3)
    assert results["FUv"]  == pytest.approx(0.7176789, rel=1e-3)

    # Assert non-numeric values directly
    assert results['shear_ok'] is True
    assert results['max_shear_ok'] is True

def test_shear_design_ACI_318_19(beam_example):
    V_u = 37.727*kip  # type: ignore
    N_u = 0*kip  # type: ignore
    A_s = 0.847*inch**2  # type: ignore
    results = beam_example.design_shear_ACI_318_19(V_u, N_u, A_s)  # type: ignore

    # Compare dictionaries with a tolerance for floating-point values, in m 
    assert results['A_v_min'].value  == pytest.approx(0.0002116, rel=1e-3)
    assert results['A_v_req'].value  == pytest.approx(0.001042, rel=1e-3)
    assert results['A_v'].value  == pytest.approx(0.0014842, rel=1e-3)
    assert results['phi_V_c'].value  == pytest.approx(56969.4, rel=1e-3)
    assert results['phi_V_s'].value  == pytest.approx(157905, rel=1e-3)
    assert results['phi_V_n'].value  == pytest.approx(214875, rel=1e-3)
    assert results['phi_V_max'].value  == pytest.approx(284847, rel=1e-3)
    assert results["FUv"]  == pytest.approx(0.7810047891, rel=1e-3)

    # Assert non-numeric values directly
    assert results['shear_ok'] is True
    assert results['max_shear_ok'] is True


# This is where pytest will collect the tests and run them
if __name__ == "__main__":
    pytest.main()