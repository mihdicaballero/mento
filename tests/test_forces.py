import pytest
import forallpeople
forallpeople.environment('structural', top_level=True)

# Import the Forces class
from structurelab.forces import Forces

def test_default_initialization():
    # Test default initialization
    f = Forces()
    assert f.Nx == 0 * kN #type: ignore
    assert f.Vz == 0 * kN #type: ignore
    assert f.My == 0 * kN * m #type: ignore

def test_custom_initialization():
    # Test initialization with custom values
    f = Forces(Nx=10 * kN, Vz=5 * kN, My=15 * kN * m)  #type: ignore
    assert f.Nx == 10 * kN #type: ignore
    assert f.Vz == 5 * kN #type: ignore
    assert f.My == 15 * kN * m #type: ignore

def test_setting_attributes():
    # Test setting attributes
    f = Forces()
    f.Nx = 20 * kN #type: ignore
    f.Vz = 10 * kN #type: ignore
    f.My = 25 * kN * m #type: ignore
    assert f.Nx == 20 * kN #type: ignore
    assert f.Vz == 10 * kN #type: ignore
    assert f.My == 25 * kN * m #type: ignore

def test_get_non_existent_attribute():
    # Test getting a non-existent attribute
    f = Forces()
    with pytest.raises(AttributeError):
        _ = f.non_existent

def test_set_non_existent_attribute():
    # Test setting a non-existent attribute
    f = Forces()
    with pytest.raises(AttributeError):
        f.non_existent = 10 * kN #type: ignore

def test_additional_kwargs():
    # Test initialization with additional keyword arguments
    f = Forces(extra_force=10 * kN) #type: ignore
    assert f.extra_force == 10 * kN #type: ignore

    # Test setting an additional keyword argument
    f.extra_force = 20 * kN #type: ignore
    assert f.extra_force == 20 * kN #type: ignore

if __name__ == "__main__":
    pytest.main()

