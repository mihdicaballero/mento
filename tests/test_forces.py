import pytest
from mento.units import kN, m
from mento.forces import Forces

def test_default_initialization() -> None:
    # Test default initialization
    f = Forces()
    assert f.Nx == 0 * kN 
    assert f.Vz == 0 * kN 
    assert f.My == 0 * kN * m 

def test_custom_initialization() -> None:
    # Test initialization with custom values
    f = Forces(Nx=10 * kN, Vz=5 * kN, My=15 * kN * m)  
    assert f.Nx == 10 * kN 
    assert f.Vz == 5 * kN 
    assert f.My == 15 * kN * m 

def test_setting_attributes() -> None:
    # Test setting attributes
    f = Forces()
    f.Nx = 20 * kN 
    f.Vz = 10 * kN 
    f.My = 25 * kN * m 
    assert f.Nx == 20 * kN 
    assert f.Vz == 10 * kN 
    assert f.My == 25 * kN * m 

def test_get_non_existent_attribute() -> None:
    # Test getting a non-existent attribute
    f = Forces()
    with pytest.raises(AttributeError):
        _ = f.non_existent

def test_set_non_existent_attribute() -> None:
    # Test setting a non-existent attribute
    f = Forces()
    with pytest.raises(AttributeError):
        f.non_existent = 10 * kN 

if __name__ == "__main__":
    pytest.main()

