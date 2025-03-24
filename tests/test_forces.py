import pytest
from mento.units import kN, m, kNm
from mento.forces import Forces


def test_default_initialization() -> None:
    # Test default initialization
    f = Forces()
    assert f.N_x == 0 * kN
    assert f.V_z == 0 * kN
    assert f.M_y == 0 * kN * m


def test_custom_initialization() -> None:
    # Test initialization with custom values
    f = Forces(N_x=10 * kN, V_z=5 * kN, M_y=15 * kN * m)
    assert f.N_x == 10 * kN
    assert f.V_z == 5 * kN
    assert f.M_y == 15 * kN * m


def test_setting_attributes() -> None:
    # Test setting attributes
    f = Forces()
    f.set_forces(N_x=20 * kN, V_z=10 * kN, M_y=25 * kNm)
    assert f.N_x == 20 * kN
    assert f.V_z == 10 * kN
    assert f.M_y == 25 * kNm


if __name__ == "__main__":
    pytest.main()
