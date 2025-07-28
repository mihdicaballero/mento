import pytest
import math
from pint import Quantity

# Import all classes from your material.py and relevant units
from mento.material import (
    Material,
    Concrete,
    Concrete_ACI_318_19,
    Concrete_CIRSOC_201_25,
    Concrete_EN_1992_2004,
    Steel,
    SteelBar,
    SteelStrand,
)
from mento.units import kg, m, MPa, ksi, GPa, psi, Pa, lb, ft, kPa, mm


# Helper to check if a quantity matches another, allowing for slight float differences
def assert_quantity_equal(
    q1: Quantity, q2: Quantity, rtol: float = 1e-6, atol: float = 1e-9
) -> None:
    """Asserts two Pint Quantities are numerically and unit-wise equal."""
    assert q1.units == q2.units, f"Units mismatch: {q1.units} vs {q2.units}"
    assert (
        pytest.approx(q1.magnitude, rel=rtol, abs=atol) == q2.magnitude
    ), f"Magnitudes mismatch: {q1.magnitude} (actual) vs {q2.magnitude} (expected) with units {q1.units}"


# --- Tests for Material (Base Class) ---


def test_material_initialization() -> None:
    """Test the basic Material class initialization."""
    mat = Material(name="Generic Material")
    assert mat.name == "Generic Material"


# --- Tests for Concrete (Base Concrete Class) ---


def test_concrete_initialization_metric_mpa() -> None:
    """Test Concrete initialization with f_c in MPa (metric)."""
    concrete = Concrete(name="H30", f_c=30 * MPa)
    assert concrete.name == "H30"  # Correctly check the name set during initialization
    assert_quantity_equal(concrete.f_c, 30 * MPa)
    assert concrete.design_code == "ACI 318-19"  # Default
    assert concrete.unit_system == "metric"
    assert_quantity_equal(concrete.density, 2500 * kg / m**3)


def test_concrete_initialization_metric_pa() -> None:
    """Test Concrete initialization with f_c in Pa (metric)."""
    concrete = Concrete(name="H25", f_c=25e6 * Pa)
    assert concrete.unit_system == "metric"
    assert_quantity_equal(concrete.density, 2500 * kg / m**3)
    assert_quantity_equal(concrete.f_c, 25e6 * Pa)  # Should convert to MPa implicitly


def test_concrete_initialization_metric_kpa() -> None:
    """Test Concrete initialization with f_c in kPa (metric)."""
    concrete = Concrete(name="H25", f_c=25000 * kPa)
    assert concrete.unit_system == "metric"
    assert_quantity_equal(concrete.density, 2500 * kg / m**3)
    assert_quantity_equal(concrete.f_c, 25000 * kPa)


def test_concrete_initialization_imperial_psi() -> None:
    """Test Concrete initialization with f_c in psi (imperial)."""
    concrete = Concrete(name="C4", f_c=4000 * psi)
    assert_quantity_equal(concrete.f_c, 4000 * psi)
    assert concrete.unit_system == "imperial"
    assert_quantity_equal(concrete.density, 155 * lb / ft**3)


def test_concrete_initialization_imperial_ksi() -> None:
    """Test Concrete initialization with f_c in ksi (imperial)."""
    concrete = Concrete(name="C4", f_c=5 * ksi)
    assert_quantity_equal(concrete.f_c, 5 * ksi)
    assert concrete.unit_system == "imperial"
    assert_quantity_equal(concrete.density, 155 * lb / ft**3)


def test_concrete_initialization_unsupported_unit() -> None:
    """Test Concrete initialization with an unsupported unit for f_c."""
    with pytest.raises(ValueError, match="Unsupported unit system for f_c"):
        Concrete(name="C25", f_c=25 * m)


def test_concrete_get_properties() -> None:
    """Test get_properties method of Concrete base class."""
    concrete = Concrete(name="C30", f_c=30 * MPa)
    props = concrete.get_properties()
    assert_quantity_equal(props["f_c"], 30 * MPa)
    assert_quantity_equal(props["density"], 2500 * kg / m**3)
    assert (
        "design_code" not in props
    )  # Base get_properties only includes f_c and density


# --- Tests for Concrete_ACI_318_19 ---


@pytest.fixture
def aci_metric_concrete() -> Concrete_ACI_318_19:
    """Fixture for a standard metric ACI concrete."""
    return Concrete_ACI_318_19(f_c=25 * MPa, name="C25_ACI")


@pytest.fixture
def aci_imperial_concrete() -> Concrete_ACI_318_19:
    """Fixture for a standard imperial ACI concrete."""
    return Concrete_ACI_318_19(f_c=4000 * psi, name="4ksi_ACI")


def test_aci_initialization_metric(aci_metric_concrete: Concrete_ACI_318_19) -> None:
    """Test ACI concrete initialization in metric."""
    assert aci_metric_concrete.name == "C25_ACI"
    assert_quantity_equal(aci_metric_concrete.f_c, 25 * MPa)
    assert aci_metric_concrete.design_code == "ACI 318-19"
    assert aci_metric_concrete.unit_system == "metric"
    assert_quantity_equal(aci_metric_concrete.density, 2500 * kg / m**3)

    # Expected values for 25 MPa concrete (metric)
    # E_c = (2500^1.5) * 0.043 * sqrt(25) * MPa = 125000 * 0.043 * 5 * MPa = 26875 MPa
    assert_quantity_equal(aci_metric_concrete.E_c, 26875 * MPa)
    # f_r = 0.625 * sqrt(25) * MPa = 0.625 * 5 * MPa = 3.125 MPa
    assert_quantity_equal(aci_metric_concrete.f_r, 3.125 * MPa)
    # beta_1 for 25 MPa (17 <= 25 <= 28) should be 0.85
    assert pytest.approx(aci_metric_concrete.beta_1) == 0.85
    assert pytest.approx(aci_metric_concrete.lambda_factor) == 1.0  # Normalweight
    assert pytest.approx(aci_metric_concrete.phi_v) == 0.75
    assert pytest.approx(aci_metric_concrete.phi_c) == 0.65
    assert (
        pytest.approx(aci_metric_concrete.phi_y) == 0.90
    )  # Note: property is phi_y but uses _phi_t


def test_aci_initialization_imperial(
    aci_imperial_concrete: Concrete_ACI_318_19,
) -> None:
    """Test ACI concrete initialization in imperial."""
    assert aci_imperial_concrete.name == "4ksi_ACI"
    assert_quantity_equal(aci_imperial_concrete.f_c, 4000 * psi)
    assert aci_imperial_concrete.design_code == "ACI 318-19"
    assert aci_imperial_concrete.unit_system == "imperial"
    assert_quantity_equal(aci_imperial_concrete.density, 155 * lb / ft**3)

    # Expected values for 4000 psi concrete (imperial)
    # E_c = (155^1.5) * 33 * sqrt(4000) * psi = 1930.5 * 33 * 63.245 * psi = 4030635.7 psi ~ 4.03e6 psi
    assert_quantity_equal(
        aci_imperial_concrete.E_c, (155**1.5) * 33 * math.sqrt(4000) * psi
    )
    # f_r = 7.5 * sqrt(4000) * psi = 7.5 * 63.245 * psi = 474.34 psi
    assert_quantity_equal(aci_imperial_concrete.f_r, 7.5 * math.sqrt(4000) * psi)
    # For beta_1, f_c must be in MPa. 4000 psi = 27.579 MPa.
    # 17 <= 27.579 <= 28, so beta_1 should be 0.85
    assert pytest.approx(aci_imperial_concrete.beta_1) == 0.85


# Test __beta_1 calculation ranges
@pytest.mark.parametrize(
    "f_c_value, expected_beta_1",
    [
        (17 * MPa, 0.85),
        (28 * MPa, 0.85),
        (20 * MPa, 0.85),  # In the first range
        (30 * MPa, 0.85 - 0.05 / 7 * (30 - 28)),  # 28 < fc <= 55 range
        (40 * MPa, 0.85 - 0.05 / 7 * (40 - 28)),
        (55 * MPa, 0.85 - 0.05 / 7 * (55 - 28)),
        (60 * MPa, 0.65),  # fc > 55 range
        (70 * MPa, 0.65),
        (10 * MPa, 0.85),  # fc < 17 MPa, falls into the else branch (0.85)
    ],
)
def test_aci_beta_1_calculation(f_c_value: Quantity, expected_beta_1: float) -> None:
    concrete = Concrete_ACI_318_19(name="C4", f_c=f_c_value)
    assert pytest.approx(concrete.beta_1) == expected_beta_1


def test_aci_get_properties(aci_metric_concrete: Concrete_ACI_318_19) -> None:
    """Test get_properties for ACI concrete."""
    props = aci_metric_concrete.get_properties()
    assert_quantity_equal(props["f_c"], aci_metric_concrete.f_c)
    assert_quantity_equal(props["density"], aci_metric_concrete.density)
    assert_quantity_equal(props["E_c"], aci_metric_concrete.E_c)
    assert_quantity_equal(props["f_r"], aci_metric_concrete.f_r)
    assert pytest.approx(props["beta_1"]) == aci_metric_concrete.beta_1
    assert pytest.approx(props["epsilon_c"]) == aci_metric_concrete._epsilon_c
    assert pytest.approx(props["lambda"]) == aci_metric_concrete.lambda_factor
    assert pytest.approx(props["phi_v"]) == aci_metric_concrete.phi_v
    assert pytest.approx(props["phi_c"]) == aci_metric_concrete.phi_c
    assert (
        pytest.approx(props["phi_t"]) == aci_metric_concrete.phi_y
    )  # Check the name mapping


def test_aci_str_representation(aci_metric_concrete: Concrete_ACI_318_19) -> None:
    """Test __str__ method for ACI concrete."""
    s = str(aci_metric_concrete)
    assert "Concrete Properties (C25_ACI):" in s
    assert f"  f_c: {aci_metric_concrete.f_c}" in s
    assert f"  Density: {aci_metric_concrete.density}" in s
    assert f"  E_c: {aci_metric_concrete.E_c}" in s
    assert f"  f_r: {aci_metric_concrete.f_r}" in s
    assert f"  beta_1: {aci_metric_concrete.beta_1:.3f}" in s
    assert f"  epsilon_c: {aci_metric_concrete._epsilon_c:.4f}" in s
    assert f"  λ: {aci_metric_concrete.lambda_factor:.2f}" in s
    assert f"  phi_v: {aci_metric_concrete.phi_v:.2f}" in s
    assert f"  phi_c: {aci_metric_concrete.phi_c:.2f}" in s
    assert (
        f"  phi_t: {aci_metric_concrete.phi_y:.2f}" in s
    )  # Using phi_y as per __str__ output


# --- Tests for Concrete_CIRSOC_201_25 ---


def test_cirsoc_initialization() -> None:
    """Test CIRSOC concrete initialization."""
    concrete = Concrete_CIRSOC_201_25(f_c=30 * MPa, name="H30_CIRSOC")
    assert concrete.name == "H30_CIRSOC"
    assert_quantity_equal(concrete.f_c, 30 * MPa)
    assert concrete.design_code == "CIRSOC 201-25"
    assert concrete.unit_system == "metric"
    assert_quantity_equal(concrete.density, 2500 * kg / m**3)
    # Verify inherited properties are calculated based on ACI logic
    assert_quantity_equal(concrete.E_c, ((2500**1.5) * 0.043 * math.sqrt(30)) * MPa)
    assert pytest.approx(concrete.beta_1) == (
        0.85 - 0.05 / 7 * (30 - 28)
    )  # Should follow ACI beta_1 logic


# --- Tests for Concrete_EN_1992_2004 ---


@pytest.fixture
def en_concrete_c30() -> Concrete_EN_1992_2004:
    """Fixture for EN concrete, C30/37 class (fck=30 MPa)."""
    return Concrete_EN_1992_2004(f_c=30 * MPa, name="C30/37")


@pytest.fixture
def en_concrete_c60() -> Concrete_EN_1992_2004:
    """Fixture for EN concrete, C60/75 class (fck=60 MPa)."""
    return Concrete_EN_1992_2004(f_c=60 * MPa, name="C60/75")


def test_en_initialization_c30(en_concrete_c30: Concrete_EN_1992_2004) -> None:
    """Test EN concrete initialization for C30."""
    assert en_concrete_c30.name == "C30/37"
    assert_quantity_equal(en_concrete_c30.f_c, 30 * MPa)  # This is f_ck
    assert en_concrete_c30.design_code == "EN 1992-2004"
    assert en_concrete_c30.unit_system == "metric"
    assert_quantity_equal(en_concrete_c30.density, 2500 * kg / m**3)

    # Expected values for C30/37 (f_ck=30 MPa)
    assert_quantity_equal(en_concrete_c30.f_ck, 30 * MPa)
    assert_quantity_equal(
        en_concrete_c30.f_cm, 30 * MPa + 8 * MPa
    )  # f_ck + 8 MPa = 38 MPa
    # E_cm = 22000 * (38/10)^0.3 * MPa = 22000 * 3.8^0.3 * MPa = 22000 * 1.488 * MPa = 32736 MPa
    assert_quantity_equal(en_concrete_c30.E_cm, 22000 * (3.8) ** 0.3 * MPa)
    # f_ctm = 0.3 * (30)^(2/3) * MPa = 0.3 * 9.6548 * MPa = 2.89644 MPa
    assert_quantity_equal(en_concrete_c30.f_ctm, 0.3 * (30) ** (2 / 3) * MPa)
    assert pytest.approx(en_concrete_c30.epsilon_cu3) == 0.0035
    assert pytest.approx(en_concrete_c30.gamma_c) == 1.5
    assert pytest.approx(en_concrete_c30.gamma_s) == 1.15
    assert pytest.approx(en_concrete_c30.alpha_cc) == 0.85  # As per _alpha_cc_calc
    assert pytest.approx(en_concrete_c30.Lambda_factor) == 0.8  # For f_ck <= 50 MPa
    assert pytest.approx(en_concrete_c30.Eta_factor) == 1.0  # For f_ck <= 50 MPa


def test_en_initialization_c60(en_concrete_c60: Concrete_EN_1992_2004) -> None:
    """Test EN concrete initialization for C60."""
    assert_quantity_equal(en_concrete_c60.f_c, 60 * MPa)  # This is f_ck
    assert en_concrete_c60.design_code == "EN 1992-2004"

    # Expected values for C60/75 (f_ck=60 MPa)
    assert_quantity_equal(en_concrete_c60.f_ck, 60 * MPa)
    assert_quantity_equal(en_concrete_c60.f_cm, 60 * MPa + 8 * MPa)  # 68 MPa
    assert pytest.approx(en_concrete_c60.Lambda_factor) == (
        0.8 - (60 - 50) / 400
    )  # 0.8 - 10/400 = 0.775
    assert pytest.approx(en_concrete_c60.Eta_factor) == (
        1.0 - (60 - 50) / 200
    )  # 1.0 - 10/200 = 0.95


def test_en_get_properties(en_concrete_c30: Concrete_EN_1992_2004) -> None:
    """Test get_properties for EN concrete."""
    props = en_concrete_c30.get_properties()
    assert_quantity_equal(props["f_ck"], en_concrete_c30.f_ck)
    assert_quantity_equal(props["f_cm"], en_concrete_c30.f_cm)
    assert_quantity_equal(props["E_cm"], en_concrete_c30.E_cm)
    assert_quantity_equal(props["f_ctm"], en_concrete_c30.f_ctm)
    assert pytest.approx(props["epsilon_cu3"]) == en_concrete_c30.epsilon_cu3
    assert pytest.approx(props["gamma_c"]) == en_concrete_c30.gamma_c
    assert pytest.approx(props["gamma_s"]) == en_concrete_c30.gamma_s
    assert pytest.approx(props["alpha_cc"]) == en_concrete_c30.alpha_cc
    assert pytest.approx(props["lambda_factor"]) == en_concrete_c30.Lambda_factor
    assert pytest.approx(props["eta_factor"]) == en_concrete_c30.Eta_factor


def test_en_str_representation(en_concrete_c30: Concrete_EN_1992_2004) -> None:
    """Test __str__ method for EN concrete."""
    s = str(en_concrete_c30)
    assert "Concrete Properties (C30/37):" in s
    assert "Design Code: EN 1992-2004" in s
    assert f"  f_c (Characteristic): {en_concrete_c30.f_ck}" in s
    assert f"  f_cm (Mean Compressive): {en_concrete_c30.f_cm}" in s
    assert f"  f_ctm (Mean Tensile): {en_concrete_c30.f_ctm}" in s
    assert f"  E_cm (Secant Modulus): {en_concrete_c30.E_cm}" in s
    assert f"  Density: {en_concrete_c30.density}" in s
    assert f"  ε_cu3: {en_concrete_c30.epsilon_cu3:.4f}" in s
    assert f"  γ_c: {en_concrete_c30.gamma_c:.2f}" in s
    assert f"  γ_s: {en_concrete_c30.gamma_s:.2f}" in s
    assert f"  α_cc: {en_concrete_c30.alpha_cc:.2f}" in s
    assert f"  λ Factor: {en_concrete_c30.Lambda_factor:.2f}" in s
    assert f"  Eta Factor: {en_concrete_c30.Eta_factor:.2f}" in s


# --- Tests for Steel (Base Class) ---


def test_steel_initialization() -> None:
    """Test the Steel base class initialization."""
    steel = Steel(name="Generic Steel", f_y=400 * MPa)
    assert steel.name == "Generic Steel"
    assert_quantity_equal(steel.f_y, 400 * MPa)
    assert_quantity_equal(steel.density, 7850 * kg / m**3)


def test_steel_initialization_custom_density() -> None:
    """Test Steel initialization with custom density."""
    custom_density = 7500 * kg / m**3
    steel = Steel(name="Light Steel", f_y=350 * MPa, density=custom_density)
    assert_quantity_equal(steel.density, custom_density)


# --- Tests for SteelBar ---


@pytest.fixture
def steel_bar_b500s() -> SteelBar:
    """Fixture for a typical SteelBar (B500S)."""
    return SteelBar(name="B500S", f_y=500 * MPa)


@pytest.fixture
def steel_bar_grade60() -> SteelBar:
    """Fixture for a typical SteelBar (Grade 60)."""
    return SteelBar(name="Grade 60", f_y=60 * ksi)


def test_steel_bar_initialization_metric(steel_bar_b500s: SteelBar) -> None:
    """Test SteelBar initialization in metric."""
    assert steel_bar_b500s.name == "B500S"
    assert_quantity_equal(steel_bar_b500s.f_y, 500 * MPa)
    assert_quantity_equal(steel_bar_b500s.E_s, 200 * GPa)
    assert_quantity_equal(steel_bar_b500s.density, 7850 * kg / m**3)
    # epsilon_y = 500 MPa / 200 GPa = 500e6 Pa / 200e9 Pa = 0.0025 (dimensionless)
    assert_quantity_equal(
        steel_bar_b500s.epsilon_y, 0.0025 * mm / mm
    )  # Use a dimensionless unit from pint if available, or just float


def test_steel_bar_initialization_imperial(steel_bar_grade60: SteelBar) -> None:
    """Test SteelBar initialization in imperial."""
    assert steel_bar_grade60.name == "Grade 60"
    assert_quantity_equal(steel_bar_grade60.f_y, 60 * ksi)
    assert_quantity_equal(
        steel_bar_grade60.E_s, 200 * GPa
    )  # E_s is hardcoded to 200 GPa
    # Convert E_s to psi for comparison: 200 GPa = 200e9 Pa = 29007550 psi ~ 29000 ksi
    # epsilon_y = 60 ksi / 29007.55 ksi = 0.002068
    expected_epsilon_y = (60 * ksi) / (200 * GPa).to("ksi")
    assert_quantity_equal(steel_bar_grade60.epsilon_y, expected_epsilon_y)


def test_steel_bar_get_properties(steel_bar_b500s: SteelBar) -> None:
    """Test get_properties for SteelBar."""
    props = steel_bar_b500s.get_properties()
    assert_quantity_equal(props["E_s"], Quantity(steel_bar_b500s.E_s.to("GPa")))
    assert_quantity_equal(props["f_y"], Quantity(steel_bar_b500s.f_y.to("MPa")))
    assert_quantity_equal(props["epsilon_y"], steel_bar_b500s.epsilon_y)


def test_steel_bar_str_representation(steel_bar_b500s: SteelBar) -> None:
    """Test __str__ method for SteelBar."""
    s = str(steel_bar_b500s)
    assert "SteelBar Properties (B500S):" in s
    assert f"  f_y: {steel_bar_b500s.f_y.to('MPa')}" in s
    assert f"  E_s: {steel_bar_b500s.E_s.to('GPa')}" in s
    assert f"  epsilon_y: {steel_bar_b500s.epsilon_y.magnitude:.4f}" in s
    assert f"  Density: {steel_bar_b500s.density}" in s


# --- Tests for SteelStrand ---


@pytest.fixture
def steel_strand_metric() -> SteelStrand:
    """Fixture for a typical SteelStrand (metric)."""
    return SteelStrand(name="Strand_1860", f_y=1600 * MPa)


def test_steel_strand_initialization_metric(steel_strand_metric: SteelStrand) -> None:
    """Test SteelStrand initialization in metric."""
    assert steel_strand_metric.name == "Strand_1860"
    assert_quantity_equal(steel_strand_metric.f_y, 1600 * MPa)
    assert_quantity_equal(steel_strand_metric.f_u, 1860 * MPa)  # Default
    assert_quantity_equal(steel_strand_metric.E_s, 190000 * MPa)  # Default
    assert_quantity_equal(steel_strand_metric.density, 7850 * kg / m**3)
    # epsilon_y = 1600 MPa / 190000 MPa = 0.008421...
    assert_quantity_equal(steel_strand_metric.epsilon_y, (1600 / 190000) * mm / mm)
    assert_quantity_equal(steel_strand_metric.prestress_stress, 0 * MPa)  # Default


def test_steel_strand_get_properties(steel_strand_metric: SteelStrand) -> None:
    """Test get_properties for SteelStrand."""
    props = steel_strand_metric.get_properties()
    assert_quantity_equal(props["E_s"], Quantity(steel_strand_metric.E_s.to("MPa")))
    assert_quantity_equal(props["f_y"], Quantity(steel_strand_metric.f_y.to("MPa")))
    assert_quantity_equal(props["f_u"], Quantity(steel_strand_metric.f_u.to("MPa")))


def test_steel_strand_str_representation(steel_strand_metric: SteelStrand) -> None:
    """Test __str__ method for SteelStrand."""
    s = str(steel_strand_metric)
    assert "SteelStrand Properties (Strand_1860):" in s
    assert f"  f_y: {steel_strand_metric.f_y.to('MPa')}" in s
    assert f"  f_u: {steel_strand_metric.f_u.to('MPa')}" in s
    assert f"  E_s: {steel_strand_metric.E_s.to('MPa')}" in s
    assert f"  Prestress Stress: {steel_strand_metric.prestress_stress}" in s
    assert f"  Density: {steel_strand_metric.density}" in s
