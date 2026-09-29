"""Beam flexure checks that reproduce a case validated outside mento.

Every test here is marked ``published_example`` and says in a ``Source:`` paragraph of
its docstring where its numbers come from: the ETABS/spreadsheet cross-check of the
flexure suite, the Calcpad beam-flexure sheets, The Concrete Centre's guide, the
examples of CRSI's Design Guide on ACI 318 and CSI's ETABS software verification examples.
``tests/architecture/test_published_examples.py`` enforces both.
"""

from tests.helpers import (
    _flexural_reinforcement_in_pint,
    _nominal_moment_double_in_pint,
    _nominal_moment_simple_in_pint,
    _required_flexural_steel_in_pint,
)

import math
import pytest
from mento.node import Node
from mento.beam import RectangularBeam
from mento.material import (
    Concrete_ACI_318_19,
    SteelBar,
    Concrete_EN_1992_2004,
)
from mento.units import psi, kip, inch, ksi, mm, cm, MPa, ft, kNm
from mento.forces import Forces


@pytest.mark.published_example
def test_flexure_EN_1992_2004_matches_concise_eurocode_closed_form() -> None:
    """The required tension steel must match the published EC2 closed form.

    Concise Eurocode 2 (The Concrete Centre), for a singly reinforced
    rectangular section with alpha_cc = 0.85::

        K   = M / (f_ck * b * d**2)
        z/d = [1 + sqrt(1 - 3.529 * K)] / 2
        A_s = M / (0.87 * f_yk * z)

    That closed form is the inversion of ``z = d - 0.4x`` with the EC2
    rectangular block (lambda = 0.8, eta = 1.0). It is computed here from
    scratch, so the test does not depend on any mento formula.

    Source: The Concrete Centre, "How to design concrete structures using Eurocode 2",
    3. Slabs (Moss and Brooker, 2006), p. 3: Figure 1 (procedure for determining
    flexural reinforcement) prints z = d/2*[1 + sqrt(1 - 3.53K)] <= 0.95d and
    A_s = M/(f_yd*z), and Table 5 (z/d for singly reinforced rectangular sections)
    tabulates it: K = 0.09 -> 0.913, K = 0.10 -> 0.902, which brackets the
    z/d = 0.9069 of K = 0.0957 below. The same closed form is tabulated in Concise
    Eurocode 2.
    """
    f_ck, f_yk = 25.0, 500.0
    b, d = 200.0, 560.0  # mm, matches beam_example_EN_1992_2004_01 with 4x16
    M = 150e6  # N*mm

    K = M / (f_ck * b * d**2)
    z_ref = d * (1 + math.sqrt(1 - 3.529 * K)) / 2
    A_s_ref = M / (0.87 * f_yk * z_ref) / 100  # cm2

    beam = RectangularBeam(
        label="concise_ec2",
        concrete=Concrete_EN_1992_2004(name="C25", f_c=f_ck * MPa),
        steel_bar=SteelBar(name="B500S", f_y=f_yk * MPa),
        width=20 * cm,
        height=60 * cm,
        c_c=2.6 * cm,
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=15 * cm)
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    beam.set_longitudinal_rebar_top(n1=0, d_b1=0 * mm)
    assert beam._d_bot.to("mm").magnitude == pytest.approx(d, rel=1e-6)

    results = Node(section=beam, forces=Forces(M_y=150 * kNm)).check_flexure()
    assert results.iloc[1]["As,req bot"] == pytest.approx(A_s_ref, rel=2e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_05() -> None:
    """
    Test_Etabs_05: b=16", h=30", fc=6000psi, fy=60ksi, Mu=200 kip.ft
    Agregado: 2026-05-03
    Singly reinforced. As_min governs over As_calc.
    Excel/ETABS: As_req = 1.7041 in² (= 11.00 cm²)

    Se testea _calculate_flexural_reinforcement_ACI_318_19 directamente
    con los inputs conocidos del Excel (d y d_prima fijos), aislando el
    cálculo de acero requerido de la selección discreta de barras y de la
    iteración del recubrimiento mecánico.

    El caso es representativo del escenario donde el momento aplicado es
    bajo respecto a la sección (As_calc < As_min), por lo que el mínimo
    normativo de ACI 318-19 gobierna el diseño. Se verifica además que
    la sección es simple (sin acero de compresión) y que el flag A_s_bool
    queda apagado: la regla del 4/3 de ACI 9.6.1.3 no alivia nada, porque
    4/3·As_calc no es menor que As_min. A mano (d = 27.5"):
      Rn = 2400/(0.9·16·27.5²) = 0.2204 ksi ; ρ = 0.085·(1 − √(1 − 2·0.2204/5.1)) = 0.003756
      As_calc = 0.003756·16·27.5 = 1.653 in² ; 4/3·As_calc = 2.204 in²
      As_min = max(3√6000, 200)/60000·16·27.5 = 1.7041 in² (§9.6.1.2) → gobierna As_min

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 31 (Test_Etabs_05), column BC (A_s,min governs: column V = column R).
    """

    concrete = Concrete_ACI_318_19(name="fc6000", f_c=6000 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_05",
        concrete=concrete,
        steel_bar=steel,
        width=16 * inch,
        height=30 * inch,
        c_c=1.5 * inch,
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)

    # Inputs directos del Excel (rec mec conocido = 2.5 in)
    d = 27.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft

    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )

    # As_min gobierna: As_final debe ser ≈ 1.7041 in² = 11.00 cm²
    assert A_s_final.to("cm**2").magnitude == pytest.approx(11.00, rel=1e-2)
    # Sección simple: sin acero de compresión
    assert A_s_comp.to("cm**2").magnitude == pytest.approx(0.0, abs=0.01)
    # Flag 4/3 apagado: 4/3·As_calc = 2.204 in² no es menor que As_min, que gobierna
    assert A_s_bool is False


@pytest.mark.published_example
def test_minimum_flexural_reinforcement_ratio_ACI_318_19_Test_Etabs_05() -> None:
    """
    Test_Etabs_05: b=16", h=30", fc=6000psi, fy=60ksi, Mu=200 kip.ft
    Agregado: 2026-05-03

    Se testea _minimum_flexural_reinforcement_ratio_ACI_318_19 directamente.
    ρ_min define el piso de armado — si esta fórmula falla, secciones con
    momento bajo quedan con menos acero del que exige la norma.

    Para fc=6000 psi gobierna 3√fc/fy sobre 200/fy:
    ρ_min = 3√6000/60000 = 0.003873

    ETABS (col BC): As = 1.7041 in², the minimum, which governs this row
    → ρ_min = As_min / (b × d) = 1.7041 / (16 × 27.5) = 0.003873

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 31 (Test_Etabs_05), column BC (A_s = A_s,min = 1.7041 in²).
    """
    from mento.codes.ACI_318_19_beam import _minimum_flexural_reinforcement_ratio_ACI_318_19

    concrete = Concrete_ACI_318_19(name="fc6000", f_c=6000 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_05",
        concrete=concrete,
        steel_bar=steel,
        width=16 * inch,
        height=30 * inch,
        c_c=1.5 * inch,
    )

    Mu = 200 * kip * ft
    rho_min = _minimum_flexural_reinforcement_ratio_ACI_318_19(beam, Mu)

    # Verificación directa contra Excel: As_min = ρ_min × b × d
    d = 27.5 * inch
    As_min = rho_min * beam.width * d
    assert As_min.to("inch**2").magnitude == pytest.approx(1.7041, rel=1e-3)


@pytest.mark.published_example
def test_determine_nominal_moment_simple_reinf_ACI_318_19_Test_Etabs_03() -> None:
    """
    Test_Etabs_03: b=12", h=24", fc=4000psi, fy=60ksi.
    Agregado: 2026-05-03

    Se testea _determine_nominal_moment_simple_reinf_ACI_318_19 directamente.
    La función calcula Mn = As·fy·(d - a/2) con a = As·fy/(0.85·fc·b).
    No aplica φ — eso lo hace la capa superior.

    Con As_req=2.2386 in² (diseñado para Mu=200 kip·ft):
      a = 3.292"  →  Mn = 222.24 kip·ft  →  φMn = 200.02 kip·ft ≈ Mu
    Verifica que la fórmula del bloque de Whitney está bien implementada.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 29 (Test_Etabs_03), columns BC (A_s =
    2.2386 in²) and E (Mu = 200 kip·ft): with ETABS's steel, phi*Mn reproduces Mu.
    """

    concrete = Concrete_ACI_318_19(name="fc4000", f_c=4000 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_03",
        concrete=concrete,
        steel_bar=steel,
        width=12 * inch,
        height=24 * inch,
        c_c=1.5 * inch,
    )

    A_s = 2.2386 * inch**2
    d = 21.5 * inch

    M_n = _nominal_moment_simple_in_pint(beam, A_s, d)

    # Mn ≈ 222.24 kip·ft
    assert M_n.to("kip * ft").magnitude == pytest.approx(222.24, rel=1e-3)
    # φMn ≈ Mu = 200 kip·ft (φ = 0.9)
    phi = 0.9
    assert (phi * M_n).to("kip * ft").magnitude == pytest.approx(200.0, rel=1e-2)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_03() -> None:
    """
    Test_Etabs_03: b=12", h=24", fc=4000psi, fy=60ksi, Mu=200 kip.ft
    Agregado: 2026-05-03

    Sección simple donde As_calc gobierna sobre As_min.
    Excel/ETABS: As_req = 2.2386 in²

    Complemento de Test_Etabs_05: mientras ese caso verifica que el mínimo
    normativo gobierna cuando el momento es bajo, este verifica que cuando
    el momento es suficientemente grande, As_calc (del bloque de Whitney)
    gobierna directamente sin intervención de la regla del 4/3.
    Se confirma además que A_s_bool es False porque la regla del 4/3
    no aplica cuando As_calc > As_min.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 29 (Test_Etabs_03), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc4000", f_c=4000 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_03",
        concrete=concrete,
        steel_bar=steel,
        width=12 * inch,
        height=24 * inch,
        c_c=1.5 * inch,
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)

    # Recubrimiento mecánico 2.5" → d = 24 - 2.5 = 21.5"
    d = 21.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft

    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )

    # As_calc governs: As_final ≈ 2.2386 in²
    assert A_s_final.to("inch**2").magnitude == pytest.approx(2.2386, rel=1e-2)
    # Sección simple: sin acero de compresión
    assert A_s_comp.to("cm**2").magnitude == pytest.approx(0.0, abs=0.01)
    # As_calc > As_min → regla del 4/3 no aplica
    assert A_s_bool is False


@pytest.mark.published_example
def test_determine_nominal_moment_double_reinf_ACI_318_19_Test_Etabs_01() -> None:
    """
    Test_Etabs_01: b=12", h=20", fc=2500psi, fy=60ksi.
    Agregado: 2026-05-03

    Se testea _determine_nominal_moment_double_reinf_ACI_318_19 directamente.
    Excel/ETABS: As=3.0045 in², As_prime=0.7628 in² (diseñado para Mu=200 kip·ft).

    El acero de compresión NO plastifica (ε_s=0.00179 < ε_y=0.00207),
    por lo que la función toma la rama cuadrática para encontrar c.
    Verifica que esa rama está correctamente implementada.

    Resultado esperado: Mn ≈ 222.6 kip·ft → φMn ≈ 200 kip·ft ≈ Mu.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 27 (Test_Etabs_01), columns BC/BE (A_s =
    3.0045, A_s' = 0.7628 in²) and E (Mu = 200 kip·ft): phi*Mn reproduces Mu.
    """

    concrete = Concrete_ACI_318_19(name="fc2500", f_c=2500 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_01",
        concrete=concrete,
        steel_bar=steel,
        width=12 * inch,
        height=20 * inch,
        c_c=1.5 * inch,
    )

    A_s = 3.0045 * inch**2
    A_s_prime = 0.7628 * inch**2
    d = 17.5 * inch
    d_prime = 2.5 * inch

    M_n = _nominal_moment_double_in_pint(beam, A_s, d, d_prime, A_s_prime)

    # Mn ≈ 222.6 kip·ft (rama cuadrática — acero compresión no plastifica)
    assert M_n.to("kip * ft").magnitude == pytest.approx(222.6, rel=1e-2)
    # φMn ≈ Mu = 200 kip·ft
    phi = 0.9
    assert (phi * M_n).to("kip * ft").magnitude == pytest.approx(200.0, rel=1e-2)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_doubly_reinforced_Test_Etabs_01() -> None:
    """
    Test_Etabs_01: b=12", h=20", fc=2500psi, fy=60ksi, Mu=200 kip.ft
    Agregado: 2026-05-03

    Sección doblemente armada. As_calc > As_max → se requiere acero de compresión.
    Excel/ETABS: As_req = 3.0045 in² (19.38 cm²), As_comp = 0.7628 in² (4.92 cm²)

    Complemento de los casos simples (Test_Etabs_03 y Test_Etabs_05):
    verifica que _calculate_flexural_reinforcement_ACI_318_19 detecta
    correctamente el caso doblemente armado y calcula ambos aceros.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 27 (Test_Etabs_01), columns BC and BE.
    """

    concrete = Concrete_ACI_318_19(name="fc2500", f_c=2500 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_01",
        concrete=concrete,
        steel_bar=steel,
        width=12 * inch,
        height=20 * inch,
        c_c=1.5 * inch,
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)

    d = 17.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft

    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )

    # Acero de tracción ≈ 3.0045 in²
    assert A_s_final.to("inch**2").magnitude == pytest.approx(3.0045, rel=1e-2)
    # Acero de compresión ≈ 0.7628 in²
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.7628, rel=1e-2)
    # As_calc > As_min → regla del 4/3 no aplica
    assert A_s_bool is False


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_doubly_reinforced_yielding_Test_Etabs_23() -> None:
    """
    Test_Etabs_23: b=12", h=26", fc=4000psi, fy=60ksi, Mu=500 kip.ft
    Agregado: 2026-05-03
    Sección doblemente armada. Acero de compresión PLASTIFICA (εs' > εy).
    c_t=8.737", εs'=0.00214 > εy=0.00207 → fsprima = fy = 60000 psi.
    Excel/ETABS: As_req = 5.5828 in², As_comp = 0.5647 in²

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 49 (Test_Etabs_23), columns BC and BE.
    """

    concrete = Concrete_ACI_318_19(name="fc4000", f_c=4000 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_23", concrete=concrete, steel_bar=steel, width=12 * inch, height=26 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 23.5 * inch
    d_prima = 2.5 * inch
    Mu = 500 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(5.5828, rel=1e-2)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.5647, rel=1e-2)
    assert A_s_bool is False


@pytest.mark.published_example
def test_determine_nominal_moment_double_reinf_ACI_318_19_Test_Etabs_23_yielding() -> None:
    """
    Test_Etabs_23: b=12", h=26", fc=4000psi, fy=60ksi.
    Agregado: 2026-05-03
    Acero compresión PLASTIFICA → rama de plastificación.
    εs' = 0.00214 > εy = 0.00207 → fsprima = fy.
    φMn ≈ Mu = 500 kip·ft (ETABS validado).

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 49 (Test_Etabs_23), columns BC/BE (A_s =
    5.5828, A_s' = 0.5647 in²) and E (Mu = 500 kip·ft): phi*Mn reproduces Mu.
    """

    concrete = Concrete_ACI_318_19(name="fc4000", f_c=4000 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_23", concrete=concrete, steel_bar=steel, width=12 * inch, height=26 * inch, c_c=1.5 * inch
    )
    A_s = 5.582781 * inch**2
    A_s_prime = 0.56469 * inch**2
    d = 23.5 * inch
    d_prime = 2.5 * inch
    M_n = _nominal_moment_double_in_pint(beam, A_s, d, d_prime, A_s_prime)
    phi = 0.9
    assert (phi * M_n).to("kip * ft").magnitude == pytest.approx(500.0, rel=1e-2)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_04() -> None:
    """
    Test_Etabs_04: b=16", h=30", fc=5000psi, fy=60ksi, Mu=200 kip.ft
    Agregado: 2026-05-03
    Sección simple. Testea β₁=0.80 (fc=5000 psi).
    As_calc gobierna sobre As_min (1.6604 > 1.5556 in²).
    Excel/ETABS: As_req = 1.6604 in², As_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 30 (Test_Etabs_04), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc5000", f_c=5000 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_04", concrete=concrete, steel_bar=steel, width=16 * inch, height=30 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 27.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(1.6604, rel=1e-2)
    assert A_s_comp.to("cm**2").magnitude == pytest.approx(0.0, abs=0.01)
    assert A_s_bool is False


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_02() -> None:
    """
    Test_Etabs_02: b=12", h=18", fc=3000psi, fy=60ksi, Mu=200 kip.ft
    Doubly reinforced (A_s_comp > 0).
    ETABS validated: A_s_final = 3.4090 in², A_s_comp = 1.1701 in².
    Note: the spreadsheet uses d_prima = 2.5" (d_prima_conocido column), not
    the 0.1*h default of 1.8".

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 28 (Test_Etabs_02), columns BC and BE.
    """

    concrete = Concrete_ACI_318_19(name="fc3000", f_c=3000 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_02", concrete=concrete, steel_bar=steel, width=12 * inch, height=18 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 15.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(3.4090, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(1.1701, rel=1e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_06() -> None:
    """
    Test_Etabs_06: b=16", h=30", fc=7000psi, fy=60ksi, Mu=200 kip.ft
    Simple section, ETABS validated: A_s_final = 1.8407 in², A_s_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 32 (Test_Etabs_06), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc7000", f_c=7000 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_06", concrete=concrete, steel_bar=steel, width=16 * inch, height=30 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 27.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(1.8407, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.0, abs=1e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_11() -> None:
    """
    Test_Etabs_11: b=24", h=20", fc=12000psi, fy=60ksi, Mu=200 kip.ft
    Simple section, ETABS validated: A_s_final = 2.5865 in², A_s_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 37 (Test_Etabs_11), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc12000", f_c=12000 * psi)
    steel = SteelBar(name="fy60", f_y=60 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_11", concrete=concrete, steel_bar=steel, width=24 * inch, height=20 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 17.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(2.5865, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.0, abs=1e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_12() -> None:
    """
    Test_Etabs_12: b=10", h=24", fc=2500psi, fy=75ksi, Mu=200 kip.ft
    Doubly reinforced (A_s_comp > 0).
    ETABS validated: A_s_final = 1.9373 in², A_s_comp = 0.1719 in².
    Note: the spreadsheet uses d_prima = 2.5" (d_prima_conocido column), not
    the 0.1*h default of 2.4".

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 38 (Test_Etabs_12), columns BC and BE.
    """

    concrete = Concrete_ACI_318_19(name="fc2500", f_c=2500 * psi)
    steel = SteelBar(name="fy75", f_y=75 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_12", concrete=concrete, steel_bar=steel, width=10 * inch, height=24 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 21.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(1.9373, rel=1e-3)
    # A_s_comp uses rel=2e-3 (0.2%) because with very low f_c (2500 psi) the
    # displaced-concrete correction (0.85 * f_c) is small, and any 4-decimal
    # rounding in the ETABS reference value amplifies to ~0.15% relative error
    # on the small (0.17 in²) A_s_comp result.
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.1719, rel=2e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_13() -> None:
    """
    Test_Etabs_13: b=10", h=24", fc=3000psi, fy=75ksi, Mu=200 kip.ft
    Simple section, ETABS validated: A_s_final = 1.9009 in², A_s_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 39 (Test_Etabs_13), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc3000", f_c=3000 * psi)
    steel = SteelBar(name="fy75", f_y=75 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_13", concrete=concrete, steel_bar=steel, width=10 * inch, height=24 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 21.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(1.9009, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.0, abs=1e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_14() -> None:
    """
    Test_Etabs_14: b=10", h=24", fc=4000psi, fy=75ksi, Mu=200 kip.ft
    Simple section, ETABS validated: A_s_final = 1.8245 in², A_s_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 40 (Test_Etabs_14), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc4000", f_c=4000 * psi)
    steel = SteelBar(name="fy75", f_y=75 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_14", concrete=concrete, steel_bar=steel, width=10 * inch, height=24 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 21.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(1.8245, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.0, abs=1e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_15() -> None:
    """
    Test_Etabs_15: b=14", h=18", fc=5000psi, fy=75ksi, Mu=200 kip.ft
    Simple section, ETABS validated: A_s_final = 2.5605 in², A_s_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 41 (Test_Etabs_15), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc5000", f_c=5000 * psi)
    steel = SteelBar(name="fy75", f_y=75 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_15", concrete=concrete, steel_bar=steel, width=14 * inch, height=18 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 15.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(2.5605, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.0, abs=1e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_16() -> None:
    """
    Test_Etabs_16: b=14", h=20", fc=6000psi, fy=75ksi, Mu=200 kip.ft
    Simple section, ETABS validated: A_s_final = 2.1735 in², A_s_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 42 (Test_Etabs_16), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc6000", f_c=6000 * psi)
    steel = SteelBar(name="fy75", f_y=75 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_16", concrete=concrete, steel_bar=steel, width=14 * inch, height=20 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 17.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(2.1735, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.0, abs=1e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_17() -> None:
    """
    Test_Etabs_17: b=14", h=19", fc=7000psi, fy=75ksi, Mu=200 kip.ft
    Simple section, ETABS validated: A_s_final = 2.2991 in², A_s_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 43 (Test_Etabs_17), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc7000", f_c=7000 * psi)
    steel = SteelBar(name="fy75", f_y=75 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_17", concrete=concrete, steel_bar=steel, width=14 * inch, height=19 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 16.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(2.2991, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.0, abs=1e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_19() -> None:
    """
    Test_Etabs_19: b=16", h=15", fc=9000psi, fy=75ksi, Mu=200 kip.ft
    Simple section, ETABS validated: A_s_final = 3.0764 in², A_s_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 45 (Test_Etabs_19), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc9000", f_c=9000 * psi)
    steel = SteelBar(name="fy75", f_y=75 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_19", concrete=concrete, steel_bar=steel, width=16 * inch, height=15 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 12.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(3.0764, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.0, abs=1e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_20() -> None:
    """
    Test_Etabs_20: b=20", h=12", fc=10000psi, fy=75ksi, Mu=200 kip.ft
    Simple section (shallow, wide). ETABS validated: A_s_final = 4.1408 in², A_s_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 46 (Test_Etabs_20), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc10000", f_c=10000 * psi)
    steel = SteelBar(name="fy75", f_y=75 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_20", concrete=concrete, steel_bar=steel, width=20 * inch, height=12 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 9.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(4.1408, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.0, abs=1e-3)


@pytest.mark.published_example
def test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_22() -> None:
    """
    Test_Etabs_22: b=26", h=12", fc=12000psi, fy=75ksi, Mu=200 kip.ft
    Simple section (shallow, very wide). ETABS validated: A_s_final = 3.9783 in², A_s_comp = 0.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 48 (Test_Etabs_22), column BC.
    """

    concrete = Concrete_ACI_318_19(name="fc12000", f_c=12000 * psi)
    steel = SteelBar(name="fy75", f_y=75 * ksi)
    beam = RectangularBeam(
        label="Test_Etabs_22", concrete=concrete, steel_bar=steel, width=26 * inch, height=12 * inch, c_c=1.5 * inch
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.375 * inch, s_l=12 * inch)
    d = 9.5 * inch
    d_prima = 2.5 * inch
    Mu = 200 * kip * ft
    A_s_min, A_s_max, A_s_final, A_s_comp, c_d, A_s_bool, _doubly = _flexural_reinforcement_in_pint(
        beam, Mu, d, d_prima
    )
    assert A_s_final.to("inch**2").magnitude == pytest.approx(3.9783, rel=1e-3)
    assert A_s_comp.to("inch**2").magnitude == pytest.approx(0.0, abs=1e-3)


# ---------------------------------------------------------------------------
# CRSI, "Design Guide on the ACI 318 Building Code Requirements for Structural
# Concrete" (ACI 318-19 edition), Chapter 6, §6.9 Examples. The book is kept
# outside the repository ("CRSI ACI318-19 Beam design examples.pdf").
#
# All the examples use f'c = 4000 psi and Grade 60 bars, and size the steel with
# a fixed d (h - 2.5 in. for one layer, h - 3.5 in. for two), which is why they
# go through the reinforcement helper with that d rather than through a design
# that picks bars. A negative section is the web b_w; a positive section is a
# T-beam whose stress block stays in the flange (a < h_f), which the book solves
# as a rectangle of width b_f. There only A_s_calc is comparable: the book keeps
# the minimum of §9.6.1.2 on b_w, and a rectangle of width b_f puts it on b_f.
#
# The book prints A_s to 0.01 in.² after rounding R_n to whole psi, so the
# areas are held to one unit of that last digit. Its A_s,t = 0.018*b*d is the
# tension-controlled limit written for eps_t = 0.005; mento reads Table 21.2.2
# exactly (eps_ty + 0.003, rho = 0.0179), 0.5 % less, hence rel=1e-2 there.
# ---------------------------------------------------------------------------

_CRSI_AREA_TOL = 0.01  # in.², one unit of the last digit the book prints


def _crsi_beam(label: str, b: float, h: float) -> RectangularBeam:
    return RectangularBeam(
        label=label,
        concrete=Concrete_ACI_318_19(name="fc 4000", f_c=4000 * psi),
        steel_bar=SteelBar(name="Grade 60", f_y=60 * ksi),
        width=b * inch,
        height=h * inch,
        c_c=1.5 * inch,
    )


def _assert_crsi_flexure(
    label: str,
    b: float,
    h: float,
    d: float,
    M_u: float,
    A_s: float | None,
    A_s_calc: float | None,
    A_s_min: float | None,
    A_s_max: float | None,
) -> None:
    """Size one CRSI section and hold each area the book prints; None is not printed."""
    beam = _crsi_beam(label, b, h)
    got_min, got_max, got_final, got_calc = _required_flexural_steel_in_pint(beam, M_u * kip * ft, d * inch, 2.5 * inch)
    if A_s is not None:
        assert got_final.to("inch**2").magnitude == pytest.approx(A_s, abs=_CRSI_AREA_TOL)
    if A_s_calc is not None:
        assert got_calc.to("inch**2").magnitude == pytest.approx(A_s_calc, abs=_CRSI_AREA_TOL)
    if A_s_min is not None:
        assert got_min.to("inch**2").magnitude == pytest.approx(A_s_min, abs=_CRSI_AREA_TOL)
    if A_s_max is not None:
        assert got_max.to("inch**2").magnitude == pytest.approx(A_s_max, rel=1e-2)


@pytest.mark.published_example
@pytest.mark.parametrize(
    ("location", "M_u", "A_s"),
    [
        # Only the minimum: A_s_calc = 0.51 in.², and mento takes the 4/3 relief
        # of §9.6.1.3 (floored at its 1.8 ‰ of b*h, 1.21 in.²), where the book
        # lays the full A_s,min. Both satisfy the code; the minimum is common.
        ("exterior_negative", 49.2, None),
        ("first_interior_negative", 165.0, 2.01),
        ("interior_negative", 156.1, 2.01),
    ],
)
def test_flexure_ACI_318_19_CRSI_example_6_2(location: str, M_u: float, A_s: float | None) -> None:
    """Beam 28 x 24 in., d = 21.5 in.: every negative section is governed by A_s,min.

    A_s,min = 200*b_w*d/f_y = 2.01 in.² (§9.6.1.2) and A_s,t = 0.018*b*d = 10.8 in.².
    For 165.0 and 156.1 kip·ft, 4/3 of A_s_calc (1.75 and 1.65 in.²) exceeds the
    minimum, so the relief of §9.6.1.3 does not apply and A_s = A_s,min. The
    positive sections are left out: the book prints only their minimum on b_w.

    Source: CRSI, Design Guide on the ACI 318 Building Code Requirements for
    Structural Concrete, §6.9.2, Example 6.2, p. 6-66, Table 6.24.
    """
    _assert_crsi_flexure(f"CRSI-6.2-{location}", 28, 24, 21.5, M_u, A_s, None, 2.01, 10.8)


@pytest.mark.published_example
@pytest.mark.parametrize(
    ("location", "b", "d", "M_u", "A_s", "A_s_calc", "A_s_min", "A_s_max"),
    [
        ("exterior_negative", 12, 21.5, 84.9, 0.91, 0.91, 0.86, 4.64),
        ("first_interior_negative", 12, 21.5, 258.3, 2.97, 2.97, 0.86, 4.64),
        # T-beam, b_f = 76.5 in., two layers (d = 20.5 in.), a = 0.54 in. < h_f.
        ("positive", 76.5, 20.5, 214.7, None, 2.36, None, None),
    ],
)
def test_flexure_ACI_318_19_CRSI_example_6_6(
    location: str,
    b: float,
    d: float,
    M_u: float,
    A_s: float | None,
    A_s_calc: float,
    A_s_min: float | None,
    A_s_max: float | None,
) -> None:
    """Beam 12 x 24 in. of the LFRS, end span: strength governs every section.

    Source: CRSI, Design Guide on the ACI 318 Building Code Requirements for
    Structural Concrete, §6.9.6, Example 6.6, p. 6-80 (b_f, d and a of the positive
    section) and p. 6-81, Table 6.28.
    """
    _assert_crsi_flexure(f"CRSI-6.6-{location}", b, 24, d, M_u, A_s, A_s_calc, A_s_min, A_s_max)


@pytest.mark.published_example
@pytest.mark.parametrize(
    ("location", "b", "M_u", "A_s", "A_s_calc", "A_s_min", "A_s_max"),
    [
        ("end_exterior_negative", 7, 100.9, 0.90, 0.90, 0.61, 3.28),
        ("end_positive", 60, 172.9, None, 1.49, None, None),
        ("end_first_interior_negative", 7, 201.3, 1.89, 1.89, 0.61, 3.28),
        ("interior_positive", 60, 102.7, None, 0.88, None, None),
        ("interior_negative", 7, 183.0, 1.71, 1.71, 0.61, 3.28),
    ],
)
def test_flexure_ACI_318_19_CRSI_example_6_11(
    location: str,
    b: float,
    M_u: float,
    A_s: float | None,
    A_s_calc: float,
    A_s_min: float | None,
    A_s_max: float | None,
) -> None:
    """Joist 7 x 28.5 in., d = 26.0 in.; positive sections as a T of b_f = 60 in.

    Source: CRSI, Design Guide on the ACI 318 Building Code Requirements for
    Structural Concrete, §6.9.11, Example 6.11, p. 6-98 (b_f) and p. 6-99, Table 6.31.
    """
    _assert_crsi_flexure(f"CRSI-6.11-{location}", b, 28.5, 26.0, M_u, A_s, A_s_calc, A_s_min, A_s_max)


@pytest.mark.published_example
@pytest.mark.parametrize(
    ("location", "b", "M_u", "A_s", "A_s_calc", "A_s_min", "A_s_max"),
    [
        ("end_exterior_negative", 30, 368.7, 3.27, 3.27, 2.60, 14.0),
        ("end_positive", 57, 421.3, None, 3.68, None, None),
        ("end_first_interior_negative", 30, 585.6, 5.33, 5.33, 2.60, 14.0),
        ("interior_positive", 57, 363.3, None, 3.17, None, None),
        ("interior_negative", 30, 532.4, 4.81, 4.81, 2.60, 14.0),
    ],
)
def test_flexure_ACI_318_19_CRSI_example_6_16(
    location: str,
    b: float,
    M_u: float,
    A_s: float | None,
    A_s_calc: float,
    A_s_min: float | None,
    A_s_max: float | None,
) -> None:
    """Edge beam 30 x 28.5 in., d = 26.0 in.; positive sections as an L of b_f = 57 in.

    Source: CRSI, Design Guide on the ACI 318 Building Code Requirements for
    Structural Concrete, §6.9.16, Example 6.16, p. 6-112 (b_f) and p. 6-113, Table 6.34.
    """
    _assert_crsi_flexure(f"CRSI-6.16-{location}", b, 28.5, 26.0, M_u, A_s, A_s_calc, A_s_min, A_s_max)


@pytest.mark.published_example
@pytest.mark.parametrize(
    ("location", "b", "M_u", "A_s", "A_s_calc", "A_s_min", "A_s_max"),
    [
        # Only the minimum: A_s_calc = 0.07 in.², and mento takes the 4/3 relief
        # of §9.6.1.3 (floored at its 1.8 ‰ of b*h, 0.36 in.²), where the book
        # lays the full A_s,min. Both satisfy the code; the minimum is common.
        ("exterior_negative", 7, 8.5, None, None, 0.61, 3.28),
        ("positive", 60, 204.5, None, 1.77, None, None),
        ("first_interior_negative", 7, 245.5, 2.37, 2.37, 0.61, 3.28),
    ],
)
def test_flexure_ACI_318_19_CRSI_example_6_18(
    location: str,
    b: float,
    M_u: float,
    A_s: float | None,
    A_s_calc: float | None,
    A_s_min: float | None,
    A_s_max: float | None,
) -> None:
    """The joist of Example 6.11 with its end moments redistributed for torsion (§22.7.3.3).

    Source: CRSI, Design Guide on the ACI 318 Building Code Requirements for
    Structural Concrete, §6.9.18, Example 6.18, p. 6-118, Table 6.35.
    """
    _assert_crsi_flexure(f"CRSI-6.18-{location}", b, 28.5, 26.0, M_u, A_s, A_s_calc, A_s_min, A_s_max)


# ---------------------------------------------------------------------------
# CSI Software Verification, "ACI 318-19 Example 001" and "EN 2-2004 Example
# 001" (ETABS): a simply supported singly reinforced rectangle designed for
# flexure and shear, checked by hand in the same documents. The PDFs are kept
# outside the repository.
# ---------------------------------------------------------------------------


@pytest.mark.published_example
def test_flexure_ACI_318_19_ETABS_example_001() -> None:
    """Beam 10 x 16 in., d = 13.5 in., f'c = 4 ksi, f_y = 60 ksi, M_u = 1460.4 kip·in.

    A_s,min = 200*b_w*d/f_y = 0.450 in.² governs over 3*sqrt(f'c)*b_w*d/f_y = 0.427 in.²;
    a = 4.183 in. < a_max = 4.266 in., so the section is singly reinforced and
    A_s = M_u/(phi*f_y*(d - a/2)) = 2.37 in.².

    Source: CSI Software Verification, ETABS, "ACI 318-19 Example 001", p. 2 (Results
    Comparison) and p. 3 (Hand Calculation, Flexural Design).
    """
    beam = RectangularBeam(
        label="ETABS-ACI-Ex001",
        concrete=Concrete_ACI_318_19(name="fc 4000", f_c=4000 * psi),
        steel_bar=SteelBar(name="Grade 60", f_y=60 * ksi),
        width=10 * inch,
        height=16 * inch,
        c_c=1.5 * inch,
    )
    A_s_min, _A_s_max, A_s, A_s_calc = _required_flexural_steel_in_pint(
        beam, 1460.4 * kip * inch, 13.5 * inch, 2.5 * inch
    )
    assert A_s_min.to("inch**2").magnitude == pytest.approx(0.450, abs=5e-4)
    assert A_s_calc.to("inch**2").magnitude == pytest.approx(2.37, abs=5e-3)
    assert A_s.to("inch**2").magnitude == pytest.approx(2.37, abs=5e-3)


@pytest.mark.published_example
def test_flexure_EN_1992_2004_ETABS_example_001() -> None:
    """Beam 230 x 550 mm, d = 490 mm, C30, f_yk = 460 MPa, M_Ed = 165.015 kN·m.

    mento takes alpha_cc = 0.85, gamma_c = 1.5 and gamma_s = 1.15, the factors of the
    national annexes that give A_s = 933 mm² in Table 2 (Finland, Germany, Ireland,
    Norway, Portugal, Singapore, UK): f_cd = 17 MPa, m = 0.1758, omega = 0.1947. The
    minimum is max(0.26*f_ctm/f_yk, 0.0013)*b*d = 184.5 mm². Stirrup Ø10 and bars Ø20
    under a cover of 40 mm put d at the example's 490 mm.

    Source: CSI Software Verification, ETABS/SAFE, "EN 2-2004 Example 001", p. 3 (Table 2),
    p. 6 (A_s,min) and p. 19 (hand calculation for the United Kingdom).
    """
    beam = RectangularBeam(
        label="ETABS-EN-Ex001",
        concrete=Concrete_EN_1992_2004(name="C30", f_c=30 * MPa),
        steel_bar=SteelBar(name="fyk 460", f_y=460 * MPa),
        width=230 * mm,
        height=550 * mm,
        c_c=40 * mm,
    )
    beam.set_transverse_rebar(n_stirrups=1, d_b=10 * mm, s_l=20 * cm)
    beam.set_longitudinal_rebar_bot(n1=3, d_b1=20 * mm)
    beam.set_longitudinal_rebar_top(n1=2, d_b1=20 * mm)
    assert beam._d_bot.to("mm").magnitude == pytest.approx(490, rel=1e-9)
    Node(section=beam, forces=Forces(label="Combo1", M_y=165.015 * kNm)).check_flexure()
    bottom = beam.flexure_checks[0].bottom
    assert bottom.A_s_req.to("mm**2").magnitude == pytest.approx(933, abs=0.5)
    assert bottom.A_s_min.to("mm**2").magnitude == pytest.approx(184.5, abs=0.05)
