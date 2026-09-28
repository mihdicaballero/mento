"""Beam flexure checks that reproduce a case validated outside mento.

Every test here is marked ``published_example`` and says in a ``Source:`` paragraph of
its docstring where its numbers come from: the ETABS/spreadsheet cross-check of the
flexure suite, the Calcpad beam-flexure sheets and The Concrete Centre's guide.
``tests/architecture/test_published_examples.py`` enforces both.
"""

from tests.helpers import (
    _flexural_reinforcement_in_pint,
    _nominal_moment_double_in_pint,
    _nominal_moment_simple_in_pint,
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
def test_flexure_check_EN_1992_2004_01(
    beam_example_EN_1992_2004_01: RectangularBeam,
) -> None:
    """Singly reinforced check: 20x60, C25/B500S, 4Ø16, M_Ed = 150 kN·m.

    Source: Calcpad "EN 1992-1-1_2004 Beam Flexure 01 - Metric v2.cpd" for the case and
    A_s,min; A_s,req and M_Rd are re-derived by hand below with the closed form of
    The Concrete Centre's "How to design concrete structures using Eurocode 2",
    3. Slabs (Moss and Brooker, 2006), p. 3, Figure 1 and Table 5 (z/d for singly
    reinforced rectangular sections) -- the v2 sheet applied lambda twice to the
    lever arm and is not the source of those two.
    """
    f = Forces(M_y=150 * kNm)
    beam_example_EN_1992_2004_01.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=15 * cm)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    beam_example_EN_1992_2004_01.set_longitudinal_rebar_top(n1=0, d_b1=0 * mm)
    # Read back before any check has run: this asserts the layout that was just
    # set, not a result, which is what `reinforcement` is for.
    rebar = beam_example_EN_1992_2004_01.reinforcement
    assert rebar.bottom.A_s.to(cm**2).magnitude == pytest.approx(8.042, rel=1e-2)
    assert beam_example_EN_1992_2004_01._d_bot.to(cm).magnitude == pytest.approx(56.0, rel=1e-2)
    assert beam_example_EN_1992_2004_01.width.to(cm).magnitude == pytest.approx(20.0, rel=1e-2)
    assert rebar.transverse.d_b.to(mm).magnitude == pytest.approx(6.0, rel=1e-2)
    node = Node(section=beam_example_EN_1992_2004_01, forces=f)
    results = node.check_flexure()
    assert results.iloc[1]["Label"] == "B_Example_EN_01"
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["As,min"] == pytest.approx(1.49, rel=1e-2)
    # Cross-checked against the Concise Eurocode 2 closed form (The Concrete
    # Centre): K = M/(f_ck*b*d^2) = 0.09566, z/d = [1+sqrt(1-3.529K)]/2 -> z =
    # 507.89 mm, A_s = M/(0.87*f_yk*z) = 6.79 cm^2. Pure ES=0/EM=0 equilibrium
    # with the EC2 block gives the same. The previous 6.656 came from applying
    # lambda twice to the lever arm.
    assert results.iloc[1]["As,req bot"] == pytest.approx(6.79, rel=1e-3)
    assert results.iloc[1]["As,req top"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["As"] == pytest.approx(8.042, rel=1e-2)
    # M_Rd of the 4x16 provided: T = 804.2*434.78 = 349.7 kN, block depth
    # u = T/(eta*f_cd*b) = 123.41 mm, z = d - u/2 = 498.3 mm -> 174.24 kNm.
    assert results.iloc[1]["MRd"] == pytest.approx(174.24, rel=1e-3)
    assert results.iloc[1]["DCR"] == pytest.approx(0.861, rel=1e-3)


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
def test_check_flexure_ACI_318_19_3(beam_example_flexure_ACI: RectangularBeam) -> None:
    """Simple bending check: Mu pequeño → sección simple, no cae en doble armadura.

    Como es caso simple no lo tocan los fixes de doble armadura ni el displaced
    concrete, y el bloque de Whitney se rehace a mano (b = 12", d = 24 - 1.5 -
    0.375 - 1.41/2 = 21.42", A_s = 2#11 = 3.1229 in²):
      a = A_s·f_y/(0.85·f'c·b) = 3.1229·60/(0.85·4·12) = 4.5925"
      Mn = 3.1229·60·(21.42 - 4.5925/2) = 3583.3 kip·in ; ØMn = 0.9·Mn = 364.37 kN·m
      A_s,min = max(3√4000, 200)/60000·12·21.42 = 0.8568 in² = 5.53 cm²  (§9.6.1.2)
      A_s,req: 60·A_s·(21.42 - 0.7353·A_s) = 200·12/0.9 → A_s = 2.249 in² = 14.51 cm²

    Source: Calcpad "ACI 318-19 Beam Flexure 03_v3 - test_3.cpd"; the numbers above
    are the hand check of that sheet.
    """
    f = Forces(label="Test_03", M_y=200 * kip * ft)
    beam_example_flexure_ACI.set_longitudinal_rebar_bot(n1=2, d_b1=1.41 * inch)
    beam_example_flexure_ACI.set_longitudinal_rebar_top(n1=2, d_b1=0.75 * inch)
    node = Node(section=beam_example_flexure_ACI, forces=f)
    results = node.check_flexure()

    assert results.iloc[1]["Label"] == "B-12x24"
    assert results.iloc[1]["Comb."] == "Test_03"
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["As,min"] == pytest.approx(5.53, rel=1e-3)
    assert results.iloc[1]["As,req bot"] == pytest.approx(14.51, rel=1e-3)
    assert results.iloc[1]["As,req top"] == pytest.approx(0, rel=1e-3)
    assert results.iloc[1]["As"] == pytest.approx(20.15, rel=1e-3)
    assert results.iloc[1]["ØMn"] == pytest.approx(364.37, rel=1e-3)


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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 31 (Test_Etabs_05), column BC (A_s,min
    governs: column V = column R).
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
def test_maximum_flexural_reinforcement_ratio_ACI_318_19_Test_Etabs_05() -> None:
    """
    Test_Etabs_05: b=16", h=30", fc=6000psi, fy=60ksi.
    Agregado: 2026-05-03

    Se testea _maximum_flexural_reinforcement_ratio_ACI_318_19 directamente.
    ρ_max determina el umbral entre sección simple y doblemente armada —
    si esta fórmula falla, todo el diseño doble puede fallar silenciosamente.

    El caso usa fc=6000 psi donde β1=0.75 (no el 0.85 default para fc≤4000 psi),
    lo que ejercita el cálculo de β1 reducido.

    Excel col S: As_max = 10.4288 in²
    → ρ_max = As_max / (b × d) = 10.4288 / (16 × 27.5) = 0.02370

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 31 (Test_Etabs_05), column S (A_s,max =
    10.4288 in²).
    """
    from mento.codes.ACI_318_19_beam import _maximum_flexural_reinforcement_ratio_ACI_318_19

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

    rho_max = _maximum_flexural_reinforcement_ratio_ACI_318_19(beam)

    # Verificación directa contra Excel: As_max = ρ_max × b × d
    d = 27.5 * inch
    As_max = rho_max * beam.width * d
    assert As_max.to("inch**2").magnitude == pytest.approx(10.4288, rel=1e-3)


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

    Excel col R: As_min = 1.7041 in²
    → ρ_min = As_min / (b × d) = 1.7041 / (16 × 27.5) = 0.003873

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 31 (Test_Etabs_05), column R (A_s,min =
    1.7041 in²).
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 29 (Test_Etabs_03), columns BC (A_s =
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 29 (Test_Etabs_03), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 27 (Test_Etabs_01), columns BC/BE (A_s =
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 27 (Test_Etabs_01), columns BC and BE.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 49 (Test_Etabs_23), columns BC and BE.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 49 (Test_Etabs_23), columns BC/BE (A_s =
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 30 (Test_Etabs_04), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 28 (Test_Etabs_02), columns BC and BE.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 32 (Test_Etabs_06), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 37 (Test_Etabs_11), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 38 (Test_Etabs_12), columns BC and BE.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 39 (Test_Etabs_13), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 40 (Test_Etabs_14), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 41 (Test_Etabs_15), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 42 (Test_Etabs_16), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 43 (Test_Etabs_17), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 45 (Test_Etabs_19), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 46 (Test_Etabs_20), column BC.
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

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 48 (Test_Etabs_22), column BC.
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


@pytest.mark.published_example
def test_check_flexure_ACI_318_19_over_reinforced_but_top_redeems(
    beam_example_flexure_ACI: RectangularBeam,
) -> None:
    """
    Cubre la rama 2 (doble armadura valida gracias al aporte del top).
    A_s_bot > A_s_max_bot pero A_s_bot <= A_s_max_bot + A_s_top·f_s'_net/f_y.
    Es el caso tipico de diseño real: la seccion sola seria sobre-armada, pero
    el aporte del acero de compresion la mantiene ductil sin necesidad de cap.

    Sobre la fixture (b=12", h=24", fc=4000psi, fy=60ksi):
      Bot: 4×#10 (5.07 in² = 32.69 cm²)  >  A_s_max_bot ≈ 29.80 cm²
      Top: 2×#6  (0.88 in² = 5.70 cm²)   → A_s_max_total ≈ 35.17 cm²
      A_s_bot (32.69) <= A_s_max_total (35.17) → rama 2 → double_reinf con A_s real.

    Compatibilidad de deformaciones escrita aparte (ACI 318-19 §22.2, con el hormigón
    desplazado descontado en el acero comprimido), d = 21.49", d' = 2.25":
      c = 7.325" ; a = 0.85·c = 6.226" ; f's = f_y (la barra de arriba plastifica)
      eps_t = 0.003·(21.49 - 7.325)/7.325 = 0.0058 ≥ 0.005 → Ø = 0.90 (Tabla 21.2.2)
      Mn = 469.19 kip·ft = 636.13 kN·m ; ØMn = 572.52 kN·m
    Sin corrida de ETABS ni spColumn todavía.

    Source: Calcpad "ACI 318-19 Beam Flexure 03_v3 - over_reinforced_but_top_redeems.cpd";
    the strain compatibility above is the hand check of that sheet.
    """
    f = Forces(label="Test_over_but_top_redeems", M_y=400 * kip * ft)
    beam_example_flexure_ACI.set_longitudinal_rebar_bot(n1=4, d_b1=1.27 * inch)
    beam_example_flexure_ACI.set_longitudinal_rebar_top(n1=2, d_b1=0.75 * inch)
    node = Node(section=beam_example_flexure_ACI, forces=f)
    results = node.check_flexure()
    assert results.iloc[1]["Position"] == "Bottom"
    assert results.iloc[1]["Mu"] == pytest.approx(542.33, rel=1e-3)
    assert results.iloc[1]["ØMn"] == pytest.approx(572.52, rel=1e-3)
    # Past A_s_max, within A_s_max_eff: complies, doubly reinforced, and the
    # report's maximum is the one it is held to, not the singly reinforced one.
    min_max = beam_example_flexure_ACI._data_min_max_flexure
    assert min_max["Ok?"][2] == "✅ D.R."
    assert min_max["Max."][2] == pytest.approx(35.17, abs=0.02)
    assert "not_tension_controlled" not in {w.code for w in node.warnings}
