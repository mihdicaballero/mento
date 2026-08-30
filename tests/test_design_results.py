"""Tests for the public design results API (mento.design_results)."""

import math

import pytest
from pandas import DataFrame

from mento import Concrete_ACI_318_19, Forces, Node, RectangularBeam, SteelBar
from mento import MPa, cm, kN, kNm, m, mm
from mento.design_results import (
    DesignNotRunError,
    FlexureDesign,
    RebarLayer,
    SectionReinforcement,
    ShearDesign,
)

pytestmark = pytest.mark.filterwarnings("ignore::UserWarning")


@pytest.fixture
def beam() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    return RectangularBeam(
        label="101",
        concrete=concrete,
        steel_bar=steel,
        width=20 * cm,
        height=60 * cm,
        c_c=25 * mm,
    )


@pytest.fixture
def designed_beam(beam: RectangularBeam) -> RectangularBeam:
    node = Node(section=beam, forces=[Forces(label="C1", V_z=80 * kN, M_y=100 * kNm)])
    node.design()
    return beam


# ============================================================================
# RebarLayer
# ============================================================================


def test_rebar_layer_area_matches_the_circle_formula() -> None:
    layer = RebarLayer(n=3, d_b=16 * mm)
    assert layer.A_s.to("cm**2").magnitude == pytest.approx(3 * math.pi * (1.6**2) / 4)


# ============================================================================
# SectionReinforcement — what the section carries, readable without a check
# ============================================================================


def test_reinforcement_is_readable_before_any_check(beam: RectangularBeam) -> None:
    """The whole point of the view: no check has run, and it still answers.

    flexure_design and shear_design raise here; reinforcement must not, because
    it describes the section rather than a result.
    """
    with pytest.raises(DesignNotRunError):
        _ = beam.flexure_design
    with pytest.raises(DesignNotRunError):
        _ = beam.shear_design

    assert isinstance(beam.reinforcement, SectionReinforcement)


def test_reinforcement_reflects_the_longitudinal_setter(beam: RectangularBeam) -> None:
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)

    bottom = beam.reinforcement.bottom
    assert bottom.layers[0].n == 4
    assert bottom.layers[0].d_b.to("mm").magnitude == pytest.approx(16)
    assert bottom.n_bars == 4
    # 4 bars of 16 mm: 4*pi*1.6^2/4 = 8.04 cm²
    assert bottom.A_s.to("cm**2").magnitude == pytest.approx(4 * math.pi * (1.6**2) / 4)


def test_reinforcement_reflects_the_transverse_setter(beam: RectangularBeam) -> None:
    beam.set_transverse_rebar(n_stirrups=2, d_b=10 * mm, s_l=15 * cm)

    stirrups = beam.reinforcement.transverse
    assert stirrups.n_stirrups == 2
    assert stirrups.n_legs == 4
    assert stirrups.d_b.to("mm").magnitude == pytest.approx(10)
    assert stirrups.s_l.to("cm").magnitude == pytest.approx(15)
    assert stirrups.A_v.to("cm**2/m").magnitude > 0


def test_reinforcement_lists_no_layers_for_a_bare_face(beam: RectangularBeam) -> None:
    beam.set_longitudinal_rebar_top(n1=0, d_b1=0 * mm)
    assert beam.reinforcement.top.layers == ()
    assert beam.reinforcement.top.n_bars == 0
    assert str(beam.reinforcement.top) == "no reinforcement"


def test_reinforcement_is_a_value_not_a_live_view(beam: RectangularBeam) -> None:
    """It is a snapshot: changing the section afterwards must not change it."""
    before = beam.reinforcement.transverse
    beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=20 * cm)
    assert beam.reinforcement.transverse != before


def test_reinforcement_matches_the_design_after_a_design_runs(
    designed_beam: RectangularBeam,
) -> None:
    """The two views agree on what the section carries; they differ only in that
    the design also reports what the check demanded."""
    rebar = designed_beam.reinforcement
    flexure = designed_beam.flexure_design
    shear = designed_beam.shear_design

    assert rebar.bottom.layers == flexure.bottom.layers
    assert rebar.bottom.A_s == flexure.bottom.A_s
    assert rebar.top.layers == flexure.top.layers
    assert rebar.transverse.n_stirrups == shear.n_stirrups
    assert rebar.transverse.s_l == shear.s_l
    assert rebar.transverse.A_v == shear.A_v


def test_reinforcement_str_describes_both_faces_and_stirrups(beam: RectangularBeam) -> None:
    beam.set_longitudinal_rebar_bot(n1=3, d_b1=20 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=20 * cm)
    text = str(beam.reinforcement)
    assert "bottom: 3Ø20" in text
    assert "stirrups: 1eØ8" in text


def test_reinforcement_str_says_so_when_there_are_no_stirrups(beam: RectangularBeam) -> None:
    beam.set_transverse_rebar(n_stirrups=0, d_b=0 * mm, s_l=0 * cm)
    assert str(beam.reinforcement.transverse) == "no stirrups"


def test_rebar_layer_str_is_the_engineering_shorthand() -> None:
    assert str(RebarLayer(n=2, d_b=12 * mm)) == "2Ø12 mm"


# ============================================================================
# Reading results before running anything
# ============================================================================


def test_flexure_design_before_design_raises(beam: RectangularBeam) -> None:
    with pytest.raises(DesignNotRunError, match="No flexure results yet"):
        beam.flexure_design


def test_shear_design_before_design_raises(beam: RectangularBeam) -> None:
    with pytest.raises(DesignNotRunError, match="No shear results yet"):
        beam.shear_design


# ============================================================================
# Flexure
# ============================================================================


def test_flexure_design_reports_the_designed_bottom_steel(designed_beam: RectangularBeam) -> None:
    flexure = designed_beam.flexure_design

    assert isinstance(flexure, FlexureDesign)
    assert flexure.bottom.A_s == designed_beam._A_s_bot
    assert flexure.bottom.A_s_req == designed_beam._A_s_req_bot
    assert flexure.bottom.A_s >= flexure.bottom.A_s_req


def test_flexure_layers_match_the_private_layer_attributes(designed_beam: RectangularBeam) -> None:
    layers = designed_beam.flexure_design.bottom.layers

    assert layers, "the designed beam should carry bottom reinforcement"
    assert layers[0].n == designed_beam._n1_b
    assert layers[0].d_b == designed_beam._d_b1_b


def test_flexure_layer_areas_add_up_to_the_face_area(designed_beam: RectangularBeam) -> None:
    bottom = designed_beam.flexure_design.bottom
    total = sum(layer.A_s.to("cm**2").magnitude for layer in bottom.layers)

    assert total == pytest.approx(bottom.A_s.to("cm**2").magnitude)
    assert bottom.n_bars == sum(layer.n for layer in bottom.layers)


def test_empty_face_reports_no_layers(designed_beam: RectangularBeam) -> None:
    """A face with no bars has no layers, and says so."""
    top = designed_beam.flexure_design.top

    assert top.layers == ()
    assert top.n_bars == 0
    assert str(top) == "no reinforcement"


def test_flexure_dcr_is_the_worst_of_both_faces(designed_beam: RectangularBeam) -> None:
    flexure = designed_beam.flexure_design
    assert flexure.DCR == max(flexure.bottom.DCR, flexure.top.DCR)


def test_flexure_str_describes_both_faces(designed_beam: RectangularBeam) -> None:
    text = str(designed_beam.flexure_design)
    assert "bottom:" in text and "top:" in text


# ============================================================================
# Shear
# ============================================================================


def test_shear_design_reports_the_designed_stirrups(designed_beam: RectangularBeam) -> None:
    shear = designed_beam.shear_design

    assert isinstance(shear, ShearDesign)
    assert shear.n_stirrups == designed_beam._stirrup_n
    assert shear.d_b == designed_beam._stirrup_d_b
    assert shear.s_l == designed_beam._stirrup_s_l
    assert shear.A_v == designed_beam._A_v


def test_shear_legs_are_twice_the_stirrup_count(designed_beam: RectangularBeam) -> None:
    shear = designed_beam.shear_design
    assert shear.n_legs == shear.n_stirrups * 2


def test_shear_provided_area_covers_the_requirement(designed_beam: RectangularBeam) -> None:
    shear = designed_beam.shear_design
    assert shear.A_v >= shear.A_v_req
    assert shear.DCR <= 1.0


def test_shear_str_is_the_engineering_shorthand(designed_beam: RectangularBeam) -> None:
    shear = designed_beam.shear_design
    assert str(shear).startswith(f"{shear.n_stirrups}eØ")


def test_shear_str_without_stirrups(beam: RectangularBeam) -> None:
    """A beam checked with no stirrups reports that, rather than an empty layout."""
    node = Node(section=beam, forces=[Forces(label="C1", V_z=5 * kN)])
    node.check_shear()

    shear = beam.shear_design
    assert shear.n_stirrups == 0
    assert shear.n_legs == 0
    assert shear.A_v.magnitude == 0
    assert str(shear) == "no stirrups"


# ============================================================================
# Several load combinations: results are the envelope, not the last one checked
# ============================================================================


def _table_dcr(results: DataFrame) -> list[float]:
    """DCR of every combination in a results table, dropping the units row."""
    return [float(value) for value in results["DCR"][1:]]


# The tables round the DCR to three decimals; the results API keeps full precision.
TABLE_ROUNDING = 5e-4


def _two_combinations() -> list[Forces]:
    """Two combinations where the *first* governs both flexure and shear.

    The second one puts the top face in tension and is mild in shear, so it leaves
    the bottom face and the stirrup requirement at zero on its way out.
    """
    return [
        Forces(label="C1", V_z=90 * kN, M_y=100 * kNm),
        Forces(label="C2", V_z=20 * kN, M_y=-30 * kNm),
    ]


@pytest.fixture
def beam_two_combinations(beam: RectangularBeam) -> RectangularBeam:
    Node(section=beam, forces=_two_combinations()).design()
    return beam


def test_flexure_reports_the_combination_that_governs_each_face(
    beam_two_combinations: RectangularBeam,
) -> None:
    beam = beam_two_combinations
    worst = max(_table_dcr(beam.check_flexure(_two_combinations())))

    flexure = beam.flexure_design

    assert flexure.DCR == pytest.approx(worst, abs=TABLE_ROUNDING)
    assert flexure.bottom.DCR == pytest.approx(worst, abs=TABLE_ROUNDING)
    # The last combination checked required nothing at the bottom; the governing one
    # is what the result has to describe.
    assert flexure.bottom.A_s_req.magnitude > 0
    assert flexure.bottom.A_s_min.magnitude > 0
    assert flexure.bottom.A_s >= flexure.bottom.A_s_req


def test_shear_reports_the_combination_that_governs(
    beam_two_combinations: RectangularBeam,
) -> None:
    beam = beam_two_combinations
    worst = max(_table_dcr(beam.check_shear(_two_combinations())))

    shear = beam.shear_design

    assert shear.DCR == pytest.approx(worst, abs=TABLE_ROUNDING)
    assert shear.A_v_req.magnitude > 0
    assert shear.A_v >= shear.A_v_req


def test_envelope_skips_quantities_the_design_code_does_not_set() -> None:
    """A design code that leaves a quantity unset must not break the envelope.

    Both codes shipped today set all of them, so this is built directly rather than
    through a check.
    """
    beam = object.__new__(RectangularBeam)
    beam._flexure_envelope = {}
    beam._shear_envelope = {}
    beam._A_s_req_bot = 4 * cm**2  # _A_s_min_bot and _A_s_max_bot deliberately absent
    beam._DCRb_bot = 0.5
    beam._A_v_req = 3 * cm**2 / m
    beam._DCRv = 0.4

    beam._update_flexure_envelope("bot")
    beam._update_shear_envelope()

    assert beam._flexure_envelope["bot"] == {"A_s_req": 4 * cm**2, "DCR": 0.5}
    assert beam._shear_envelope == {"A_v_req": 3 * cm**2 / m, "DCR": 0.4}


def test_envelope_resets_between_checks(beam_two_combinations: RectangularBeam) -> None:
    """Checking a milder set of forces afterwards must not keep the earlier worst case."""
    beam = beam_two_combinations
    mild = [Forces(label="C3", V_z=10 * kN, M_y=10 * kNm)]

    beam.check_flexure(mild)
    beam.check_shear(mild)

    assert beam.flexure_design.DCR == pytest.approx(max(_table_dcr(beam.check_flexure(mild))), abs=TABLE_ROUNDING)
    assert beam.shear_design.DCR == pytest.approx(max(_table_dcr(beam.check_shear(mild))), abs=TABLE_ROUNDING)


# ============================================================================
# Results are a snapshot, not a live view
# ============================================================================


def test_results_are_frozen(designed_beam: RectangularBeam) -> None:
    from dataclasses import FrozenInstanceError

    with pytest.raises(FrozenInstanceError):
        designed_beam.flexure_design.bottom.A_s = 0 * cm**2  # type: ignore[misc]
