"""Tests for the public design results API (mento.design_results)."""

import math

import pytest

from mento import Concrete_ACI_318_19, Forces, Node, RectangularBeam, SteelBar
from mento import MPa, cm, kN, kNm, mm
from mento.design_results import (
    DesignNotRunError,
    FlexureDesign,
    RebarLayer,
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
# Results are a snapshot, not a live view
# ============================================================================


def test_results_are_frozen(designed_beam: RectangularBeam) -> None:
    from dataclasses import FrozenInstanceError

    with pytest.raises(FrozenInstanceError):
        designed_beam.flexure_design.bottom.A_s = 0 * cm**2  # type: ignore[misc]
