"""Design alternatives and repeatable designs.

``flexure_design.bottom.options`` / ``.top.options`` and ``shear_design.options``
expose what the rebar search ranked, with the layout it applied first. A design
also has to be a function of its inputs: running it again, or after the bars
were changed by hand, gives the same reinforcement.
"""

import pytest

from mento import (
    BeamSettings,
    Concrete_ACI_318_19,
    Concrete_EN_1992_2004,
    Forces,
    Node,
    OneWaySlab,
    RebarOption,
    RectangularBeam,
    SteelBar,
    StirrupOption,
)
from mento.rebar import Rebar
from mento.units import MPa, cm, kN, kNm, m, mm


def _aci_beam(settings: BeamSettings | None = None) -> RectangularBeam:
    return RectangularBeam(
        label="V101",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=20 * cm,
        height=60 * cm,
        c_c=25 * mm,
        settings=settings,
    )


def _designed(forces: list[Forces], settings: BeamSettings | None = None) -> RectangularBeam:
    beam = _aci_beam(settings)
    Node(section=beam, forces=forces).design()
    return beam


HIGH_SHEAR = [Forces(label="1.2D+1.6L", V_z=250 * kN, M_y=150 * kNm)]
TWO_FACES = [
    Forces(label="1.2D+1.6L", V_z=250 * kN, M_y=150 * kNm),
    Forces(label="1.4D", V_z=120 * kN, M_y=-40 * kNm),
]


# ---------------------------------------------------------------------------
# Longitudinal options
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("face", ["bottom", "top"])
def test_the_first_option_is_the_layout_applied(face: str) -> None:
    beam = _designed(TWO_FACES)
    options = getattr(beam.flexure_design, face).options
    applied = getattr(beam.reinforcement, face)

    assert options, "a designed face reports at least the layout it carries"
    assert all(isinstance(option, RebarOption) for option in options)
    assert options[0].layers == applied.layers
    assert options[0].A_s.to("cm**2").magnitude == pytest.approx(applied.A_s.to("cm**2").magnitude)


def test_options_are_ranked_by_their_functional_and_distinct() -> None:
    options = _designed(TWO_FACES).flexure_design.bottom.options

    assert len(options) == 3
    functionals = [option.functional for option in options]
    assert None not in functionals
    assert functionals == sorted(functionals)  # type: ignore[type-var]
    assert len({option.layers for option in options}) == len(options)


def test_an_option_area_is_the_area_of_its_bars() -> None:
    for option in _designed(TWO_FACES).flexure_design.bottom.options:
        expected = sum(layer.A_s.to("cm**2").magnitude for layer in option.layers)
        assert option.A_s.to("cm**2").magnitude == pytest.approx(expected)


@pytest.mark.parametrize("count", [1, 2, 5])
def test_the_number_of_options_is_a_setting(count: int) -> None:
    beam = _designed(TWO_FACES, BeamSettings(design_options=count))

    assert len(beam.flexure_design.bottom.options) == count
    assert 1 <= len(beam.shear_design.options) <= count


@pytest.mark.parametrize("bad", [0, -1, 2.5, True])
def test_the_number_of_options_must_be_a_positive_integer(bad: object) -> None:
    with pytest.raises(ValueError, match="design_options"):
        BeamSettings(design_options=bad)


def test_options_are_empty_before_any_design() -> None:
    beam = _aci_beam()
    beam.set_longitudinal_rebar_bot(n1=3, d_b1=16 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=10 * mm, s_l=20 * cm)
    Node(section=beam, forces=HIGH_SHEAR).check()

    assert beam.flexure_design.bottom.options == ()
    assert beam.shear_design.options == ()


def test_options_are_dropped_once_the_bars_are_changed_by_hand() -> None:
    beam = _designed(TWO_FACES)
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=25 * mm)
    beam.set_transverse_rebar(n_stirrups=2, d_b=12 * mm, s_l=10 * cm)

    assert beam.flexure_design.bottom.options == ()
    assert beam.flexure_design.top.options != ()  # untouched face keeps its own
    assert beam.shear_design.options == ()


def test_slab_options_read_as_spacings() -> None:
    slab = OneWaySlab(
        label="L1",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=SteelBar(name="ADN 420", f_y=420 * MPa),
        width=1 * m,
        height=20 * cm,
        c_c=25 * mm,
    )
    Node(section=slab, forces=[Forces(label="1.2D+1.6L", M_y=40 * kNm, V_z=50 * kN)]).design()

    options = slab.flexure_design.bottom.options
    assert options[0].layers == slab.reinforcement.bottom.layers
    assert all(layer.s is not None for option in options for layer in option.layers)
    assert len({option.layers for option in options}) == len(options)


# ---------------------------------------------------------------------------
# Stirrup options
# ---------------------------------------------------------------------------


def test_the_first_stirrup_option_is_the_one_applied() -> None:
    beam = _designed(TWO_FACES)
    design = beam.shear_design
    first = design.options[0]

    assert isinstance(first, StirrupOption)
    assert (first.n_stirrups, first.d_b, first.s_l) == (design.n_stirrups, design.d_b, design.s_l)
    assert first.A_v.to("cm**2/m").magnitude == pytest.approx(design.A_v.to("cm**2/m").magnitude)


def test_the_alternatives_are_the_same_cage_in_a_heavier_bar() -> None:
    """Where s_max governs, every diameter lands on the same spacing: the options are the bars."""
    options = _designed(HIGH_SHEAR).shear_design.options
    diameters = [option.d_b.to("mm").magnitude for option in options]

    assert diameters == [10, 12, 16]
    assert diameters[1:] == sorted(diameters[1:])
    assert {option.s_l for option in options} == {13 * cm}
    assert {option.n_stirrups for option in options} == {1}
    # The functional says what each heavier bar costs in steel.
    assert [round(option.functional, 2) for option in options] == [0.18, 0.70, 2.02]


def test_the_alternatives_are_ordered_by_diameter_when_the_demand_governs() -> None:
    """A heavier bar buys a wider spacing there, so the options are not all the same cage."""
    beam = _designed([Forces(label="1.2D+1.6L", V_z=450 * kN, M_y=150 * kNm)])
    options = beam.shear_design.options

    assert [(o.d_b.to("mm").magnitude, o.s_l.to("cm").magnitude) for o in options] == [(10, 7), (12, 10), (16, 13)]


def test_the_stirrup_table_offers_one_layout_per_diameter() -> None:
    rebar = Rebar(_aci_beam())

    def row(n: int, d_b: float, s_l: float) -> dict:
        A_v = n * 2 * 3.14159 * (d_b * mm) ** 2 / 4 / (s_l * cm)
        return {"n_stir": n, "d_b": d_b * mm, "s_l": s_l * cm, "s_w": 14 * cm, "A_v": A_v.to("cm**2/m")}

    table = rebar._rank_stirrup_options([row(1, 10, 13), row(1, 12, 13), row(2, 10, 20)], 10 * cm**2 / m)
    kept = [(r.n_stir, r.d_b.to("mm").magnitude, r.s_l.to("cm").magnitude) for r in table.itertuples()]

    # Nothing is dropped: a heavier bar at the same spacing stays on offer.
    assert kept == [(1, 10, 13), (1, 12, 13), (2, 10, 20)]
    # Fewest stirrups first, least steel among those: the row the design applies.
    assert kept[0] == (1, 10, 13)
    # The functional of the applied row is the smallest excess, and two stirrups carry +1.
    functional = list(table["functional"])
    assert functional[0] == min(functional)
    assert functional[2] > 1


# ---------------------------------------------------------------------------
# Repeatable designs
# ---------------------------------------------------------------------------


def _signature(beam: RectangularBeam) -> tuple:
    return (str(beam.reinforcement), beam.shear_design.s_w, beam.flexure_design.bottom.options)


@pytest.mark.parametrize(
    "forces",
    [
        [Forces(label="V", V_z=120 * kN)],
        [Forces(label="VM", V_z=120 * kN, M_y=100 * kNm)],
        TWO_FACES,
    ],
    ids=["shear-only", "shear-and-moment", "two-faces"],
)
def test_design_gives_the_same_result_every_time(forces: list[Forces]) -> None:
    beam = _aci_beam()
    node = Node(section=beam, forces=forces)
    node.design()
    first = _signature(beam)
    node.design()
    node.design()

    assert _signature(beam) == first
    assert _signature(_designed(forces)) == first


def test_design_ignores_the_reinforcement_set_before_it() -> None:
    beam = _aci_beam()
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=25 * mm)
    beam.set_longitudinal_rebar_top(n1=2, d_b1=20 * mm)
    beam.set_transverse_rebar(n_stirrups=2, d_b=16 * mm, s_l=8 * cm)
    Node(section=beam, forces=TWO_FACES).design()

    assert _signature(beam) == _signature(_designed(TWO_FACES))


def test_design_shear_alone_is_repeatable() -> None:
    beam = _aci_beam()
    node = Node(section=beam, forces=[Forces(label="V", V_z=120 * kN)])
    node.design_shear()
    first = str(beam.shear_design)
    node.design_shear()

    assert str(beam.shear_design) == first


def test_the_first_design_respects_the_spacing_limit_of_the_finished_beam() -> None:
    """The spacing used to be chosen at the depth of the 8 mm starter stirrup: 28 cm past a 27.95 cm limit."""
    beam = _designed([Forces(label="V", V_z=120 * kN)])

    assert beam.shear_design.s_l <= beam._d_shear / 2
    assert beam.shear_design.s_l == 27 * cm
    assert "stirrup_spacing_exceeds_max" not in [w.code for w in beam.warnings]


def test_en_design_is_repeatable_and_reports_options() -> None:
    def make() -> RectangularBeam:
        return RectangularBeam(
            label="V",
            concrete=Concrete_EN_1992_2004(name="C25", f_c=25 * MPa),
            steel_bar=SteelBar(name="B500S", f_y=500 * MPa),
            width=20 * cm,
            height=60 * cm,
            c_c=25 * mm,
        )

    beam = make()
    node = Node(section=beam, forces=TWO_FACES)
    node.design()
    first = _signature(beam)
    node.design()

    assert _signature(beam) == first
    assert beam.flexure_design.bottom.options[0].layers == beam.reinforcement.bottom.layers
    assert beam.shear_design.options[0].s_l == beam.shear_design.s_l


def test_options_read_like_the_reinforcement() -> None:
    beam = _designed(TWO_FACES)
    bottom = beam.flexure_design.bottom.options[0]
    stirrups = beam.shear_design.options[0]

    assert str(bottom) == str(beam.reinforcement.bottom)
    assert bottom.n_bars == beam.reinforcement.bottom.n_bars
    assert str(stirrups) == str(beam.shear_design)
    assert stirrups.n_legs == beam.shear_design.n_legs
    assert str(RebarOption(layers=(), A_s=0 * cm**2)) == "no reinforcement"
