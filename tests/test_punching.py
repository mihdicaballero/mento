import inspect
import math
from unittest.mock import patch

import pytest

from typing import Generator

from mento.codes.aci_318_19.equations import punching as aci_punching_eq
from mento.codes.en_1992_2004.equations import punching as en_punching_eq
from mento.codes.registry import design_code
from mento.column import Column
from mento.punching import Capital, Opening, PunchingNode, PunchingSlab
from mento.punching_results import PunchingCheck, PunchingCheckNotRunError, envelope_punching
from mento.forces import Forces
from mento.material import Concrete_ACI_318_19, Concrete_CIRSOC_201_25, Concrete_EN_1992_2004
from mento.units import cm, mm, kN, kNm, MPa, inch, psi, deg

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def reset_punching_node_id() -> Generator[None, None, None]:
    original = PunchingNode._last_id
    PunchingNode._last_id = 0
    yield
    PunchingNode._last_id = original


@pytest.fixture
def conc_aci() -> Concrete_ACI_318_19:
    return Concrete_ACI_318_19(name="C25", f_c=25 * MPa)


@pytest.fixture
def conc_en() -> Concrete_EN_1992_2004:
    return Concrete_EN_1992_2004(name="C25/30", f_c=25 * MPa)


@pytest.fixture
def slab(conc_aci, steel) -> PunchingSlab:
    """A 25 cm slab with Ø12/15 base mat plus Ø16/15 extra bars over the column in x."""
    slab = PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25 * mm)
    slab.set_rebar_x(d_b1=12 * mm, s_b1=15 * cm, d_b3=16 * mm, s_b3=15 * cm)
    slab.set_rebar_y(d_b1=12 * mm, s_b1=15 * cm)
    return slab


@pytest.fixture
def bare_slab(conc_aci, steel) -> PunchingSlab:
    """The same slab with no reinforcement declared — depths from the estimate."""
    return PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25 * mm)


@pytest.fixture
def col_interior() -> Column:
    return Column(shape="rectangular", position="interior", b=40 * cm, h=40 * cm)


@pytest.fixture
def col_edge() -> Column:
    return Column(shape="rectangular", position="edge", b=40 * cm, h=40 * cm, edge_distance_x=20 * cm)


@pytest.fixture
def col_corner() -> Column:
    return Column(
        shape="rectangular", position="corner", b=35 * cm, h=35 * cm, edge_distance_x=18 * cm, edge_distance_y=18 * cm
    )


@pytest.fixture
def f1() -> Forces:
    return Forces(label="ELU 1", V_z=500 * kN)


# ---------------------------------------------------------------------------
# Column
# ---------------------------------------------------------------------------


class TestColumn:
    def test_rectangular_interior(self):
        col = Column(shape="rectangular", position="interior", b=40 * cm, h=40 * cm)
        assert col.shape == "rectangular"
        assert col.position == "interior"
        assert col.b.to("cm").magnitude == pytest.approx(40)
        assert col.h.to("cm").magnitude == pytest.approx(40)
        assert col.edge_distance_x is None
        assert col.edge_distance_y is None

    def test_circular_interior(self):
        col = Column(shape="circular", position="interior", b=50 * cm)
        assert col.shape == "circular"
        assert col.b.to("cm").magnitude == pytest.approx(50)

    def test_edge_requires_edge_distance_x(self):
        with pytest.raises(ValueError, match="edge_distance_x is required"):
            Column(shape="rectangular", position="edge", b=40 * cm, h=40 * cm)

    def test_edge_stores_edge_distance_x(self, col_edge):
        assert col_edge.edge_distance_x.to("cm").magnitude == pytest.approx(20)
        assert col_edge.edge_distance_y is None

    def test_corner_requires_both_distances(self):
        with pytest.raises(ValueError, match="edge_distance_y is required"):
            Column(shape="rectangular", position="corner", b=35 * cm, h=35 * cm, edge_distance_x=18 * cm)

    def test_corner_requires_edge_distance_x(self):
        with pytest.raises(ValueError, match="edge_distance_x is required"):
            Column(shape="rectangular", position="corner", b=35 * cm, h=35 * cm)

    def test_corner_stores_both_distances(self, col_corner):
        assert col_corner.edge_distance_x.to("cm").magnitude == pytest.approx(18)
        assert col_corner.edge_distance_y.to("cm").magnitude == pytest.approx(18)

    def test_invalid_shape(self):
        with pytest.raises(ValueError, match="shape must be"):
            Column(shape="triangular", position="interior", b=40 * cm, h=40 * cm)  # type: ignore[arg-type]

    def test_invalid_position(self):
        with pytest.raises(ValueError, match="position must be"):
            Column(shape="rectangular", position="central", b=40 * cm, h=40 * cm)  # type: ignore[arg-type]

    def test_repr_rectangular(self, col_interior):
        r = repr(col_interior)
        assert "Rectangular" in r
        assert "Interior" in r

    def test_repr_circular(self):
        col = Column(shape="circular", position="interior", b=50 * cm)
        r = repr(col)
        assert "Circular" in r
        assert "Interior" in r

    def test_default_b_h_zero(self):
        col = Column(shape="rectangular", position="interior")
        assert col.b.to("cm").magnitude == pytest.approx(0)
        assert col.h.to("cm").magnitude == pytest.approx(0)


# ---------------------------------------------------------------------------
# PunchingSlab
# ---------------------------------------------------------------------------


class TestPunchingSlab:
    # --- Effective depths derived from the declared bars -------------------

    def test_d_avg_estimate_without_rebar(self, bare_slab):
        # Two mats of the Ø16 estimate: d_avg = h - c_c - 16 mm = 209 mm
        assert bare_slab.d_avg.to("mm").magnitude == pytest.approx(209.0)

    def test_depths_without_rebar_straddle_the_estimate(self, bare_slab):
        # x outside y: 250 - 25 - 8 and 250 - 25 - 16 - 8
        assert bare_slab.d_x.to("mm").magnitude == pytest.approx(217.0)
        assert bare_slab.d_y.to("mm").magnitude == pytest.approx(201.0)

    def test_depths_from_declared_bars(self, slab):
        # Ø16 governs in x, Ø12 in y; x is the outer mat.
        assert slab.d_x.to("mm").magnitude == pytest.approx(217.0)  # 250 - 25 - 8
        assert slab.d_y.to("mm").magnitude == pytest.approx(203.0)  # 250 - 25 - 16 - 6
        assert slab.d_avg.to("mm").magnitude == pytest.approx(210.0)

    def test_outer_direction_swaps_the_layers(self, conc_aci, steel):
        slab = PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25 * mm, outer_direction="y")
        slab.set_rebar_x(d_b1=16 * mm, s_b1=15 * cm)
        slab.set_rebar_y(d_b1=12 * mm, s_b1=15 * cm)
        assert slab.d_y.to("mm").magnitude == pytest.approx(219.0)  # 250 - 25 - 6
        assert slab.d_x.to("mm").magnitude == pytest.approx(205.0)  # 250 - 25 - 12 - 8

    def test_largest_bar_governs_the_depth_of_a_mat(self, conc_aci, steel):
        slab = PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25 * mm)
        slab.set_rebar_x(d_b1=10 * mm, s_b1=20 * cm, d_b3=20 * mm, s_b3=20 * cm)
        assert slab.d_x.to("mm").magnitude == pytest.approx(215.0)  # 250 - 25 - 10

    def test_invalid_outer_direction(self, conc_aci, steel):
        with pytest.raises(ValueError, match="outer_direction must be 'x' or 'y'"):
            PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25 * mm, outer_direction="z")  # type: ignore[arg-type]

    # --- Reinforcement ratios derived from the same bars -------------------

    def test_A_s_sums_base_mat_and_extra_bars(self, slab):
        # (pi*12^2/4 + pi*16^2/4) / 150 mm = 2.0944 mm2/mm = 20.944 cm2/m
        assert slab.A_s_x.to("cm**2/m").magnitude == pytest.approx(20.944, rel=1e-3)
        assert slab.A_s_y.to("cm**2/m").magnitude == pytest.approx(7.540, rel=1e-3)

    def test_rho_is_A_s_over_d(self, slab):
        assert slab.rho_x == pytest.approx(20.944e-4 / 0.217, rel=1e-3)
        assert slab.rho_y == pytest.approx(7.540e-4 / 0.203, rel=1e-3)

    def test_rho_l_is_the_geometric_mean(self, slab):
        assert slab.rho_l == pytest.approx(math.sqrt(slab.rho_x * slab.rho_y))

    def test_rho_is_none_without_rebar(self, bare_slab):
        assert bare_slab.rho_x is None
        assert bare_slab.rho_y is None
        assert bare_slab.rho_l is None
        assert bare_slab.A_s_x is None
        assert bare_slab.has_rebar is False

    def test_rho_l_needs_both_directions(self, conc_aci, steel):
        slab = PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25 * mm)
        slab.set_rebar_x(d_b1=12 * mm, s_b1=15 * cm)
        assert slab.rho_x is not None
        assert slab.rho_y is None
        assert slab.rho_l is None
        assert slab.has_rebar is False

    def test_extra_bars_raise_rho(self, conc_aci, steel):
        base_only = PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25 * mm)
        base_only.set_rebar_x(d_b1=12 * mm, s_b1=15 * cm)
        reinforced = PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25 * mm)
        reinforced.set_rebar_x(d_b1=12 * mm, s_b1=15 * cm, d_b3=16 * mm, s_b3=15 * cm)
        assert reinforced.rho_x > base_only.rho_x

    def test_set_rebar_replaces_rather_than_merges(self, slab):
        slab.set_rebar_x(d_b1=12 * mm, s_b1=15 * cm)
        assert slab.A_s_x.to("cm**2/m").magnitude == pytest.approx(7.540, rel=1e-3)

    # --- No back door between rho and d ------------------------------------

    def test_assigning_d_avg_is_refused(self, slab):
        with pytest.raises(AttributeError, match="set_effective_depth"):
            slab.d_avg = 19.8 * cm

    def test_set_effective_depth_rederives_rho(self, slab):
        rho_x_before = slab.rho_x
        slab.set_effective_depth(d_x=20.7 * cm, d_y=19.3 * cm)
        assert slab.d_avg.to("cm").magnitude == pytest.approx(20.0)
        assert slab.rho_x == pytest.approx(slab.A_s_x.to("cm**2/m").magnitude * 1e-4 / 0.207, rel=1e-3)
        assert slab.rho_x > rho_x_before  # a shallower d raises the ratio

    def test_set_effective_depth_one_direction_only(self, slab):
        slab.set_effective_depth(d_y=19.0 * cm)
        assert slab.d_x.to("mm").magnitude == pytest.approx(217.0)
        assert slab.d_y.to("cm").magnitude == pytest.approx(19.0)

    def test_set_effective_depth_needs_a_depth(self, slab):
        with pytest.raises(ValueError, match="needs d_x, d_y, or both"):
            slab.set_effective_depth()

    def test_set_effective_depth_rejects_non_length(self, slab):
        with pytest.raises(TypeError, match="d_x must be a length Quantity"):
            slab.set_effective_depth(d_x=20.7)

    def test_set_effective_depth_rejects_non_positive(self, slab):
        with pytest.raises(ValueError, match="d_y must be greater than zero"):
            slab.set_effective_depth(d_y=0 * cm)

    def test_rho_is_not_a_constructor_argument(self, conc_aci, steel):
        with pytest.raises(TypeError):
            PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25 * mm, rho_x=0.01)  # type: ignore[call-arg]

    # --- Imperial ----------------------------------------------------------

    def test_d_avg_imperial_estimate(self, steel):
        conc_imp = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
        slab = PunchingSlab(concrete=conc_imp, steel_bar=steel, h=10 * inch, c_c=1 * inch)
        # d_avg = 10 - 1 - 0.625 = 8.375 inch
        assert slab.d_avg.to("inch").magnitude == pytest.approx(8.375)
        assert slab.d_x.to("inch").magnitude == pytest.approx(8.6875)
        assert slab.d_y.to("inch").magnitude == pytest.approx(8.0625)

    def test_imperial_rho_from_bars(self, steel):
        conc_imp = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
        slab = PunchingSlab(concrete=conc_imp, steel_bar=steel, h=10 * inch, c_c=1 * inch)
        slab.set_rebar_x(d_b1=0.625 * inch, s_b1=12 * inch)
        assert slab.A_s_x.to("inch**2/ft").magnitude == pytest.approx(0.3068, rel=1e-3)
        assert slab.rho_x == pytest.approx(0.0255663 / 8.6875, rel=1e-3)

    # --- Bar declaration validation ----------------------------------------

    def test_diameter_without_spacing_is_refused(self, bare_slab):
        with pytest.raises(ValueError, match="position 1 in x need a spacing"):
            bare_slab.set_rebar_x(d_b1=12 * mm)

    def test_spacing_without_diameter_is_refused(self, bare_slab):
        with pytest.raises(ValueError, match="position 3 in y need a diameter"):
            bare_slab.set_rebar_y(d_b1=12 * mm, s_b1=15 * cm, s_b3=15 * cm)

    def test_rejects_non_length_diameter(self, bare_slab):
        with pytest.raises(TypeError, match="d_b1 must be a length Quantity"):
            bare_slab.set_rebar_x(d_b1=12, s_b1=15 * cm)

    def test_rejects_negative_spacing(self, bare_slab):
        with pytest.raises(ValueError, match="s_b1 must be greater than or equal to zero"):
            bare_slab.set_rebar_x(d_b1=12 * mm, s_b1=-15 * cm)

    def test_a_refused_declaration_leaves_the_slab_untouched(self, slab):
        with pytest.raises(ValueError):
            slab.set_rebar_x(d_b1=25 * mm)
        assert slab.A_s_x.to("cm**2/m").magnitude == pytest.approx(20.944, rel=1e-3)

    # --- Geometry validation (unchanged) ------------------------------------

    def test_unit_system_metric(self, slab):
        assert slab.unit_system == "metric"

    def test_unit_system_from_concrete(self, conc_en, steel):
        slab_en = PunchingSlab(concrete=conc_en, steel_bar=steel, h=25 * cm, c_c=25 * mm)
        assert slab_en.unit_system == "metric"

    def test_h_and_cc_stored(self, slab):
        assert slab.h.to("cm").magnitude == pytest.approx(25)
        assert slab.c_c.to("mm").magnitude == pytest.approx(25)

    def test_d_avg_thicker_slab(self, conc_aci, steel):
        slab = PunchingSlab(concrete=conc_aci, steel_bar=steel, h=30 * cm, c_c=30 * mm)
        # 300mm - 30mm - 16mm = 254mm
        assert slab.d_avg.to("mm").magnitude == pytest.approx(254.0)

    def test_rejects_non_positive_thickness(self, conc_aci, steel):
        with pytest.raises(ValueError, match="h must be greater than zero"):
            PunchingSlab(concrete=conc_aci, steel_bar=steel, h=0 * cm, c_c=25 * mm)

    def test_rejects_negative_cover(self, conc_aci, steel):
        with pytest.raises(ValueError, match="c_c must be greater than or equal to zero"):
            PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=-10 * mm)

    def test_rejects_cover_that_consumes_the_depth(self, conc_aci, steel):
        with pytest.raises(ValueError, match="d_avg must be greater than zero"):
            PunchingSlab(concrete=conc_aci, steel_bar=steel, h=10 * cm, c_c=9 * cm)

    def test_rejects_rebar_that_consumes_the_depth(self, conc_aci, steel):
        slab = PunchingSlab(concrete=conc_aci, steel_bar=steel, h=12 * cm, c_c=8 * cm)
        with pytest.raises(ValueError, match="d_avg must be greater than zero"):
            slab.set_rebar_x(d_b1=32 * mm, s_b1=15 * cm)

    def test_rejects_non_length_thickness(self, conc_aci, steel):
        with pytest.raises(TypeError, match="h must be a length Quantity"):
            PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25, c_c=25 * mm)

    def test_rejects_non_length_cover(self, conc_aci, steel):
        with pytest.raises(TypeError, match="c_c must be a length Quantity"):
            PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25)


class TestOpening:
    def test_rectangular_opening(self):
        op = Opening(shape="rectangular", x=80 * cm, y=60 * cm, b=40 * cm, h=40 * cm)
        assert op.shape == "rectangular"
        assert op.x.to("cm").magnitude == pytest.approx(80)
        assert op.y.to("cm").magnitude == pytest.approx(60)
        assert op.b.to("cm").magnitude == pytest.approx(40)
        assert op.h.to("cm").magnitude == pytest.approx(40)

    def test_circular_opening(self):
        op = Opening(shape="circular", x=50 * cm, y=50 * cm, diameter=30 * cm)
        assert op.shape == "circular"
        assert op.diameter.to("cm").magnitude == pytest.approx(30)

    def test_negative_offsets(self):
        op = Opening(shape="rectangular", x=-80 * cm, y=-60 * cm, b=40 * cm, h=40 * cm)
        assert op.x.to("cm").magnitude == pytest.approx(-80)
        assert op.y.to("cm").magnitude == pytest.approx(-60)

    def test_invalid_shape(self):
        with pytest.raises(ValueError, match="Opening shape must be"):
            Opening(shape="oval", x=50 * cm, y=50 * cm)  # type: ignore[arg-type]

    def test_default_dimensions_zero(self):
        op = Opening(shape="rectangular", x=50 * cm, y=50 * cm)
        assert op.b.to("cm").magnitude == pytest.approx(0)
        assert op.h.to("cm").magnitude == pytest.approx(0)
        assert op.diameter.to("cm").magnitude == pytest.approx(0)


# ---------------------------------------------------------------------------
# Capital
# ---------------------------------------------------------------------------


class TestCapital:
    def test_basic_construction(self):
        cap = Capital(b=100 * cm, h=100 * cm, thickness=25 * cm)
        assert cap.b.to("cm").magnitude == pytest.approx(100)
        assert cap.h.to("cm").magnitude == pytest.approx(100)
        assert cap.thickness.to("cm").magnitude == pytest.approx(25)

    def test_rectangular_capital(self):
        cap = Capital(b=120 * cm, h=80 * cm, thickness=30 * cm)
        assert cap.b.to("cm").magnitude == pytest.approx(120)
        assert cap.h.to("cm").magnitude == pytest.approx(80)

    def test_zero_thickness_raises(self):
        with pytest.raises(ValueError, match="Capital thickness must be positive"):
            Capital(b=100 * cm, h=100 * cm, thickness=0 * cm)

    def test_negative_thickness_raises(self):
        with pytest.raises(ValueError, match="Capital thickness must be positive"):
            Capital(b=100 * cm, h=100 * cm, thickness=-5 * cm)


# ---------------------------------------------------------------------------
# PunchingNode
# ---------------------------------------------------------------------------


class TestPunchingNode:
    def test_id_auto_increments(self, slab, col_interior, f1):
        n1 = PunchingNode(slab=slab, column=col_interior, forces=f1)
        n2 = PunchingNode(slab=slab, column=col_interior, forces=f1)
        assert n1.id == 1
        assert n2.id == 2

    def test_id_read_only(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        with pytest.raises(AttributeError):
            node.id = 99  # type: ignore[misc]

    def test_single_force_wrapped_in_list(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        assert isinstance(node.forces, list)
        assert len(node.forces) == 1
        assert node.forces[0] is f1

    def test_list_of_forces(self, slab, col_interior):
        fa = Forces(label="1.4D", V_z=400 * kN)
        fb = Forces(label="1.2D+1.6L", V_z=500 * kN)
        node = PunchingNode(slab=slab, column=col_interior, forces=[fa, fb])
        assert len(node.forces) == 2

    def test_default_no_openings(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        assert node.openings == []

    def test_default_no_capital(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        assert node.capital is None

    def test_with_openings(self, slab, col_interior, f1):
        op = Opening(shape="rectangular", x=80 * cm, y=60 * cm, b=40 * cm, h=40 * cm)
        node = PunchingNode(slab=slab, column=col_interior, forces=f1, openings=[op])
        assert len(node.openings) == 1
        assert node.openings[0] is op

    def test_with_capital(self, slab, col_interior, f1):
        cap = Capital(b=100 * cm, h=100 * cm, thickness=25 * cm)
        node = PunchingNode(slab=slab, column=col_interior, forces=f1, capital=cap)
        assert node.capital is cap

    def test_edge_column_node(self, slab, col_edge, f1):
        node = PunchingNode(slab=slab, column=col_edge, forces=f1)
        assert node.column.position == "edge"

    def test_corner_column_node(self, slab, col_corner, f1):
        cap = Capital(b=100 * cm, h=100 * cm, thickness=25 * cm)
        op = Opening(shape="rectangular", x=-80 * cm, y=-60 * cm, b=40 * cm, h=40 * cm)
        node = PunchingNode(slab=slab, column=col_corner, forces=f1, openings=[op], capital=cap)
        assert node.column.position == "corner"
        assert node.capital is not None
        assert len(node.openings) == 1

    def test_check_raises_not_implemented(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        with pytest.raises(NotImplementedError):
            node.check()

    def test_design_raises_not_implemented(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        with pytest.raises(NotImplementedError):
            node.design()

    def test_repr_contains_position_and_shape(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        r = repr(node)
        assert "PunchingNode" in r
        assert "Interior" in r
        assert "Rectangular" in r

    def test_repr_with_capital_and_openings(self, slab, col_corner, f1):
        cap = Capital(b=100 * cm, h=100 * cm, thickness=25 * cm)
        op = Opening(shape="circular", x=60 * cm, y=60 * cm, diameter=30 * cm)
        node = PunchingNode(slab=slab, column=col_corner, forces=f1, openings=[op], capital=cap)
        r = repr(node)
        assert "openings=1" in r
        assert "capital=" in r

    def test_biaxial_forces(self, slab, col_interior):
        f = Forces(label="ELU", V_z=300 * kN, M_y=50 * kNm, M_x=30 * kNm)
        node = PunchingNode(slab=slab, column=col_interior, forces=f)
        assert node.forces[0]._M_x.to("kN*m").magnitude == pytest.approx(30)

    # --- Plot (smoke tests — just verify no exception is raised) ---

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_plot_interior_rectangular(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        node.plot()

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_plot_edge_column(self, slab, col_edge, f1):
        node = PunchingNode(slab=slab, column=col_edge, forces=f1)
        node.plot()

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_plot_corner_with_capital_and_opening(self, slab, col_corner, f1):
        cap = Capital(b=100 * cm, h=100 * cm, thickness=25 * cm)
        op = Opening(shape="rectangular", x=-80 * cm, y=-60 * cm, b=40 * cm, h=40 * cm)
        node = PunchingNode(slab=slab, column=col_corner, forces=f1, openings=[op], capital=cap)
        node.plot()

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_plot_circular_column(self, slab, f1):
        col = Column(shape="circular", position="interior", b=50 * cm)
        node = PunchingNode(slab=slab, column=col, forces=f1)
        node.plot()

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_plot_circular_opening(self, slab, col_interior, f1):
        op = Opening(shape="circular", x=70 * cm, y=70 * cm, diameter=30 * cm)
        node = PunchingNode(slab=slab, column=col_interior, forces=f1, openings=[op])
        node.plot()


# ---------------------------------------------------------------------------
# Markdown views
# ---------------------------------------------------------------------------


class TestPunchingData:
    def test_slab_data_reports_bars_depths_and_ratios(self, slab):
        with patch("mento.reports.punching.display") as mock_display:
            assert slab.data is None
        mock_display.assert_called_once()
        md = slab._md_data
        assert "Ø12/15.0 cm ++ Ø16/15.0 cm" in md
        assert "21.7 cm" in md and "20.3 cm" in md  # d_x, d_y
        assert f"{slab.rho_x:.4f}" in md
        assert slab.concrete.name in md

    def test_slab_data_says_when_no_rebar_is_declared(self, bare_slab):
        with patch("mento.reports.punching.display"):
            bare_slab.data
        md = bare_slab._md_data
        assert "not declared" in md
        assert "set_rebar_x()" in md
        assert "—" in md  # rho renders as a dash rather than a zero

    def test_slab_data_imperial(self, steel):
        conc_imp = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
        slab = PunchingSlab(concrete=conc_imp, steel_bar=steel, h=10 * inch, c_c=1 * inch)
        slab.set_rebar_x(d_b1=0.625 * inch, s_b1=12 * inch)
        with patch("mento.reports.punching.display"):
            slab.data
        assert "in" in slab._md_data
        assert "cm" not in slab._md_data

    def test_node_data_shows_column_slab_and_forces(self, slab, col_edge):
        f = Forces(label="ELU", V_z=300 * kN, M_y=50 * kNm, M_x=30 * kNm)
        node = PunchingNode(slab=slab, column=col_edge, forces=f)
        with patch("mento.reports.punching.display") as mock_display:
            assert node.data is None
        rendered = " ".join(str(call.args[0].data) for call in mock_display.call_args_list)
        assert "edge" in rendered and "rectangular" in rendered
        assert "free edge" in rendered
        assert "ELU" in rendered
        assert "M_{x}" in rendered and "M_{y}" in rendered
        # The punching demand is Vu / VEd, the vertical load at the connection.
        assert "$V$=300.0 kN" in rendered
        assert "$N$" not in rendered
        assert slab._md_data in rendered

    def test_node_data_omits_absent_moments(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        with patch("mento.reports.punching.display") as mock_display:
            node.data
        rendered = " ".join(str(call.args[0].data) for call in mock_display.call_args_list)
        assert "M_{x}" not in rendered
        assert "M_{y}" not in rendered

    def test_node_data_with_capital_and_openings(self, slab, col_corner, f1):
        cap = Capital(b=100 * cm, h=100 * cm, thickness=25 * cm)
        op = Opening(shape="circular", x=60 * cm, y=60 * cm, diameter=30 * cm)
        node = PunchingNode(slab=slab, column=col_corner, forces=f1, openings=[op], capital=cap)
        with patch("mento.reports.punching.display") as mock_display:
            node.data
        rendered = " ".join(str(call.args[0].data) for call in mock_display.call_args_list)
        assert "Capital" in rendered
        assert "1 opening(s)" in rendered

    def test_node_data_circular_column(self, slab, f1):
        col = Column(shape="circular", position="interior", b=50 * cm)
        node = PunchingNode(slab=slab, column=col, forces=f1)
        with patch("mento.reports.punching.display") as mock_display:
            node.data
        rendered = " ".join(str(call.args[0].data) for call in mock_display.call_args_list)
        assert "D$=50.0 cm" in rendered
        assert "free edge" not in rendered


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------


def _check(label: str, dcr: float, v_c: float = 1.0) -> PunchingCheck:
    return PunchingCheck(
        label=label,
        b_0=200 * cm,
        d=21 * cm,
        v_u=dcr * v_c * MPa,
        v_c=v_c * MPa,
        DCR=dcr,
    )


class TestPunchingCheckEnvelope:
    def test_envelope_is_the_worst_combination(self):
        env = envelope_punching([_check("a", 0.40), _check("b", 0.90), _check("c", 0.55)])
        assert env.DCR == pytest.approx(0.90)
        assert env.label == "envelope"

    def test_envelope_carries_the_capacity_its_dcr_was_formed_from(self):
        env = envelope_punching([_check("a", 0.90, v_c=1.0), _check("b", 0.40, v_c=2.0)])
        assert env.v_c.to("MPa").magnitude == pytest.approx(1.0)
        assert (env.v_u / env.DCR).to("MPa").magnitude == pytest.approx(env.v_c.to("MPa").magnitude)

    def test_a_tie_on_dcr_breaks_to_the_lowest_resistance(self):
        env = envelope_punching([_check("a", 0.0, v_c=2.0), _check("b", 0.0, v_c=1.0)])
        assert env.v_c.to("MPa").magnitude == pytest.approx(1.0)

    def test_envelope_needs_a_combination(self):
        with pytest.raises(ValueError, match="at least one checked combination"):
            envelope_punching([])

    def test_result_is_frozen(self):
        with pytest.raises(AttributeError):
            _check("a", 0.5).DCR = 0.9  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Check skeleton: dispatch and preconditions (the equations are Phase 2 / 5)
# ---------------------------------------------------------------------------


class TestPunchingCheckWiring:
    @pytest.mark.parametrize(
        "concrete",
        [
            Concrete_ACI_318_19(name="C25", f_c=25 * MPa),
            Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa),
            Concrete_EN_1992_2004(name="C25/30", f_c=25 * MPa),
        ],
        ids=lambda c: c.design_code,
    )
    def test_every_code_registers_a_punching_check(self, concrete):
        assert design_code(concrete).check_punching is not None

    def test_check_needs_a_force(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        node.forces = []
        with pytest.raises(ValueError, match="at least one Forces"):
            node.check()

    def test_results_are_not_readable_before_check(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        with pytest.raises(PunchingCheckNotRunError, match=r"call check\(\) first"):
            node.punching_checks

    def test_aci_check_reaches_the_missing_equations(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        with pytest.raises(NotImplementedError, match="Phase 2"):
            node.check()

    def test_en_check_reaches_the_missing_equations(self, conc_en, steel, col_interior, f1):
        slab_en = PunchingSlab(concrete=conc_en, steel_bar=steel, h=25 * cm, c_c=25 * mm)
        slab_en.set_rebar_x(d_b1=12 * mm, s_b1=15 * cm)
        slab_en.set_rebar_y(d_b1=12 * mm, s_b1=15 * cm)
        node = PunchingNode(slab=slab_en, column=col_interior, forces=f1)
        with pytest.raises(NotImplementedError, match="Phase 5"):
            node.check()

    def test_a_force_without_vertical_load_is_refused_by_name(self, slab, col_interior):
        node = PunchingNode(slab=slab, column=col_interior, forces=Forces(label="ELU", M_y=50 * kNm))
        with pytest.raises(ValueError, match="V_z"):
            node.check()

    def test_en_refuses_a_slab_with_no_declared_rho(self, conc_en, steel, col_interior, f1):
        bare_en = PunchingSlab(concrete=conc_en, steel_bar=steel, h=25 * cm, c_c=25 * mm)
        node = PunchingNode(slab=bare_en, column=col_interior, forces=f1)
        with pytest.raises(ValueError, match="set_rebar_x"):
            node.check()

    def test_aci_does_not_need_rho(self, bare_slab, col_interior, f1):
        """ACI's v_c never reads rho, so an undeclared one must not block the check."""
        node = PunchingNode(slab=bare_slab, column=col_interior, forces=f1)
        with pytest.raises(NotImplementedError, match="Phase 2"):
            node.check()

    def test_a_capital_is_refused_until_phase_3(self, slab, col_interior, f1):
        cap = Capital(b=100 * cm, h=100 * cm, thickness=25 * cm)
        node = PunchingNode(slab=slab, column=col_interior, forces=f1, capital=cap)
        with pytest.raises(NotImplementedError, match="Phase 3"):
            node.check()

    def test_an_opening_is_refused_until_phase_3(self, slab, col_interior, f1):
        op = Opening(shape="rectangular", x=60 * cm, y=0 * cm, b=40 * cm, h=40 * cm)
        node = PunchingNode(slab=slab, column=col_interior, forces=f1, openings=[op])
        with pytest.raises(NotImplementedError, match="Phase 3"):
            node.check()

    def test_the_node_keeps_every_combination_and_returns_their_envelope(self, slab, col_interior, monkeypatch):
        """The wiring around the equations, exercised without them.

        Stands in for the checker so that storing the per-combination results and
        enveloping them is covered now rather than waiting on Phase 2.
        """
        results = iter([_check("a", 0.40), _check("b", 0.90)])

        class _StubCode:
            def requires(self, hook: str):
                assert hook == "check_punching"
                return lambda node, force: next(results)

        monkeypatch.setattr("mento.punching.design_code", lambda concrete: _StubCode())

        fa = Forces(label="a", V_z=400 * kN)
        fb = Forces(label="b", V_z=900 * kN)
        node = PunchingNode(slab=slab, column=col_interior, forces=[fa, fb])

        envelope = node.check()
        assert envelope.DCR == pytest.approx(0.90)
        assert [c.label for c in node.punching_checks] == ["a", "b"]

    def test_design_names_the_code_that_lacks_it(self, slab, col_interior, f1):
        node = PunchingNode(slab=slab, column=col_interior, forces=f1)
        with pytest.raises(NotImplementedError, match="design punching is not implemented"):
            node.design()


# ---------------------------------------------------------------------------
# Equation stubs
# ---------------------------------------------------------------------------


def _public_equations(module) -> list:
    return [
        getattr(module, name)
        for name in dir(module)
        if not name.startswith("_") and inspect.isfunction(getattr(module, name))
    ]


class TestEquationStubs:
    """The equations are recreated from a validated Calcpad sheet, one at a time.

    Until each lands, its stub must raise rather than return a wrong number, and
    these tests shrink as the module fills in.
    """

    @pytest.mark.parametrize("module", [aci_punching_eq, en_punching_eq], ids=["aci", "en"])
    def test_the_module_declares_equations(self, module):
        assert _public_equations(module), f"{module.__name__} declares no equations"

    @pytest.mark.parametrize("module", [aci_punching_eq, en_punching_eq], ids=["aci", "en"])
    def test_no_stub_returns_a_number(self, module):
        for function in _public_equations(module):
            n_args = len(inspect.signature(function).parameters)
            with pytest.raises(NotImplementedError, match="punching-roadmap"):
                function(*([0.0] * n_args))

def test_opening_rotation_defaults_to_zero() -> None:
    opening = Opening(
        shape="rectangular",
        x=40 * cm,
        y=0 * cm,
        b=40 * cm,
        h=20 * cm,
    )

    assert opening.rotation.to("degree").magnitude == pytest.approx(0.0)


def test_opening_rotation_is_passed_to_geometry() -> None:
    opening = Opening(
        shape="rectangular",
        x=40 * cm,
        y=0 * cm,
        b=40 * cm,
        h=20 * cm,
        rotation=90 * deg,
    )

    angles = opening.shadow_angles()

    expected = math.atan2(20.0, 30.0)
    assert angles == pytest.approx((-expected, expected))

def test_opening_circular_shadow_converts_units_and_diameter() -> None:
    opening = Opening(
        shape="circular",
        x=400 * mm,
        y=0 * cm,
        diameter=40 * cm,
    )

    angles = opening.shadow_angles()

    # Center distance: 400 mm. Radius: 200 mm.
    # Half-angle: 30 degrees.
    assert angles == pytest.approx((-math.pi / 6, math.pi / 6))

def test_node_interior_section_properties(
    bare_slab: PunchingSlab,
    f1: Forces,
) -> None:
    bare_slab.set_effective_depth(
        d_x=190 * mm,
        d_y=210 * mm,
    )
    column = Column(
        shape="rectangular",
        position="interior",
        b=30 * cm,
        h=500 * mm,
    )
    node = PunchingNode(
        slab=bare_slab,
        column=column,
        forces=f1,
    )

    result = node.interior_section_properties(offset=10 * cm)

    # Critical contour: 500 x 700 mm. Effective depth: 200 mm.
    assert result.b_0 == pytest.approx(2400.0)
    assert result.A_c == pytest.approx(480_000.0)
    assert (result.x_g, result.y_g) == pytest.approx((0.0, 0.0))
    assert result.extents == pytest.approx((500.0, 700.0))

    assert result.J_x == pytest.approx(36_866_666_666.6667)
    assert result.J_y == pytest.approx(22_333_333_333.3333)
    assert result.J_xy == pytest.approx(0.0)

    assert result.parts == 1
    assert len(result.segments) == 4

@pytest.fixture
def interior_geometry_node(
    bare_slab: PunchingSlab, f1: Forces,
) -> PunchingNode:
    bare_slab.set_effective_depth(d_x=190 * mm, d_y=210 * mm)

    column = Column(
        shape="rectangular",
        position="interior",
        b=30 * cm,
        h=50 * cm,
    )

    return PunchingNode(slab=bare_slab, column=column, forces=f1)

@pytest.mark.parametrize("unit", ["mm", "cm", "inch"])
def test_node_properties_normalize_length_units(
    interior_geometry_node: PunchingNode, unit: str,
) -> None:
    node = interior_geometry_node
    node.column.b = (30 * cm).to(unit)
    node.column.h = (50 * cm).to(unit)
    node.slab.set_effective_depth(
        d_x=(190 * mm).to(unit),
        d_y=(210 * mm).to(unit),
    )

    result = node.interior_section_properties((10 * cm).to(unit))

    assert result.b_0 == pytest.approx(2400.0)
    assert result.A_c == pytest.approx(480_000.0)
    assert result.extents == pytest.approx((500.0, 700.0))
    assert result.J_x == pytest.approx(36_866_666_666.6667)
    assert result.J_y == pytest.approx(22_333_333_333.3333)


@pytest.mark.parametrize("quad_segs", [1, 8, 32])
def test_node_properties_pass_rounding_options(
    interior_geometry_node: PunchingNode, quad_segs: int,
) -> None:
    result = interior_geometry_node.interior_section_properties(
        10 * cm,
        corner_style="round",
        quad_segs=quad_segs,
    )

    expected = (
        1600.0
        + 8 * quad_segs * 100.0 * math.sin(math.pi / (4 * quad_segs))
    )

    assert result.b_0 == pytest.approx(expected)
    assert result.A_c == pytest.approx(expected * 200.0)
    assert len(result.segments) == 4 + 4 * quad_segs


def test_node_properties_accept_zero_offset(
    interior_geometry_node: PunchingNode,
) -> None:
    result = interior_geometry_node.interior_section_properties(0 * mm)

    assert result.b_0 == pytest.approx(1600.0)
    assert result.extents == pytest.approx((300.0, 500.0))


@pytest.mark.parametrize("offset", [100.0, None, 1 * MPa])
def test_node_properties_reject_invalid_offset_type(
    interior_geometry_node: PunchingNode, offset: object,
) -> None:
    with pytest.raises(TypeError, match="offset must be a length Quantity"):
        interior_geometry_node.interior_section_properties(
            offset,  # type: ignore[arg-type]
        )


@pytest.mark.parametrize(
    "value",
    [-1.0, float("nan"), float("inf"), float("-inf")],
)
def test_node_properties_forward_invalid_offset_value(
    interior_geometry_node: PunchingNode, value: float,
) -> None:
    with pytest.raises(ValueError, match="Offset must be finite and nonnegative"):
        interior_geometry_node.interior_section_properties(value * cm)


@pytest.mark.parametrize(
    "column",
    [
        Column(
            shape="circular",
            position="interior",
            b=30 * cm,
        ),
        Column(
            shape="rectangular",
            position="edge",
            b=30 * cm,
            h=50 * cm,
            edge_distance_x=20 * cm,
        ),
        Column(
            shape="rectangular",
            position="corner",
            b=30 * cm,
            h=50 * cm,
            edge_distance_x=20 * cm,
            edge_distance_y=30 * cm,
        ),
        Column(
            shape="rectangular",
            position="interior",
            b=30 * cm,
            h=50 * cm,
            edge_distance_x=20 * cm,
        ),
        Column(
            shape="rectangular",
            position="interior",
            b=30 * cm,
            h=50 * cm,
            edge_distance_y=30 * cm,
        ),
    ],
)
def test_node_properties_reject_unsupported_columns(
    interior_geometry_node: PunchingNode, column: Column,
) -> None:
    interior_geometry_node.column = column

    with pytest.raises(NotImplementedError):
        interior_geometry_node.interior_section_properties(10 * cm)


def test_node_properties_reject_openings(
    interior_geometry_node: PunchingNode,
) -> None:
    interior_geometry_node.openings = [
        Opening(
            shape="rectangular",
            x=60 * cm,
            y=0 * cm,
            b=20 * cm,
            h=20 * cm,
        ),
    ]

    with pytest.raises(NotImplementedError, match="Openings"):
        interior_geometry_node.interior_section_properties(10 * cm)


def test_node_properties_reject_capital(
    interior_geometry_node: PunchingNode,
) -> None:
    interior_geometry_node.capital = Capital(
        b=100 * cm,
        h=100 * cm,
        thickness=10 * cm,
    )

    with pytest.raises(NotImplementedError, match="Capitals"):
        interior_geometry_node.interior_section_properties(10 * cm)


def test_node_properties_read_current_dimensions_and_depth(
    interior_geometry_node: PunchingNode,
) -> None:
    node = interior_geometry_node
    first = node.interior_section_properties(10 * cm)

    node.column.b = 40 * cm
    second = node.interior_section_properties(10 * cm)

    assert second.b_0 == pytest.approx(2600.0)
    assert second.A_c == pytest.approx(520_000.0)

    node.slab.set_effective_depth(d_x=210 * mm, d_y=210 * mm)
    third = node.interior_section_properties(10 * cm)

    assert third.b_0 == pytest.approx(2600.0)
    assert third.A_c == pytest.approx(546_000.0)

    # Previously returned results must remain unchanged.
    assert first.b_0 == pytest.approx(2400.0)
    assert first.A_c == pytest.approx(480_000.0)
    assert second.A_c == pytest.approx(520_000.0)


def test_node_geometry_does_not_mark_node_checked(
    interior_geometry_node: PunchingNode,
) -> None:
    node = interior_geometry_node
    node.forces = []  # Geometry does not require a load combination.

    node.interior_section_properties(10 * cm)

    with pytest.raises(PunchingCheckNotRunError):
        _ = node.punching_checks


def test_node_properties_forward_invalid_corner_style(
    interior_geometry_node: PunchingNode,
) -> None:
    with pytest.raises(ValueError, match="Corner style"):
        interior_geometry_node.interior_section_properties(
            10 * cm,
            corner_style="other",  # type: ignore[arg-type]
        )