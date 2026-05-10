import pytest
import matplotlib

matplotlib.use("Agg")

from typing import Generator

from mento.column import Column
from mento.punching import Capital, Opening, PunchingNode, PunchingSlab
from mento.forces import Forces
from mento.material import Concrete_ACI_318_19, Concrete_EN_1992_2004, SteelBar
from mento.units import cm, mm, kN, kNm, MPa, inch, psi


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
def steel() -> SteelBar:
    return SteelBar(name="ADN 420", f_y=420 * MPa)


@pytest.fixture
def slab(conc_aci, steel) -> PunchingSlab:
    return PunchingSlab(
        concrete=conc_aci,
        steel_bar=steel,
        h=25 * cm,
        c_c=25 * mm,
        rho_x=0.010,
        rho_y=0.010,
    )


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
    def test_d_avg_metric_estimate(self, slab):
        # d_avg = h - c_c - 16mm = 250mm - 25mm - 16mm = 209mm = 20.9cm
        assert slab.d_avg.to("mm").magnitude == pytest.approx(209.0)

    def test_d_avg_override(self, slab):
        slab.d_avg = 19.8 * cm
        assert slab.d_avg.to("cm").magnitude == pytest.approx(19.8)

    def test_unit_system_metric(self, slab):
        assert slab.unit_system == "metric"

    def test_unit_system_from_concrete(self, conc_en, steel):
        slab_en = PunchingSlab(concrete=conc_en, steel_bar=steel, h=25 * cm, c_c=25 * mm)
        assert slab_en.unit_system == "metric"

    def test_rho_defaults_to_zero(self, conc_aci, steel):
        slab = PunchingSlab(concrete=conc_aci, steel_bar=steel, h=25 * cm, c_c=25 * mm)
        assert slab.rho_x == pytest.approx(0.0)
        assert slab.rho_y == pytest.approx(0.0)

    def test_rho_stored(self, slab):
        assert slab.rho_x == pytest.approx(0.010)
        assert slab.rho_y == pytest.approx(0.010)

    def test_h_and_cc_stored(self, slab):
        assert slab.h.to("cm").magnitude == pytest.approx(25)
        assert slab.c_c.to("mm").magnitude == pytest.approx(25)

    def test_d_avg_thicker_slab(self, conc_aci, steel):
        slab = PunchingSlab(concrete=conc_aci, steel_bar=steel, h=30 * cm, c_c=30 * mm)
        # 300mm - 30mm - 16mm = 254mm
        assert slab.d_avg.to("mm").magnitude == pytest.approx(254.0)

    def test_d_avg_imperial_estimate(self, steel):
        conc_imp = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
        slab = PunchingSlab(concrete=conc_imp, steel_bar=steel, h=10 * inch, c_c=1 * inch)
        # d_avg = 10 - 1 - 0.625 = 8.375 inch
        assert slab.d_avg.to("inch").magnitude == pytest.approx(8.375)


# ---------------------------------------------------------------------------
# Opening
# ---------------------------------------------------------------------------


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
