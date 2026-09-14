from math import pi

import pytest
from shapely.geometry import box

from mento.punching_geometry import (
    PunchingPerimeterCandidate,
    candidate_fits_length_limit,
    candidate_fits_slab,
    clip_perimeter_to_slab,
    column_offset_perimeter,
    corner_perimeter_candidate,
    edge_perimeter_candidate,
    section_properties,
    segments_of,
    rectangular_slab_perimeter_candidates,
)

def test_interior_rectangular_section() -> None:
    # All lengths in mm.
    # Column: 300 x 500; effective depth: 200.
    # Perimeter at d/2 from each face: 500 x 700.
    d = 200.0
    p1 = (-250.0, -350.0)
    p2 = (250.0, -350.0)
    p3 = (250.0, 350.0)
    p4 = (-250.0, 350.0)

    segments = (
        (p1, p2),
        (p2, p3),
        (p3, p4),
        (p4, p1),
    )

    result = section_properties(segments, d)

    # Perimeter [mm], resisting area [mm²] and centroid [mm].
    assert result.b_0 == pytest.approx(2400.0)
    assert result.A_c == pytest.approx(480_000.0)
    assert result.x_g == pytest.approx(0.0)
    assert result.y_g == pytest.approx(0.0)
    assert result.extents == pytest.approx((500.0, 700.0))

    # In-plane contributions [mm⁴].
    assert result.J_x_plan == pytest.approx(35_933_333_333.3333)
    assert result.J_y_plan == pytest.approx(21_666_666_666.6667)

    # Through-depth contributions [mm⁴].
    assert result.dJ_x == pytest.approx(933_333_333.3333)
    assert result.dJ_y == pytest.approx(666_666_666.6667)

    # Total properties [mm⁴].
    assert result.J_x == pytest.approx(36_866_666_666.6667)
    assert result.J_y == pytest.approx(22_333_333_333.3333)
    assert result.J_xy == pytest.approx(0.0)


def test_column_offset_sharp() -> None:
    perimeter = column_offset_perimeter(
        c_x=300.0, c_y=500.0, offset=100.0, corner_style="sharp"
    )
    segments = segments_of(perimeter)
    result = section_properties(segments, d=200.0)

    assert perimeter.is_ring
    assert len(segments) == 4
    assert result.b_0 == pytest.approx(2400.0)
    assert result.A_c == pytest.approx(480_000.0)
    assert (result.x_g, result.y_g) == pytest.approx((0.0, 0.0), abs=1e-10)
    assert result.extents == pytest.approx((500.0, 700.0))
    assert result.J_x == pytest.approx(36_866_666_666.6667)
    assert result.J_y == pytest.approx(22_333_333_333.3333)
    assert result.J_xy == pytest.approx(0.0)

def test_column_offset_round() -> None:
    perimeter = column_offset_perimeter(
        c_x=300.0, c_y=500.0, offset=100.0,
        corner_style="round", quad_segs=32,
    )
    segments = segments_of(perimeter)
    result = section_properties(segments, d=200.0)

    # Four straight sides plus four quarter circles of radius 100 mm.
    exact_b_0 = 2 * (300.0 + 500.0) + 2 * pi * 100.0

    assert perimeter.is_ring
    assert len(segments) == 4 + 4 * 32
    assert result.b_0 == pytest.approx(exact_b_0, rel=1e-4)
    assert result.A_c == pytest.approx(exact_b_0 * 200.0, rel=1e-4)
    assert (result.x_g, result.y_g) == pytest.approx((0.0, 0.0), abs=1e-10)
    assert result.extents == pytest.approx((500.0, 700.0))

def test_round_offset_convergence() -> None:
    coarse = column_offset_perimeter(
        c_x=300.0, c_y=500.0, offset=100.0,
        corner_style="round", quad_segs=32,
    )
    fine = column_offset_perimeter(
        c_x=300.0, c_y=500.0, offset=100.0,
        corner_style="round", quad_segs=64,
    )
    exact_b_0 = 2 * (300.0 + 500.0) + 2 * pi * 100.0

    assert coarse.length < fine.length < exact_b_0
    assert fine.length == pytest.approx(exact_b_0, rel=1e-5)

def test_clip_perimeter_at_right_edge() -> None:
    # All lengths in mm. Column: 300 x 500, centered at the origin.
    perimeter = column_offset_perimeter(
        c_x=300.0, c_y=500.0, offset=100.0, corner_style="sharp"
    )

    # Right free edge at x = 150: coincident with the column face.
    # The other slab edges are well beyond the critical contour.
    slab_outline = box(-2000.0, -2000.0, 150.0, 2000.0)

    clipped = clip_perimeter_to_slab(perimeter, slab_outline)
    segments = segments_of(clipped)
    result = section_properties(segments, d=200.0)

    assert not clipped.is_closed
    assert len(segments) == 3
    assert slab_outline.covers(clipped)
    assert clipped.intersection(slab_outline.boundary).length == 0.0

    # One vertical side of 700 mm and two horizontal sides of 400 mm.
    assert result.b_0 == pytest.approx(1500.0)
    assert result.A_c == pytest.approx(300_000.0)
    assert result.extents == pytest.approx((400.0, 700.0))
    assert result.x_g == pytest.approx(-143.333333333333)
    assert result.y_g == pytest.approx(0.0, abs=1e-10)

def test_clip_perimeter_at_corner() -> None:
    # All lengths in mm. Column: 300 x 500, centered at the origin.
    perimeter = column_offset_perimeter(
        c_x=300.0, c_y=500.0, offset=100.0, corner_style="sharp"
    )

    # Free edges at x = 150 and y = 250, coincident with column faces.
    slab_outline = box(-2000.0, -2000.0, 150.0, 250.0)

    clipped = clip_perimeter_to_slab(perimeter, slab_outline)
    segments = segments_of(clipped)
    result = section_properties(segments, d=200.0)

    assert not clipped.is_closed
    assert len(segments) == 2
    assert slab_outline.covers(clipped)
    assert clipped.intersection(slab_outline.boundary).length == 0.0

    # One vertical side of 600 mm and one horizontal side of 400 mm.
    assert result.b_0 == pytest.approx(1000.0)
    assert result.A_c == pytest.approx(200_000.0)
    assert result.extents == pytest.approx((400.0, 600.0))
    assert result.x_g == pytest.approx(-170.0)
    assert result.y_g == pytest.approx(-170.0)
    assert result.J_xy == pytest.approx(-2_880_000_000.0)

def test_length_filter_preserves_corner_combination() -> None:
    # All lengths in cm. Closed candidate: 192 cm.
    max_length = column_offset_perimeter(30.0, 50.0, 4.0).length

    # Right edge alone: 246 cm, so discard this candidate.
    assert not candidate_fits_length_limit(
        30.0, 50.0, 4.0,
        max_length=max_length,
        face_gap_x=60.0,
    )

    # The same right edge with the top edge: 153 cm, so retain it.
    assert candidate_fits_length_limit(
        30.0, 50.0, 4.0,
        max_length=max_length,
        face_gap_x=60.0,
        face_gap_y=5.0,
    )

    # A more distant corner: 198 cm, so discard it.
    assert not candidate_fits_length_limit(
        30.0, 50.0, 4.0,
        max_length=max_length,
        face_gap_x=100.0,
        face_gap_y=10.0,
    )

def test_edge_candidate_extends_beyond_offset() -> None:
    # All lengths in cm.
    # Column: 30 x 50. Right face: x = 15.
    # Free edge: x = 25, therefore 10 cm from the column face.
    edge = ((25.0, -100.0), (25.0, 100.0))

    candidate = edge_perimeter_candidate(
        c_x=30.0,
        c_y=50.0,
        offset=4.0,
        edge=edge,
        corner_style="sharp",
    )

    assert candidate.reached_edges == (edge,)
    assert not candidate.perimeter.is_closed
    assert len(segments_of(candidate.perimeter)) == 3

    # Two horizontal legs of 44 cm and one vertical side of 58 cm.
    assert candidate.perimeter.length == pytest.approx(146.0)
    assert candidate.perimeter.bounds == pytest.approx(
        (-19.0, -29.0, 25.0, 29.0)
    )

def test_corner_candidate_extends_to_two_edges() -> None:
    # All lengths in cm. Column: 30 x 50, centered at the origin.
    # Right edge: x = 75, therefore 60 cm from the column face.
    # Top edge: y = 30, therefore 5 cm from the column face.
    vertical_edge = ((75.0, -100.0), (75.0, 30.0))
    horizontal_edge = ((-100.0, 30.0), (75.0, 30.0))

    candidate = corner_perimeter_candidate(
        c_x=30.0,
        c_y=50.0,
        offset=4.0,
        vertical_edge=vertical_edge,
        horizontal_edge=horizontal_edge,
    )

    assert candidate.reached_edges == (vertical_edge, horizontal_edge)
    assert not candidate.perimeter.is_closed
    assert len(segments_of(candidate.perimeter)) == 2

    # Horizontal leg: 94 cm. Vertical leg: 59 cm.
    assert candidate.perimeter.length == pytest.approx(153.0)
    assert candidate.perimeter.bounds == pytest.approx(
        (-19.0, -29.0, 75.0, 30.0)
    )

    slab_outline = box(-100.0, -100.0, 75.0, 30.0)
    assert slab_outline.covers(candidate.perimeter)

    closed = column_offset_perimeter(30.0, 50.0, 4.0)
    assert candidate.perimeter.length < closed.length

def test_candidate_fits_slab_at_corner() -> None:
    # All lengths in cm.
    # The right and top slab edges coincide with the column faces.
    slab_outline = box(-100.0, -100.0, 15.0, 25.0)

    closed = PunchingPerimeterCandidate(
        perimeter=column_offset_perimeter(30.0, 50.0, 4.0)
    )

    corner = corner_perimeter_candidate(
        c_x=30.0,
        c_y=50.0,
        offset=4.0,
        vertical_edge=((15.0, -100.0), (15.0, 25.0)),
        horizontal_edge=((-100.0, 25.0), (15.0, 25.0)),
    )

    assert not candidate_fits_slab(closed, slab_outline)
    assert candidate_fits_slab(corner, slab_outline)

def test_rectangular_slab_perimeter_candidates() -> None:
    # All lengths in cm.
    slab_outline = box(-100.0, -100.0, 75.0, 30.0)

    # Without a length limit, retain all nine fitting candidates.
    candidates = rectangular_slab_perimeter_candidates(
        30.0, 50.0, 4.0, slab_outline
    )

    assert len(candidates) == 9
    assert sorted(c.perimeter.length for c in candidates) == pytest.approx(
        [153.0, 156.0, 178.0, 192.0, 223.0, 246.0, 248.0, 296.0, 296.0]
    )

    # Apply a length limit only when explicitly requested.
    limited = rectangular_slab_perimeter_candidates(
        30.0, 50.0, 4.0, slab_outline, max_length=192.0
    )

    assert len(limited) == 4
    assert sorted(c.perimeter.length for c in limited) == pytest.approx(
        [153.0, 156.0, 178.0, 192.0]
    )