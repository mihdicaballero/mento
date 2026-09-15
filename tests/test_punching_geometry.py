from math import atan2, pi

import pytest
from shapely.geometry import LineString, MultiLineString, Polygon, box
from typing import Literal

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
    rectangular_opening_shadow_angles,
    subtract_opening_shadow,
    circular_opening_shadow_angles,
    Segment,
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

def test_rectangular_opening_shadow_rotated_90_degrees() -> None:
    # Opening centered at (40, 0), with local dimensions 40 x 20.
    # After rotating 90 degrees, its global dimensions are 20 x 40.
    # Therefore, its nearest corners are (30, -20) and (30, 20).
    angles = rectangular_opening_shadow_angles(
        x=40.0,
        y=0.0,
        b=40.0,
        h=20.0,
        rotation=pi / 2,
    )

    expected = atan2(20.0, 30.0)
    assert angles == pytest.approx((-expected, expected))

def test_subtract_opening_shadow_from_interior_perimeter() -> None:
    # All lengths in cm.
    # Critical rectangle: 38 x 58, total perimeter 192.
    original = column_offset_perimeter(30.0, 50.0, 4.0)

    start, end = rectangular_opening_shadow_angles(
        x=60.0,
        y=0.0,
        b=20.0,
        h=20.0,
    )

    reduced = subtract_opening_shadow(original, start, end)

    # Tangents have slopes +/-10/50.
    # On the right side, x=19, they intersect at y=+/-3.8.
    # Removed length: 7.6 cm.
    assert reduced.length == pytest.approx(192.0 - 7.6)
    assert len(segments_of(reduced)) == 5
    assert isinstance(reduced, LineString)
    assert not reduced.is_closed

    # The original perimeter remains unchanged.
    assert original.length == pytest.approx(192.0)
    assert original.is_closed

def test_overlapping_opening_shadows_are_not_counted_twice() -> None:
    # All lengths in cm. Original perimeter: 192.
    original = column_offset_perimeter(30.0, 50.0, 4.0)

    # Opening 1: x=[50, 70], y=[-10, 10].
    start_1, end_1 = rectangular_opening_shadow_angles(
        x=60.0, y=0.0, b=20.0, h=20.0,
    )

    # Opening 2: x=[100, 140], y=[0, 40].
    # The openings are separate, but their angular shadows overlap.
    start_2, end_2 = rectangular_opening_shadow_angles(
        x=120.0, y=20.0, b=40.0, h=40.0,
    )

    first_only = subtract_opening_shadow(original, start_1, end_1)
    second_only = subtract_opening_shadow(original, start_2, end_2)
    both = subtract_opening_shadow(first_only, start_2, end_2)

    # On x=19, shadow 1 covers y=[-3.8, 3.8], shadow 2 y=[0, 7.6].
    # Each removes 7.6, but they share 3.8. Their union removes 11.4.
    assert first_only.length == pytest.approx(184.4)
    assert second_only.length == pytest.approx(184.4)
    assert both.length == pytest.approx(180.6)
    assert first_only.length - both.length == pytest.approx(3.8)

    assert isinstance(both, LineString)
    assert not both.is_closed
    assert len(segments_of(both)) == 5
    endpoints_y = sorted((both.coords[0][1], both.coords[-1][1]))
    assert endpoints_y == pytest.approx([-3.8, 7.6])

    # Reversing the order must retain the same length.
    reversed_order = subtract_opening_shadow(second_only, start_1, end_1)
    assert reversed_order.length == pytest.approx(both.length)

    # The original perimeter remains unchanged.
    assert original.length == pytest.approx(192.0)
    assert original.is_closed

def test_circular_opening_shadow_angles() -> None:
    # Distance to the opening center: 40. Radius: 20.
    # sin(half_angle) = 20 / 40 -> half_angle = 30 degrees.
    start, end = circular_opening_shadow_angles(
        x=40.0,
        y=0.0,
        radius=20.0,
    )

    assert start == pytest.approx(-pi / 6)
    assert end == pytest.approx(pi / 6)


def test_subtract_circular_opening_shadow() -> None:
    # All lengths in cm.
    original = column_offset_perimeter(30.0, 50.0, 4.0)

    start, end = circular_opening_shadow_angles(
        x=50.0,
        y=0.0,
        radius=30.0,
    )
    reduced = subtract_opening_shadow(original, start, end)

    # Tangent slopes are +/-3/4.
    # At x=19, they intersect the perimeter at y=+/-14.25.
    # Removed length: 28.5 cm.
    assert reduced.length == pytest.approx(163.5)
    assert isinstance(reduced, LineString)
    assert not reduced.is_closed
    assert sorted(
        (reduced.coords[0][1], reduced.coords[-1][1])
    ) == pytest.approx([-14.25, 14.25])

    # The original perimeter remains unchanged.
    assert original.length == pytest.approx(192.0)
    assert original.is_closed

def test_section_properties_arms_from_shifted_centroid() -> None:
    # Rectangle 20 x 40, with centroid at (20, 40), not at the origin.
    p1 = (10.0, 20.0)
    p2 = (30.0, 20.0)
    p3 = (30.0, 60.0)
    p4 = (10.0, 60.0)
    segments = ((p1, p2), (p2, p3), (p3, p4), (p4, p1))

    result = section_properties(segments, d=5.0)

    assert (result.x_g, result.y_g) == pytest.approx((20.0, 40.0))
    assert result.b_0 == pytest.approx(120.0)
    assert result.A_c == pytest.approx(600.0)
    assert result.extents == pytest.approx((20.0, 40.0))
    assert result.arms((25.0, 30.0)) == pytest.approx((5.0, -10.0))
    assert result.arms((15.0, 50.0)) == pytest.approx((-5.0, 10.0))
    assert result.arms((20.0, 40.0)) == pytest.approx((0.0, 0.0))


def test_section_properties_rejects_empty_segments() -> None:
    with pytest.raises(ValueError, match="at least one segment"):
        section_properties([], d=5.0)


@pytest.mark.parametrize("d", [0.0, -1.0, float("nan"), float("inf"), float("-inf")])
def test_section_properties_rejects_invalid_depth(d: float) -> None:
    segments = [((0.0, 0.0), (10.0, 0.0))]

    with pytest.raises(ValueError, match="Effective depth"):
        section_properties(segments, d=d)


@pytest.mark.parametrize("coordinate_index", range(4))
@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_section_properties_rejects_nonfinite_coordinates(
    coordinate_index: int,
    value: float,
) -> None:
    # Test each coordinate: x1, y1, x2, y2.
    coords = [0.0, 0.0, 10.0, 0.0]
    coords[coordinate_index] = value
    segments = [((coords[0], coords[1]), (coords[2], coords[3]))]

    with pytest.raises(ValueError, match="coordinates must be finite"):
        section_properties(segments, d=5.0)


def test_section_properties_rejects_zero_total_length() -> None:
    segments = [((10.0, 20.0), (10.0, 20.0))]

    with pytest.raises(ValueError, match="positive length"):
        section_properties(segments, d=5.0)

@pytest.mark.parametrize("axis", ["x", "y"])
@pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), float("inf"), float("-inf")])
def test_column_offset_rejects_invalid_dimensions(axis: str, value: float) -> None:
    dimensions = {"x": 30.0, "y": 50.0}
    dimensions[axis] = value

    with pytest.raises(ValueError, match="Column dimensions"):
        column_offset_perimeter(dimensions["x"], dimensions["y"], 4.0)


@pytest.mark.parametrize("offset", [-1.0, float("nan"), float("inf"), float("-inf")])
def test_column_offset_rejects_invalid_offset(offset: float) -> None:
    with pytest.raises(ValueError, match="Offset must be"):
        column_offset_perimeter(30.0, 50.0, offset)


def test_column_offset_rejects_unknown_corner_style() -> None:
    with pytest.raises(ValueError, match="Corner style"):
        column_offset_perimeter(30.0, 50.0, 4.0, corner_style="other")  # type: ignore[arg-type]


@pytest.mark.parametrize("quad_segs", [0, -1, 1.5, True, False])
def test_column_offset_rejects_invalid_quad_segs(quad_segs: int | float) -> None:
    with pytest.raises(ValueError, match="quad_segs"):
        column_offset_perimeter(30.0, 50.0, 4.0, quad_segs=quad_segs)  # type: ignore[arg-type]


def test_column_offset_zero_preserves_column_boundary() -> None:
    perimeter = column_offset_perimeter(30.0, 50.0, 0.0)

    assert perimeter.equals(box(-15.0, -25.0, 15.0, 25.0).boundary)
    assert perimeter.is_ring
    assert perimeter.length == pytest.approx(160.0)
    assert len(segments_of(perimeter)) == 4

def test_segments_of_does_not_join_or_close_separate_lines() -> None:
    perimeter = MultiLineString([
        [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)],
        [(20.0, 0.0), (30.0, 0.0)],
    ])

    assert segments_of(perimeter) == [
        ((0.0, 0.0), (10.0, 0.0)),
        ((10.0, 0.0), (10.0, 10.0)),
        ((20.0, 0.0), (30.0, 0.0)),
    ]


def test_segments_of_skips_repeated_points() -> None:
    perimeter = LineString([
        (0.0, 0.0), (0.0, 0.0),
        (10.0, 0.0), (10.0, 0.0), (10.0, 10.0),
    ])

    assert segments_of(perimeter) == [
        ((0.0, 0.0), (10.0, 0.0)),
        ((10.0, 0.0), (10.0, 10.0)),
    ]


@pytest.mark.parametrize("perimeter", [LineString(), MultiLineString([])])
def test_segments_of_empty_perimeter(perimeter: LineString | MultiLineString) -> None:
    assert segments_of(perimeter) == []


def test_segments_of_rejects_polygon() -> None:
    with pytest.raises(TypeError, match="LineString or MultiLineString"):
        segments_of(box(0.0, 0.0, 10.0, 10.0))  # type: ignore[arg-type]


@pytest.mark.parametrize("perimeter", [
    LineString([(0.0, 0.0, 0.0), (10.0, 0.0, 1.0)]),
    MultiLineString([[(0.0, 0.0, 0.0), (10.0, 0.0, 1.0)]]),
])
def test_segments_of_rejects_3d_coordinates(perimeter: LineString | MultiLineString) -> None:
    with pytest.raises(ValueError, match="two-dimensional"):
        segments_of(perimeter)

def test_clip_preserves_perimeter_fully_inside_slab() -> None:
    perimeter = column_offset_perimeter(30.0, 50.0, 4.0)
    slab_outline = box(-100.0, -100.0, 100.0, 100.0)

    clipped = clip_perimeter_to_slab(perimeter, slab_outline)

    assert clipped.equals(perimeter)
    assert isinstance(clipped, LineString)
    assert clipped.is_ring
    assert clipped.length == pytest.approx(192.0)


@pytest.mark.parametrize("perimeter", [
    LineString([(-5.0, 2.0), (15.0, 2.0), (15.0, 8.0), (-5.0, 8.0)]),
    MultiLineString([
        [(-5.0, 2.0), (15.0, 2.0)],
        [(-5.0, 8.0), (15.0, 8.0)],
    ]),
])
def test_clip_keeps_disconnected_parts(perimeter: LineString | MultiLineString) -> None:
    slab_outline = box(0.0, 0.0, 10.0, 10.0)
    original_wkt = perimeter.wkt
    expected = MultiLineString([
        [(0.0, 2.0), (10.0, 2.0)],
        [(0.0, 8.0), (10.0, 8.0)],
    ])

    clipped = clip_perimeter_to_slab(perimeter, slab_outline)

    assert isinstance(clipped, MultiLineString)
    assert len(clipped.geoms) == 2
    assert clipped.equals(expected)
    assert clipped.length == pytest.approx(20.0)
    assert clipped.intersection(slab_outline.boundary).length == 0.0
    assert perimeter.wkt == original_wkt


def test_clip_removes_segments_on_free_edges() -> None:
    perimeter = LineString([(0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (5.0, 5.0)])
    slab_outline = box(0.0, 0.0, 10.0, 10.0)

    clipped = clip_perimeter_to_slab(perimeter, slab_outline)

    assert clipped.equals(LineString([(10.0, 5.0), (5.0, 5.0)]))
    assert clipped.length == pytest.approx(5.0)


@pytest.mark.parametrize("perimeter", [
    LineString(),
    MultiLineString([]),
    LineString([(11.0, 11.0), (12.0, 12.0)]),  # Outside the slab.
    LineString([(0.0, 0.0), (10.0, 0.0)]),  # On a free edge.
    LineString([(-1.0, 1.0), (1.0, -1.0)]),  # Touches only at (0, 0).
])
def test_clip_returns_empty_line_when_nothing_remains(perimeter: LineString | MultiLineString) -> None:
    clipped = clip_perimeter_to_slab(perimeter, box(0.0, 0.0, 10.0, 10.0))

    assert isinstance(clipped, LineString)
    assert clipped.is_empty
    assert clipped.length == 0.0

def test_clip_rejects_polygon_as_perimeter() -> None:
    with pytest.raises(TypeError, match="Perimeter must be"):
        clip_perimeter_to_slab(box(1.0, 1.0, 2.0, 2.0), box(0.0, 0.0, 10.0, 10.0))  # type: ignore[arg-type]


def test_clip_rejects_line_as_slab() -> None:
    perimeter = LineString([(1.0, 1.0), (9.0, 1.0)])

    with pytest.raises(TypeError, match="Slab outline must be a Polygon"):
        clip_perimeter_to_slab(perimeter, perimeter)  # type: ignore[arg-type]


@pytest.mark.parametrize("slab_outline", [
    Polygon(),
    Polygon([(0.0, 0.0), (10.0, 10.0), (0.0, 10.0), (10.0, 0.0)]),
])
def test_clip_rejects_empty_or_invalid_slab(slab_outline: Polygon) -> None:
    perimeter = LineString([(1.0, 1.0), (9.0, 1.0)])

    with pytest.raises(ValueError, match="nonempty valid polygon"):
        clip_perimeter_to_slab(perimeter, slab_outline)


def test_clip_rejects_slab_with_holes() -> None:
    slab_outline = Polygon(
        shell=[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
        holes=[[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0)]],
    )
    perimeter = LineString([(1.0, 1.0), (9.0, 1.0)])
    assert slab_outline.is_valid

    with pytest.raises(ValueError, match="Define openings separately"):
        clip_perimeter_to_slab(perimeter, slab_outline)


@pytest.mark.parametrize("perimeter,slab_outline", [
    (LineString([(1.0, 1.0, 0.0), (9.0, 1.0, 0.0)]), box(0.0, 0.0, 10.0, 10.0)),
    (
        LineString([(1.0, 1.0), (9.0, 1.0)]),
        Polygon([(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (10.0, 10.0, 0.0), (0.0, 10.0, 0.0)]),
    ),
])
def test_clip_rejects_3d_geometry(perimeter: LineString, slab_outline: Polygon) -> None:
    with pytest.raises(ValueError, match="two-dimensional"):
        clip_perimeter_to_slab(perimeter, slab_outline)


def test_clip_rejects_invalid_perimeter() -> None:
    perimeter = LineString([(1.0, 1.0), (1.0, 1.0)])
    assert not perimeter.is_valid

    with pytest.raises(ValueError, match="Perimeter geometry must be valid"):
        clip_perimeter_to_slab(perimeter, box(0.0, 0.0, 10.0, 10.0))


def test_length_filter_has_no_implicit_limit() -> None:
    # Much longer than the interior perimeter, but no limit was requested.
    assert candidate_fits_length_limit(30.0, 50.0, 4.0, face_gap_x=1000.0)


@pytest.mark.parametrize("style,gap_x,gap_y,length", [
    ("sharp", None, None, 192.0),
    ("sharp", 10.0, None, 146.0),
    ("sharp", None, 5.0, 156.0),
    ("sharp", 10.0, 5.0, 103.0),
    ("sharp", 0.0, None, 126.0),
    ("sharp", None, 0.0, 146.0),
    ("sharp", 0.0, 0.0, 88.0),
    # With quad_segs=1, each rounded quarter is one chord: 4 * sqrt(2).
    ("round", None, None, 160.0 + 16.0 * 2.0**0.5),
    ("round", 10.0, None, 130.0 + 8.0 * 2.0**0.5),
    ("round", None, 5.0, 140.0 + 8.0 * 2.0**0.5),
    ("round", 10.0, 5.0, 95.0 + 4.0 * 2.0**0.5),
])
def test_length_filter_known_lengths(
    style: Literal["sharp", "round"],
    gap_x: float | None,
    gap_y: float | None,
    length: float,
) -> None:
    assert candidate_fits_length_limit(
        30.0, 50.0, 4.0, max_length=length,
        face_gap_x=gap_x, face_gap_y=gap_y, corner_style=style, quad_segs=1,
    )
    assert not candidate_fits_length_limit(
        30.0, 50.0, 4.0, max_length=length - 1.0,
        face_gap_x=gap_x, face_gap_y=gap_y, corner_style=style, quad_segs=1,
    )


@pytest.mark.parametrize("quad_segs", [1, 8, 32])
def test_length_filter_matches_generated_round_perimeter(quad_segs: int) -> None:
    perimeter = column_offset_perimeter(30.0, 50.0, 4.0, corner_style="round", quad_segs=quad_segs)

    assert candidate_fits_length_limit(
        30.0, 50.0, 4.0, max_length=perimeter.length,
        corner_style="round", quad_segs=quad_segs,
    )
    assert not candidate_fits_length_limit(
        30.0, 50.0, 4.0, max_length=perimeter.length - 1.0,
        corner_style="round", quad_segs=quad_segs,
    )


def test_length_filter_allows_only_tiny_roundoff() -> None:
    assert candidate_fits_length_limit(30.0, 50.0, 4.0, max_length=192.0 - 1e-11)
    assert not candidate_fits_length_limit(30.0, 50.0, 4.0, max_length=192.0 - 1e-6)


def test_length_filter_accepts_zero_offset_and_zero_limit() -> None:
    assert candidate_fits_length_limit(30.0, 50.0, 0.0, max_length=160.0)
    assert not candidate_fits_length_limit(30.0, 50.0, 0.0, max_length=159.0)
    # A zero limit is valid, but this nonzero perimeter does not fit it.
    assert not candidate_fits_length_limit(30.0, 50.0, 0.0, max_length=0.0)

@pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), float("inf"), float("-inf")])
def test_length_filter_rejects_invalid_dimensions(value: float) -> None:
    with pytest.raises(ValueError, match="Column dimensions"):
        candidate_fits_length_limit(value, 50.0, 4.0)
    with pytest.raises(ValueError, match="Column dimensions"):
        candidate_fits_length_limit(30.0, value, 4.0)


@pytest.mark.parametrize("offset", [-1.0, float("nan"), float("inf"), float("-inf")])
def test_length_filter_rejects_invalid_offset(offset: float) -> None:
    with pytest.raises(ValueError, match="Offset must be"):
        candidate_fits_length_limit(30.0, 50.0, offset)


@pytest.mark.parametrize("gap", [-1.0, float("nan"), float("inf"), float("-inf")])
def test_length_filter_rejects_invalid_gaps(gap: float) -> None:
    with pytest.raises(ValueError, match="Face gaps"):
        candidate_fits_length_limit(30.0, 50.0, 4.0, face_gap_x=gap)
    with pytest.raises(ValueError, match="Face gaps"):
        candidate_fits_length_limit(30.0, 50.0, 4.0, face_gap_y=gap)


@pytest.mark.parametrize("limit", [-1.0, float("nan"), float("inf"), float("-inf")])
def test_length_filter_rejects_invalid_limit(limit: float) -> None:
    with pytest.raises(ValueError, match="Maximum length"):
        candidate_fits_length_limit(30.0, 50.0, 4.0, max_length=limit)


def test_length_filter_rejects_unknown_corner_style() -> None:
    with pytest.raises(ValueError, match="Corner style"):
        candidate_fits_length_limit(30.0, 50.0, 4.0, corner_style="other")  # type: ignore[arg-type]


@pytest.mark.parametrize("quad_segs", [0, -1, 1.5, True, False])
def test_length_filter_rejects_invalid_quad_segs(quad_segs: int | float) -> None:
    with pytest.raises(ValueError, match="quad_segs"):
        candidate_fits_length_limit(30.0, 50.0, 4.0, quad_segs=quad_segs)  # type: ignore[arg-type]


def test_length_filter_rejects_numeric_overflow() -> None:
    # Inputs are finite, but their resulting perimeter exceeds float range.
    with pytest.raises(ValueError, match="finite numeric range"):
        candidate_fits_length_limit(1e308, 1e308, 4.0)

@pytest.mark.parametrize("style", ["sharp", "round"])
@pytest.mark.parametrize("edge,bounds,straight_length", [
    (((25.0, -29.0), (25.0, 29.0)), (-19.0, -29.0, 25.0, 29.0), 130.0),
    (((-25.0, -29.0), (-25.0, 29.0)), (-25.0, -29.0, 19.0, 29.0), 130.0),
    (((-19.0, 35.0), (19.0, 35.0)), (-19.0, -29.0, 19.0, 35.0), 150.0),
    (((-19.0, -35.0), (19.0, -35.0)), (-19.0, -35.0, 19.0, 29.0), 150.0),
])
def test_edge_candidate_all_four_sides(
    style: Literal["sharp", "round"],
    edge: Segment,
    bounds: tuple[float, float, float, float],
    straight_length: float,
) -> None:
    candidate = edge_perimeter_candidate(30.0, 50.0, 4.0, edge, corner_style=style, quad_segs=1)
    perimeter = candidate.perimeter

    assert isinstance(perimeter, LineString)
    assert perimeter.is_valid and perimeter.is_simple
    assert not perimeter.is_closed
    assert perimeter.bounds == pytest.approx(bounds)
    assert candidate.reached_edges == (edge,)
    assert perimeter.boundary.equals(LineString(edge).boundary)
    assert perimeter.intersection(LineString(edge)).length == 0.0

    # Two retained corners: sharp or one diagonal chord per rounded quarter.
    corners = 16.0 if style == "sharp" else 8.0 * 2.0**0.5
    assert perimeter.length == pytest.approx(straight_length + corners)
    assert len(segments_of(perimeter)) == (3 if style == "sharp" else 5)

    reversed_edge = (edge[1], edge[0])
    reversed_candidate = edge_perimeter_candidate(
        30.0, 50.0, 4.0, reversed_edge, corner_style=style, quad_segs=1,
    )
    assert reversed_candidate.perimeter.equals(perimeter)
    assert reversed_candidate.reached_edges == (reversed_edge,)


@pytest.mark.parametrize("offset", [0.0, 4.0])
@pytest.mark.parametrize("gap", [0.0, 2.0, 4.0, 10.0])
def test_edge_candidate_reaches_near_and_far_edges(offset: float, gap: float) -> None:
    edge_x = 15.0 + gap
    reach_y = 25.0 + offset
    edge = ((edge_x, -reach_y), (edge_x, reach_y))

    candidate = edge_perimeter_candidate(30.0, 50.0, offset, edge)
    expected = LineString([
        (edge_x, -reach_y), (-15.0 - offset, -reach_y),
        (-15.0 - offset, reach_y), (edge_x, reach_y),
    ])

    assert candidate.perimeter.equals(expected)
    assert not candidate.perimeter.is_closed


def test_edge_candidate_round_with_fine_discretization() -> None:
    edge = ((25.0, -29.0), (25.0, 29.0))
    candidate = edge_perimeter_candidate(30.0, 50.0, 4.0, edge, corner_style="round", quad_segs=32)

    # Three straight sections and two quarter circles, each with 32 chords.
    exact_arc_length = 130.0 + 4.0 * pi
    assert candidate.perimeter.length < exact_arc_length
    assert candidate.perimeter.length == pytest.approx(exact_arc_length, rel=1e-5)
    assert len(segments_of(candidate.perimeter)) == 67
    assert candidate.perimeter.boundary.equals(LineString(edge).boundary)

@pytest.mark.parametrize("edge,message", [
    ((), "two 2D endpoints"),
    (((25.0, -29.0),), "two 2D endpoints"),
    (((25.0, -29.0), (25.0, 0.0), (25.0, 29.0)), "two 2D endpoints"),
    (((25.0, -29.0, 0.0), (25.0, 29.0, 0.0)), "two 2D endpoints"),
    (((float("nan"), -29.0), (25.0, 29.0)), "coordinates must be finite"),
    (((25.0, float("inf")), (25.0, 29.0)), "coordinates must be finite"),
    (((25.0, -29.0), (25.0, float("-inf"))), "coordinates must be finite"),
    (((25.0, 0.0), (25.0, 0.0)), "nonzero and axis-aligned"),
    (((25.0, -29.0), (26.0, 29.0)), "nonzero and axis-aligned"),
    (((10.0, -29.0), (10.0, 29.0)), "column interior"),
    (((-10.0, -29.0), (-10.0, 29.0)), "column interior"),
    (((-19.0, 20.0), (19.0, 20.0)), "column interior"),
    (((-19.0, -20.0), (19.0, -20.0)), "column interior"),
    (((25.0, -28.0), (25.0, 29.0)), "contain both contour endpoints"),
    (((25.0, -29.0), (25.0, 28.0)), "contain both contour endpoints"),
])
def test_edge_candidate_rejects_invalid_edge(
    edge: tuple[tuple[float, ...], ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        edge_perimeter_candidate(30.0, 50.0, 4.0, edge)  # type: ignore[arg-type]


@pytest.mark.filterwarnings("ignore:overflow encountered in buffer:RuntimeWarning")
def test_edge_candidate_rejects_numeric_overflow() -> None:
    # Finite inputs, but the temporary clipping box exceeds float range.
    edge = ((6.5e307, -1.0), (6.5e307, 1.0))
    with pytest.raises(ValueError, match="Geometry exceeds the finite numeric range"):
        edge_perimeter_candidate(1.3e308, 1.0, 0.0, edge)

@pytest.mark.parametrize("opened", [
    LineString(),
    LineString([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]),
    MultiLineString([[(0.0, 0.0), (1.0, 0.0)], [(2.0, 0.0), (3.0, 0.0)]]),
])
def test_edge_candidate_rejects_unexpected_clip_result(
    monkeypatch: pytest.MonkeyPatch,
    opened: LineString | MultiLineString,
) -> None:
    # Simulate the helper returning an empty, closed or disconnected contour.
    monkeypatch.setattr("mento.punching_geometry.clip_perimeter_to_slab", lambda *_args: opened)

    with pytest.raises(ValueError, match="single open contour"):
        edge_perimeter_candidate(30.0, 50.0, 4.0, ((25.0, -29.0), (25.0, 29.0)))

@pytest.mark.parametrize("sign_x,sign_y", [(1, 1), (-1, 1), (-1, -1), (1, -1)])
@pytest.mark.parametrize("style", ["sharp", "round"])
def test_corner_candidate_all_four_corners(
    sign_x: int,
    sign_y: int,
    style: Literal["sharp", "round"],
) -> None:
    vertical = ((25.0 * sign_x, -29.0 * sign_y), (25.0 * sign_x, 30.0 * sign_y))
    horizontal = ((-19.0 * sign_x, 30.0 * sign_y), (25.0 * sign_x, 30.0 * sign_y))
    candidate = corner_perimeter_candidate(
        30.0, 50.0, 4.0, vertical, horizontal, corner_style=style, quad_segs=1,
    )
    perimeter = candidate.perimeter

    points = [(25.0, -29.0), (-19.0, -29.0), (-19.0, 30.0)]
    if style == "round":
        points = [(25.0, -29.0), (-15.0, -29.0), (-19.0, -25.0), (-19.0, 30.0)]
    expected = LineString([(sign_x * x, sign_y * y) for x, y in points])

    assert isinstance(perimeter, LineString)
    assert perimeter.is_valid and perimeter.is_simple
    assert not perimeter.is_closed
    assert perimeter.equals(expected)
    assert candidate.reached_edges == (vertical, horizontal)
    assert len(segments_of(perimeter)) == (2 if style == "sharp" else 3)
    corner_length = 8.0 if style == "sharp" else 4.0 * 2.0**0.5
    assert perimeter.length == pytest.approx(95.0 + corner_length)
    for edge in (vertical, horizontal):
        assert perimeter.intersection(LineString(edge)).length == 0.0

    reversed_vertical = (vertical[1], vertical[0])
    reversed_horizontal = (horizontal[1], horizontal[0])
    reversed_candidate = corner_perimeter_candidate(
        30.0, 50.0, 4.0, reversed_vertical, reversed_horizontal,
        corner_style=style, quad_segs=1,
    )
    assert reversed_candidate.perimeter.equals(perimeter)
    assert reversed_candidate.reached_edges == (reversed_vertical, reversed_horizontal)


@pytest.mark.parametrize("offset", [0.0, 4.0])
@pytest.mark.parametrize("gap_x,gap_y", [(0.0, 0.0), (2.0, 3.0), (4.0, 4.0), (10.0, 5.0)])
def test_corner_candidate_reaches_near_and_far_edges(offset: float, gap_x: float, gap_y: float) -> None:
    x, y = 15.0 + gap_x, 25.0 + gap_y
    left, bottom = -15.0 - offset, -25.0 - offset
    vertical = ((x, bottom), (x, y))
    horizontal = ((left, y), (x, y))

    candidate = corner_perimeter_candidate(30.0, 50.0, offset, vertical, horizontal)

    assert candidate.perimeter.equals(LineString([(x, bottom), (left, bottom), (left, y)]))
    assert len(segments_of(candidate.perimeter)) == 2


def test_corner_candidate_round_with_fine_discretization() -> None:
    candidate = corner_perimeter_candidate(
        30.0, 50.0, 4.0, ((25.0, -29.0), (25.0, 30.0)), ((-19.0, 30.0), (25.0, 30.0)),
        corner_style="round", quad_segs=32,
    )
    # Two straight sections and one quarter circle of radius 4.
    exact_arc_length = 95.0 + 2.0 * pi
    assert candidate.perimeter.length < exact_arc_length
    assert candidate.perimeter.length == pytest.approx(exact_arc_length, rel=1e-5)
    assert len(segments_of(candidate.perimeter)) == 34
    assert candidate.perimeter.bounds == pytest.approx((-19.0, -29.0, 25.0, 30.0))

@pytest.mark.parametrize("bad_edge,message", [
    ((), "two 2D endpoints"),
    (((0.0, 0.0),), "two 2D endpoints"),
    (((0.0, 0.0), (1.0, 1.0), (2.0, 2.0)), "two 2D endpoints"),
    (((0.0, 0.0, 0.0), (1.0, 1.0, 0.0)), "two 2D endpoints"),
    (((float("nan"), 0.0), (1.0, 1.0)), "coordinates must be finite"),
    (((0.0, float("inf")), (1.0, 1.0)), "coordinates must be finite"),
    (((0.0, 0.0), (float("-inf"), 1.0)), "coordinates must be finite"),
])
def test_corner_candidate_rejects_bad_edge_data(bad_edge: tuple[tuple[float, ...], ...], message: str) -> None:
    vertical = ((25.0, -29.0), (25.0, 30.0))
    horizontal = ((-19.0, 30.0), (25.0, 30.0))
    with pytest.raises(ValueError, match=message):
        corner_perimeter_candidate(30.0, 50.0, 4.0, bad_edge, horizontal)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match=message):
        corner_perimeter_candidate(30.0, 50.0, 4.0, vertical, bad_edge)  # type: ignore[arg-type]


@pytest.mark.parametrize("vertical,horizontal,message", [
    (((25.0, 30.0), (25.0, 30.0)), ((-19.0, 30.0), (25.0, 30.0)), "nonzero vertical"),
    (((25.0, -29.0), (25.0, 30.0)), ((25.0, 30.0), (25.0, 30.0)), "nonzero vertical"),
    (((0.0, 30.0), (25.0, 30.0)), ((-19.0, 30.0), (25.0, 30.0)), "nonzero vertical"),
    (((25.0, -29.0), (25.0, 30.0)), ((25.0, 0.0), (25.0, 30.0)), "nonzero vertical"),
    (((25.0, -29.0), (26.0, 30.0)), ((-19.0, 30.0), (25.0, 30.0)), "nonzero vertical"),
    (((25.0, -29.0), (25.0, 30.0)), ((-19.0, 31.0), (25.0, 30.0)), "nonzero vertical"),
    (((10.0, -29.0), (10.0, 30.0)), ((-19.0, 30.0), (25.0, 30.0)), "column interior"),
    (((25.0, -29.0), (25.0, 30.0)), ((-19.0, 20.0), (25.0, 20.0)), "column interior"),
    (((25.0, -29.0), (25.0, 31.0)), ((-19.0, 30.0), (25.0, 30.0)), "share their corner endpoint"),
    (((25.0, -29.0), (25.0, 30.0)), ((-19.0, 30.0), (26.0, 30.0)), "share their corner endpoint"),
    (((25.0, -28.0), (25.0, 30.0)), ((-19.0, 30.0), (25.0, 30.0)), "contain its contour endpoint"),
    (((25.0, -29.0), (25.0, 30.0)), ((-18.0, 30.0), (25.0, 30.0)), "contain its contour endpoint"),
])
def test_corner_candidate_rejects_incompatible_edges(vertical: Segment, horizontal: Segment, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        corner_perimeter_candidate(30.0, 50.0, 4.0, vertical, horizontal)


@pytest.mark.filterwarnings("ignore:overflow encountered in buffer:RuntimeWarning")
def test_corner_candidate_rejects_numeric_overflow() -> None:
    vertical = ((6.5e307, -0.5), (6.5e307, 1.0))
    horizontal = ((-6.5e307, 1.0), (6.5e307, 1.0))
    with pytest.raises(ValueError, match="Geometry exceeds the finite numeric range"):
        corner_perimeter_candidate(1.3e308, 1.0, 0.0, vertical, horizontal)

@pytest.mark.parametrize("opened", [
    LineString(),
    LineString([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]),
    MultiLineString([[(0.0, 0.0), (1.0, 0.0)], [(2.0, 0.0), (3.0, 0.0)]]),
])
def test_corner_candidate_rejects_unexpected_clip_result(
    monkeypatch: pytest.MonkeyPatch,
    opened: LineString | MultiLineString,
) -> None:
    monkeypatch.setattr("mento.punching_geometry.clip_perimeter_to_slab", lambda *_args: opened)
    with pytest.raises(ValueError, match="single open contour"):
        corner_perimeter_candidate(
            30.0, 50.0, 4.0, ((25.0, -29.0), (25.0, 30.0)), ((-19.0, 30.0), (25.0, 30.0)),
        )


@pytest.mark.parametrize("opened", [
    LineString([(0.0, 0.0), (1.0, 1.0)]),  # No end on x=15.
    LineString([(15.0, 0.0), (15.0, 1.0)]),  # Both ends on x=15.
    LineString([(15.0, -29.0), (-19.0, 24.0)]),  # No end on y=25.
    LineString([(15.0, 25.0), (-19.0, 25.0)]),  # Both ends on y=25.
])
def test_corner_candidate_rejects_ambiguous_ends(monkeypatch: pytest.MonkeyPatch, opened: LineString) -> None:
    monkeypatch.setattr("mento.punching_geometry.clip_perimeter_to_slab", lambda *_args: opened)
    with pytest.raises(ValueError, match="end facing each edge"):
        corner_perimeter_candidate(
            30.0, 50.0, 4.0, ((25.0, -29.0), (25.0, 30.0)), ((-19.0, 30.0), (25.0, 30.0)),
        )

@pytest.mark.parametrize("style", ["sharp", "round"])
@pytest.mark.parametrize("bounds,expected_keys", [
    ((-100.0, -100.0, 100.0, 100.0), {
        (0, 0), (-1, 0), (1, 0), (0, -1), (0, 1),
        (-1, -1), (-1, 1), (1, -1), (1, 1),
    }),
    ((-100.0, -100.0, 15.0, 100.0), {(1, 0), (1, -1), (1, 1)}),
    ((-15.0, -100.0, 100.0, 100.0), {(-1, 0), (-1, -1), (-1, 1)}),
    ((-100.0, -100.0, 100.0, 25.0), {(0, 1), (-1, 1), (1, 1)}),
    ((-100.0, -25.0, 100.0, 100.0), {(0, -1), (-1, -1), (1, -1)}),
    ((-100.0, -100.0, 15.0, 25.0), {(1, 1)}),
    ((-15.0, -100.0, 100.0, 25.0), {(-1, 1)}),
    ((-15.0, -25.0, 100.0, 100.0), {(-1, -1)}),
    ((-100.0, -25.0, 15.0, 100.0), {(1, -1)}),
    ((-15.0, -100.0, 15.0, 100.0), set()),
    ((-100.0, -25.0, 100.0, 25.0), set()),
    ((-19.0, -29.0, 19.0, 29.0), set()),
])
def test_generator_selects_geometrically_compatible_candidates(
    style: Literal["sharp", "round"],
    bounds: tuple[float, float, float, float],
    expected_keys: set[tuple[int, int]],
) -> None:
    slab_outline = box(*bounds)
    candidates = rectangular_slab_perimeter_candidates(
        30.0, 50.0, 4.0, slab_outline, corner_style=style, quad_segs=8,
    )
    found: set[tuple[int, int]] = set()
    for candidate in candidates:
        # -1: left/bottom; +1: right/top; 0: no edge on that axis.
        sign_x = sign_y = 0
        for edge in candidate.reached_edges:
            assert slab_outline.boundary.covers(LineString(edge))
            if edge[0][0] == edge[1][0]:
                sign_x = 1 if edge[0][0] > 0 else -1
            else:
                sign_y = 1 if edge[0][1] > 0 else -1
        key = (sign_x, sign_y)
        assert key not in found
        found.add(key)

        assert not candidate.perimeter.is_empty
        assert slab_outline.covers(candidate.perimeter)
        assert candidate.perimeter.intersection(slab_outline.boundary).length == 0.0
        assert candidate.perimeter.is_closed == (key == (0, 0))

    assert found == expected_keys

@pytest.mark.parametrize("limit,expected_lengths", [
    (0.0, []),
    (152.0, []),
    (153.0, [153.0]),
    (156.0, [153.0, 156.0]),
    (192.0, [153.0, 156.0, 178.0, 192.0]),
])
def test_generator_applies_explicit_length_limit(limit: float, expected_lengths: list[float]) -> None:
    candidates = rectangular_slab_perimeter_candidates(
        30.0, 50.0, 4.0, box(-100.0, -100.0, 75.0, 30.0), max_length=limit,
    )
    assert sorted(c.perimeter.length for c in candidates) == pytest.approx(expected_lengths)


def test_generator_keeps_corner_when_single_edge_is_too_long() -> None:
    candidates = rectangular_slab_perimeter_candidates(
        30.0, 50.0, 4.0, box(-100.0, -100.0, 75.0, 30.0), max_length=192.0,
    )
    right = ((75.0, -100.0), (75.0, 30.0))
    top = ((-100.0, 30.0), (75.0, 30.0))
    assert not any(c.reached_edges == (right,) for c in candidates)
    assert any(c.reached_edges == (right, top) for c in candidates)


@pytest.mark.parametrize("quad_segs", [1, 8, 32])
def test_generator_passes_rounding_resolution(quad_segs: int) -> None:
    candidates = rectangular_slab_perimeter_candidates(
        30.0, 50.0, 4.0, box(-100.0, -100.0, 100.0, 100.0),
        corner_style="round", quad_segs=quad_segs,
    )
    assert len(candidates) == 9
    expected_segments = {0: 4 + 4 * quad_segs, 1: 3 + 2 * quad_segs, 2: 2 + quad_segs}
    for candidate in candidates:
        assert len(segments_of(candidate.perimeter)) == expected_segments[len(candidate.reached_edges)]

@pytest.mark.parametrize("slab_outline,message", [
    (Polygon(), "nonempty valid polygon"),
    (Polygon([(-100, -100), (100, 100), (-100, 100), (100, -100)]), "nonempty valid polygon"),
    (Polygon([(-100, -100), (100, -100), (0, 100)]), "axis-aligned rectangular slab"),
    (Polygon([(-100, 0), (0, -100), (100, 0), (0, 100)]), "axis-aligned rectangular slab"),
    (
        Polygon(box(-100, -100, 100, 100).exterior.coords, holes=[box(60, 60, 80, 80).exterior.coords]),
        "Define openings separately",
    ),
    (
        Polygon([(-100, -100, 0), (100, -100, 0), (100, 100, 0), (-100, 100, 0)]),
        "two-dimensional",
    ),
    (box(-10, -100, 100, 100), "complete column footprint"),
    (box(-100, -100, 10, 100), "complete column footprint"),
    (box(-100, -20, 100, 100), "complete column footprint"),
    (box(-100, -100, 100, 20), "complete column footprint"),
])
def test_generator_rejects_invalid_or_unsupported_slab(slab_outline: Polygon, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        rectangular_slab_perimeter_candidates(30.0, 50.0, 4.0, slab_outline)


def test_generator_rejects_line_as_slab() -> None:
    with pytest.raises(TypeError, match="Slab outline must be a Polygon"):
        rectangular_slab_perimeter_candidates(
            30.0, 50.0, 4.0, LineString([(0.0, 0.0), (100.0, 100.0)]),  # type: ignore[arg-type]
        )


@pytest.mark.parametrize("limit", [-1.0, float("nan"), float("inf"), float("-inf")])
def test_generator_rejects_invalid_length_limit(limit: float) -> None:
    with pytest.raises(ValueError, match="Maximum length"):
        rectangular_slab_perimeter_candidates(
            30.0, 50.0, 4.0, box(-100.0, -100.0, 100.0, 100.0), max_length=limit,
        )


def test_generator_accepts_extra_collinear_slab_vertices() -> None:
    slab_outline = Polygon([
        (-100.0, -100.0), (0.0, -100.0), (75.0, -100.0),
        (75.0, 0.0), (75.0, 30.0), (-100.0, 30.0),
    ])
    actual = rectangular_slab_perimeter_candidates(30.0, 50.0, 4.0, slab_outline)
    expected = rectangular_slab_perimeter_candidates(30.0, 50.0, 4.0, box(-100.0, -100.0, 75.0, 30.0))

    assert len(actual) == len(expected) == 9
    for candidate, reference in zip(actual, expected):
        assert candidate.perimeter.equals(reference.perimeter)
        assert candidate.reached_edges == reference.reached_edges

@pytest.mark.parametrize("scale", [1e-6, 1.0, 1e6])
@pytest.mark.parametrize("x,y,direction", [
    (40.0, 0.0, 0.0),
    (0.0, 40.0, pi / 2),
    (-40.0, 0.0, pi),
    (0.0, -40.0, -pi / 2),
])
def test_rectangular_shadow_directions_and_scale(
    x: float, y: float, direction: float, scale: float,
) -> None:
    angles = rectangular_opening_shadow_angles(
        x * scale, y * scale, 20.0 * scale, 20.0 * scale,
    )
    half_angle = atan2(10.0, 30.0)

    assert angles == pytest.approx(
        (direction - half_angle, direction + half_angle)
    )
    assert 0 < angles[1] - angles[0] < pi


def test_rectangular_shadow_rotated_45_degrees() -> None:
    # Rotate both the center (40, 0) and a 40 x 20 opening by 45 degrees.
    coordinate = 40.0 / (2.0 ** 0.5)
    angles = rectangular_opening_shadow_angles(
        coordinate, coordinate, 40.0, 20.0, rotation=pi / 4,
    )
    half_angle = atan2(10.0, 20.0)

    assert angles == pytest.approx(
        (pi / 4 - half_angle, pi / 4 + half_angle)
    )


@pytest.mark.parametrize("name", ["x", "y", "b", "h", "rotation"])
@pytest.mark.parametrize("value", [
    float("nan"), float("inf"), float("-inf"),
])
def test_rectangular_shadow_rejects_nonfinite_values(
    name: str, value: float,
) -> None:
    arguments = {
        "x": 40.0, "y": 0.0, "b": 20.0, "h": 20.0, "rotation": 0.0,
    }
    arguments[name] = value

    with pytest.raises(ValueError, match="must be finite"):
        rectangular_opening_shadow_angles(**arguments)


@pytest.mark.parametrize("b,h", [
    (0.0, 20.0),
    (-1.0, 20.0),
    (20.0, 0.0),
    (20.0, -1.0),
])
def test_rectangular_shadow_rejects_nonpositive_dimensions(
    b: float, h: float,
) -> None:
    with pytest.raises(ValueError, match="must be positive"):
        rectangular_opening_shadow_angles(40.0, 0.0, b, h)


@pytest.mark.parametrize("x,y,rotation", [
    (0.0, 0.0, 0.0),           # Origin at the opening center.
    (1.0, 0.0, 0.0),           # Origin inside the opening.
    (2.0, 0.0, 0.0),           # Origin on a side.
    (2.0, 1.0, 0.0),           # Origin on a corner.
    (0.0, 2.0, pi / 2),        # Contact after rotation.
    (2.0 + 1e-15, 0.0, 0.0),   # Numerically indistinguishable contact.
])
def test_rectangular_shadow_rejects_origin_contact(
    x: float, y: float, rotation: float,
) -> None:
    with pytest.raises(
        ValueError, match="contains or touches the column center",
    ):
        rectangular_opening_shadow_angles(
            x, y, 4.0, 2.0, rotation=rotation,
        )


@pytest.mark.parametrize("x,y,b,h", [
    (1e20, 0.0, 1.0, 20.0),       # Both x bounds round to the same float.
    (0.0, 1e20, 20.0, 1.0),       # Both y bounds round to the same float.
    (1e308, 0.0, 1.7e308, 20.0),  # A computed bound overflows.
])
def test_rectangular_shadow_rejects_unrepresentable_bounds(
    x: float, y: float, b: float, h: float,
) -> None:
    with pytest.raises(ValueError, match="bounds cannot be represented"):
        rectangular_opening_shadow_angles(x, y, b, h)


def test_rectangular_shadow_rejects_unresolvable_angles() -> None:
    # Bounds remain distinct, but both tangent angles round to the same value.
    with pytest.raises(ValueError, match="tangents cannot be resolved"):
        rectangular_opening_shadow_angles(0.0, 1e20, 20.0, 1e5)

@pytest.mark.parametrize("scale", [1e-6, 1.0, 1e6])
@pytest.mark.parametrize("x,y", [
    (40.0, 0.0),
    (0.0, 40.0),
    (-40.0, 0.0),
    (0.0, -40.0),
    (24.0, 32.0),
])
def test_circular_shadow_directions_and_scale(
    x: float, y: float, scale: float,
) -> None:
    # Every center is 40 units from the origin; the radius is 20.
    angles = circular_opening_shadow_angles(
        x * scale, y * scale, 20.0 * scale,
    )
    direction = atan2(y, x)
    half_angle = pi / 6  # sin(30 degrees) = 20 / 40.

    assert angles == pytest.approx(
        (direction - half_angle, direction + half_angle)
    )
    assert 0 < angles[1] - angles[0] < pi


@pytest.mark.parametrize("name", ["x", "y", "radius"])
@pytest.mark.parametrize("value", [
    float("nan"), float("inf"), float("-inf"),
])
def test_circular_shadow_rejects_nonfinite_values(
    name: str, value: float,
) -> None:
    arguments = {"x": 40.0, "y": 0.0, "radius": 20.0}
    arguments[name] = value

    with pytest.raises(ValueError, match="must be finite"):
        circular_opening_shadow_angles(**arguments)


@pytest.mark.parametrize("radius", [0.0, -1.0])
def test_circular_shadow_rejects_nonpositive_radius(radius: float) -> None:
    with pytest.raises(ValueError, match="radius must be positive"):
        circular_opening_shadow_angles(40.0, 0.0, radius)


@pytest.mark.parametrize("x,y", [
    (0.0, 0.0),           # Origin at the opening center.
    (10.0, 0.0),          # Origin inside the opening.
    (20.0, 0.0),          # Origin on the circumference.
    (12.0, 16.0),         # Contact away from the coordinate axes.
    (20.0 + 1e-14, 0.0),  # Numerically indistinguishable contact.
])
def test_circular_shadow_rejects_origin_contact(x: float, y: float) -> None:
    with pytest.raises(
        ValueError, match="contains or touches the column center",
    ):
        circular_opening_shadow_angles(x, y, 20.0)


def test_circular_shadow_accepts_resolvable_gap_from_origin() -> None:
    start, end = circular_opening_shadow_angles(
        20.0 + 1e-8, 0.0, 20.0,
    )

    assert -pi / 2 < start < 0.0 < end < pi / 2
    assert start == pytest.approx(-end)
    assert 0 < end - start < pi


def test_circular_shadow_rejects_distance_overflow() -> None:
    # The inputs are finite, but their combined distance exceeds float range.
    with pytest.raises(ValueError, match="distance cannot be represented"):
        circular_opening_shadow_angles(1.3e308, 1.3e308, 20.0)


@pytest.mark.parametrize("x,y,radius", [
    (0.0, 1e20, 1.0),      # Both tangent angles round to pi/2.
    (1e308, 0.0, 1e-100),  # The radius/distance ratio underflows to zero.
])
def test_circular_shadow_rejects_unresolvable_angles(
    x: float, y: float, radius: float,
) -> None:
    with pytest.raises(ValueError, match="tangents cannot be resolved"):
        circular_opening_shadow_angles(x, y, radius)

@pytest.mark.parametrize("scale", [1e-6, 1.0, 1e6])
def test_shadow_cut_splits_line_without_bridging_gap(scale: float) -> None:
    original = LineString([
        (20 * scale, -30 * scale),
        (20 * scale, 30 * scale),
    ])
    before = original.wkb
    reduced = subtract_opening_shadow(original, -pi / 4, pi / 4)

    assert isinstance(reduced, MultiLineString)
    parts = sorted(reduced.geoms, key=lambda part: part.bounds[1])

    assert len(parts) == 2
    assert parts[0].bounds == pytest.approx(
        tuple(v * scale for v in (20, -30, 20, -20))
    )
    assert parts[1].bounds == pytest.approx(
        tuple(v * scale for v in (20, 20, 20, 30))
    )
    assert reduced.length == pytest.approx(20 * scale)
    assert all(not part.is_closed for part in parts)
    assert original.wkb == before


def test_shadow_cut_preserves_disconnected_input() -> None:
    original = MultiLineString([
        [(20, -30), (20, 30)],
        [(-20, -30), (-20, 30)],
    ])
    before = original.wkb
    reduced = subtract_opening_shadow(original, -pi / 4, pi / 4)

    assert isinstance(reduced, MultiLineString)
    assert len(reduced.geoms) == 3
    assert reduced.length == pytest.approx(80.0)

    left = LineString([(-20, -30), (-20, 30)])
    assert any(part.equals(left) for part in reduced.geoms)
    assert original.wkb == before


def test_shadow_cut_leaves_geometry_outside_shadow_unchanged() -> None:
    original = LineString([(-20, -30), (-20, 30)])
    reduced = subtract_opening_shadow(original, -pi / 4, pi / 4)

    assert reduced.equals(original)


@pytest.mark.parametrize("width", [pi / 6, 2 * pi / 3, pi - 1e-8])
def test_shadow_cut_removes_entire_line_even_for_wide_shadow(
    width: float,
) -> None:
    original = LineString([(100, -1), (100, 1)])
    reduced = subtract_opening_shadow(original, -width / 2, width / 2)

    assert isinstance(reduced, LineString)
    assert reduced.is_empty
    assert original.length == pytest.approx(2.0)


@pytest.mark.parametrize("original", [LineString(), MultiLineString()])
def test_shadow_cut_accepts_empty_input(
    original: LineString | MultiLineString,
) -> None:
    reduced = subtract_opening_shadow(original, -pi / 4, pi / 4)

    assert isinstance(reduced, LineString)
    assert reduced.is_empty


def test_shadow_cut_is_idempotent() -> None:
    original = LineString([(20, -30), (20, 30)])
    once = subtract_opening_shadow(original, -pi / 4, pi / 4)
    twice = subtract_opening_shadow(once, -pi / 4, pi / 4)

    assert twice.length == pytest.approx(once.length)
    assert twice.hausdorff_distance(once) < 1e-10


@pytest.mark.parametrize("turns", [-2, 1, 3])
def test_shadow_cut_accepts_equivalent_angles(turns: int) -> None:
    original = LineString([(20, -30), (20, 30)])
    reference = subtract_opening_shadow(original, -pi / 4, pi / 4)
    shift = turns * 2 * pi
    reduced = subtract_opening_shadow(
        original, -pi / 4 + shift, pi / 4 + shift,
    )

    assert reduced.length == pytest.approx(reference.length)
    assert reduced.hausdorff_distance(reference) < 1e-10


@pytest.mark.parametrize("original", [None, box(-1, -1, 1, 1)])
def test_shadow_cut_rejects_non_linear_input(original: object) -> None:
    with pytest.raises(
        TypeError, match="Perimeter must be a LineString or MultiLineString",
    ):
        subtract_opening_shadow(
            original, -pi / 4, pi / 4,  # type: ignore[arg-type]
        )


@pytest.mark.parametrize("original,message", [
    (LineString([(1, 1, 0), (2, 2, 0)]), "two-dimensional"),
    (LineString([(1, 1), (1, 1)]), "geometry must be valid"),
])
def test_shadow_cut_rejects_invalid_geometry(
    original: LineString, message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        subtract_opening_shadow(original, -pi / 4, pi / 4)


@pytest.mark.parametrize("index", [0, 1])
@pytest.mark.parametrize("value", [
    float("nan"), float("inf"), float("-inf"),
])
def test_shadow_cut_rejects_nonfinite_angles(
    index: int, value: float,
) -> None:
    angles = [-pi / 4, pi / 4]
    angles[index] = value

    with pytest.raises(ValueError, match="angles must be finite"):
        subtract_opening_shadow(
            LineString([(20, -30), (20, 30)]), *angles,
        )


@pytest.mark.parametrize("start,end", [
    (0.0, 0.0),
    (1.0, 0.0),
    (0.0, pi),
    (0.0, 2 * pi),
    (-1e308, 1e308),
])
def test_shadow_cut_rejects_invalid_angular_width(
    start: float, end: float,
) -> None:
    with pytest.raises(
        ValueError, match="angular width between zero and pi",
    ):
        subtract_opening_shadow(
            LineString([(20, -30), (20, 30)]), start, end,
        )


def test_shadow_cut_rejects_extent_overflow() -> None:
    original = LineString([(1e308, 0.0), (1.1e308, 0.0)])

    with pytest.raises(ValueError, match="extent cannot be represented"):
        subtract_opening_shadow(original, -pi / 4, pi / 4)


def test_shadow_cut_rejects_numerically_degenerate_polygon() -> None:
    # Tiny coordinates and angular width make every shadow vertex lie on y=0.
    original = LineString([(1e-160, 0.0), (2e-160, 0.0)])

    with pytest.raises(
        ValueError, match="shadow polygon could not be constructed",
    ):
        subtract_opening_shadow(original, 0.0, 1e-200)