"""Geometric properties of punching shear critical sections.

All lengths must use the same unit. Inputs and outputs are plain floats.
The caller defines the perimeter and supplies the effective depth.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import asin, atan2, cos, hypot, isclose, isfinite, pi, sin, ulp
from typing import List, Literal, Sequence, Tuple
from shapely import line_merge
from shapely.geometry import LineString, MultiLineString, Polygon, box

XY = Tuple[float, float]
Segment = Tuple[XY, XY]

@dataclass(frozen=True)
class PunchingPerimeterCandidate:
    """A candidate contour and its reached free slab edges.

    perimeter is the contour before reductions due to openings.
    reached_edges stores the slab edge segments where the contour ends,
    using the same coordinates and length unit as the perimeter.
    """

    perimeter: LineString | MultiLineString
    reached_edges: Tuple[Segment, ...] = ()

@dataclass(frozen=True)
class PunchingSectionProperties:
    """Properties of a punching shear critical section.

    J_x and J_y include in-plane and through-depth contributions.
    J_xy has no through-depth contribution.
    """

    b_0: float
    A_c: float
    x_g: float
    y_g: float
    J_x_plan: float
    J_y_plan: float
    dJ_x: float
    dJ_y: float
    J_x: float
    J_y: float
    J_xy: float
    vertices: Tuple[XY, ...]
    segments: Tuple[Segment, ...]
    parts: int

    def arms(self, point: XY) -> XY:
        """Signed distances from the section centroid."""
        return (point[0] - self.x_g, point[1] - self.y_g)

    @property
    def extents(self) -> XY:
        """Overall widths along x and y."""
        xs = [p[0] for p in self.vertices]
        ys = [p[1] for p in self.vertices]
        return (max(xs) - min(xs), max(ys) - min(ys))


def section_properties(segments: Sequence[Segment], d: float, parts: int = 1) -> PunchingSectionProperties:
    """Calculate properties by integrating along straight segments.

    Lengths use the input unit, A_c its square, and J its fourth power.
    Through-depth terms retain the projection rule of the local development:
    dJ_x = (d**3 / 12) * sum(abs(delta_y)), and likewise for dJ_y.
    """
    segments = tuple(segments)
    if not segments:
        raise ValueError("A critical section needs at least one segment.")
    if not isfinite(d) or d <= 0:
        raise ValueError("Effective depth must be finite and greater than zero.")

    # 1. Perimeter length and centroid.
    length = sx = sy = 0.0
    for (x1, y1), (x2, y2) in segments:
        if not all(isfinite(value) for value in (x1, y1, x2, y2)):
            raise ValueError("Segment coordinates must be finite.")
        ell = hypot(x2 - x1, y2 - y1)
        length += ell
        sx += ell * (x1 + x2) / 2
        sy += ell * (y1 + y2) / 2

    if length <= 0:
        raise ValueError("The critical perimeter must have positive length.")

    x_g, y_g = sx / length, sy / length

    # 2. In-plane integrals about the section centroid.
    j_x = j_y = j_xy = proj_y = proj_x = 0.0
    for (x1, y1), (x2, y2) in segments:
        ell = hypot(x2 - x1, y2 - y1)
        u1, v1 = x1 - x_g, y1 - y_g
        u2, v2 = x2 - x_g, y2 - y_g

        j_x += ell / 3 * (v1 * v1 + v1 * v2 + v2 * v2)
        j_y += ell / 3 * (u1 * u1 + u1 * u2 + u2 * u2)
        j_xy += ell / 6 * (2 * u1 * v1 + u1 * v2 + u2 * v1 + 2 * u2 * v2)
        proj_y += abs(y2 - y1)
        proj_x += abs(x2 - x1)

    j_x, j_y, j_xy = d * j_x, d * j_y, d * j_xy

    # 3. Through-depth contributions, kept separate for review.
    dj_x = d**3 / 12 * proj_y
    dj_y = d**3 / 12 * proj_x

    # 4. Unique vertices for subsequent stress evaluations.
    vertices: List[XY] = []
    for a, b in segments:
        for point in (a, b):
            rounded = (round(point[0], 9), round(point[1], 9))
            if rounded not in vertices:
                vertices.append(rounded)

    return PunchingSectionProperties(
        b_0=length,
        A_c=length * d,
        x_g=x_g,
        y_g=y_g,
        J_x_plan=j_x,
        J_y_plan=j_y,
        dJ_x=dj_x,
        dJ_y=dj_y,
        J_x=j_x + dj_x,
        J_y=j_y + dj_y,
        J_xy=j_xy,
        vertices=tuple(vertices),
        segments=segments,
        parts=parts,
    )

def column_offset_perimeter(
    c_x: float,
    c_y: float,
    offset: float,
    *,
    corner_style: Literal["sharp", "round"] = "sharp",
    quad_segs: int = 32,
) -> LineString:
    """Offset a rectangular column centered at the origin.

    All lengths use the same unit. The result is a closed contour,
    without cuts for slab edges or openings.
    Rounded corners use quad_segs straight segments per quarter circle.
    """
    if not all(isfinite(value) and value > 0 for value in (c_x, c_y)):
        raise ValueError("Column dimensions must be finite and positive.")

    if not isfinite(offset) or offset < 0:
        raise ValueError("Offset must be finite and nonnegative.")

    if corner_style not in ("sharp", "round"):
        raise ValueError("Corner style must be 'sharp' or 'round'.")

    if isinstance(quad_segs, bool) or not isinstance(quad_segs, int) or quad_segs < 1:
        raise ValueError("quad_segs must be a positive integer.")

    column = box(-c_x / 2, -c_y / 2, c_x / 2, c_y / 2)

    expanded = column.buffer(
        offset,
        join_style="mitre" if corner_style == "sharp" else "round",
        quad_segs=quad_segs,
    )

    return LineString(expanded.exterior.coords)

def segments_of(perimeter: LineString | MultiLineString) -> List[Segment]:
    """Extract 2D segments without closing or joining separate lines."""
    lines: Sequence[LineString]

    if isinstance(perimeter, LineString):
        lines = [perimeter]
    elif isinstance(perimeter, MultiLineString):
        lines = list(perimeter.geoms)
    else:
        raise TypeError("Expected a LineString or MultiLineString.")

    segments: List[Segment] = []

    for line in lines:
        coords = list(line.coords)

        if any(len(point) != 2 for point in coords):
            raise ValueError("Perimeter coordinates must be two-dimensional.")

        for p1, p2 in zip(coords, coords[1:]):
            a = (float(p1[0]), float(p1[1]))
            b = (float(p2[0]), float(p2[1]))

            if a != b:
                segments.append((a, b))

    return segments

def clip_perimeter_to_slab(
    perimeter: LineString | MultiLineString,
    slab_outline: Polygon,
) -> LineString | MultiLineString:
    """Keep perimeter lines inside a slab with free exterior edges.

    Both geometries must be 2D and use the same origin and length unit.
    The slab polygon must have no holes; openings are handled separately.
    This only clips: it does not extend lines or select critical candidates.
    """
    if not isinstance(perimeter, (LineString, MultiLineString)):
        raise TypeError("Perimeter must be a LineString or MultiLineString.")

    if not isinstance(slab_outline, Polygon):
        raise TypeError("Slab outline must be a Polygon.")

    if slab_outline.is_empty or not slab_outline.is_valid:
        raise ValueError("Slab outline must be a nonempty valid polygon.")

    if len(slab_outline.interiors) > 0:
        raise ValueError("Define openings separately from the slab outline.")

    if any(g.has_z or g.has_m for g in (perimeter, slab_outline)):
        raise ValueError("Both geometries must be two-dimensional.")

    if not perimeter.is_valid:
        raise ValueError("Perimeter geometry must be valid.")

    clipped = perimeter.intersection(slab_outline)
    clipped = clipped.difference(slab_outline.boundary)
    merged = line_merge(clipped)

    return LineString() if merged.is_empty else merged

def candidate_fits_length_limit(
    c_x: float,
    c_y: float,
    offset: float,
    *,
    max_length: float | None = None,
    face_gap_x: float | None = None,
    face_gap_y: float | None = None,
    corner_style: Literal["sharp", "round"] = "sharp",
    quad_segs: int = 32,
) -> bool:
    """Check a rectangular candidate's length before building its contour.

    Selected edges are parallel to the column faces, at most one per axis.
    Gaps are measured from column faces; None means no edge on that axis
    is selected for this candidate. Round corners use polygonal arcs.
    No length limit is applied when max_length is None. This does not
    check whether a candidate fits the slab or is admissible under a code.
    """
    if not all(isfinite(value) and value > 0 for value in (c_x, c_y)):
        raise ValueError("Column dimensions must be finite and positive.")

    if not isfinite(offset) or offset < 0:
        raise ValueError("Offset must be finite and nonnegative.")

    for gap in (face_gap_x, face_gap_y):
        if gap is not None and (not isfinite(gap) or gap < 0):
            raise ValueError("Face gaps must be finite and nonnegative.")

    if max_length is not None and (not isfinite(max_length) or max_length < 0):
        raise ValueError("Maximum length must be finite and nonnegative.")

    if corner_style not in ("sharp", "round"):
        raise ValueError("Corner style must be 'sharp' or 'round'.")

    if isinstance(quad_segs, bool) or not isinstance(quad_segs, int) or quad_segs < 1:
        raise ValueError("quad_segs must be a positive integer.")

    corner_length = 2 * offset
    if corner_style == "round":
        corner_length = offset * (2 * quad_segs * sin(pi / (4 * quad_segs)))

    if face_gap_x is not None and face_gap_y is not None:
        length = c_x + c_y + face_gap_x + face_gap_y + corner_length
    elif face_gap_x is not None:
        length = 2 * c_x + c_y + 2 * face_gap_x + 2 * corner_length
    elif face_gap_y is not None:
        length = c_x + 2 * c_y + 2 * face_gap_y + 2 * corner_length
    else:
        length = 2 * (c_x + c_y) + 4 * corner_length

    if not isfinite(length):
        raise ValueError("Candidate length exceeds the finite numeric range.")

    return (
        max_length is None
        or length <= max_length
        or isclose(length, max_length, rel_tol=1e-12, abs_tol=0.0)
    )

def edge_perimeter_candidate(
    c_x: float,
    c_y: float,
    offset: float,
    edge: Segment,
    *,
    corner_style: Literal["sharp", "round"] = "sharp",
    quad_segs: int = 32,
) -> PunchingPerimeterCandidate:
    """Build an open contour reaching one free slab edge.

    The rectangular column is centered at the origin. The edge must be
    horizontal or vertical, outside or on a column face, and long enough
    to contain both contour endpoints. All coordinates use the same unit.
    Other slab edges, openings and candidate selection are handled later.
    """
    closed = column_offset_perimeter(
        c_x, c_y, offset, corner_style=corner_style, quad_segs=quad_segs
    )

    if len(edge) != 2 or any(len(point) != 2 for point in edge):
        raise ValueError("An edge needs two 2D endpoints.")
    if not all(isfinite(value) for point in edge for value in point):
        raise ValueError("Edge coordinates must be finite.")

    (x1, y1), (x2, y2) = edge
    if x1 == x2 and y1 != y2:
        axis = 0
    elif y1 == y2 and x1 != x2:
        axis = 1
    else:
        raise ValueError("The edge must be nonzero and axis-aligned.")

    edge_coordinate = edge[0][axis]
    half_size = (c_x, c_y)[axis] / 2
    if abs(edge_coordinate) < half_size:
        raise ValueError("The edge cannot cross the column interior.")

    # Open the closed offset at the column face, not at the actual edge.
    sign = 1 if edge_coordinate > 0 else -1
    xmin, ymin, xmax, ymax = closed.bounds
    margin = max(c_x, c_y, offset)
    limits = [xmin - margin, ymin - margin, xmax + margin, ymax + margin]
    limits[axis + 2 if sign > 0 else axis] = sign * half_size
    if not all(isfinite(value) for value in limits):
        raise ValueError("Geometry exceeds the finite numeric range.")

    opened = clip_perimeter_to_slab(closed, box(*limits))
    if not isinstance(opened, LineString) or opened.is_empty or opened.is_closed:
        raise ValueError("Could not construct a single open contour.")

    # Extend the two terminal straight legs from the face to the free edge.
    coords = [list(point) for point in opened.coords]
    low, high = sorted((edge[0][1 - axis], edge[1][1 - axis]))
    for point in (coords[0], coords[-1]):
        if not low <= point[1 - axis] <= high:
            raise ValueError("The edge must contain both contour endpoints.")
        point[axis] = edge_coordinate

    return PunchingPerimeterCandidate(
        perimeter=LineString(coords),
        reached_edges=(edge,),
    )

def corner_perimeter_candidate(
    c_x: float,
    c_y: float,
    offset: float,
    vertical_edge: Segment,
    horizontal_edge: Segment,
    *,
    corner_style: Literal["sharp", "round"] = "sharp",
    quad_segs: int = 32,
) -> PunchingPerimeterCandidate:
    """Build an open contour reaching two adjacent free slab edges.

    The rectangular column is centered at the origin. The edges must
    share a corner endpoint and lie outside or on the column faces.
    Each edge must contain its contour endpoint. All lengths use the
    same unit. Other boundaries, openings and selection are not handled.
    """
    closed = column_offset_perimeter(
        c_x, c_y, offset, corner_style=corner_style, quad_segs=quad_segs
    )
    edges = (vertical_edge, horizontal_edge)
    faces = []

    for axis, edge in enumerate(edges):
        if len(edge) != 2 or any(len(point) != 2 for point in edge):
            raise ValueError("Each edge needs two 2D endpoints.")
        if not all(isfinite(value) for point in edge for value in point):
            raise ValueError("Edge coordinates must be finite.")
        if edge[0][axis] != edge[1][axis] or edge[0][1 - axis] == edge[1][1 - axis]:
            raise ValueError("Provide a nonzero vertical edge and a horizontal edge.")

        coordinate = edge[0][axis]
        half_size = (c_x, c_y)[axis] / 2
        if abs(coordinate) < half_size:
            raise ValueError("Edges cannot cross the column interior.")
        faces.append(half_size if coordinate > 0 else -half_size)

    corner = (vertical_edge[0][0], horizontal_edge[0][1])
    if corner not in vertical_edge or corner not in horizontal_edge:
        raise ValueError("The edges must share their corner endpoint.")

    # Open the offset at both column faces, before extending its ends.
    xmin, ymin, xmax, ymax = closed.bounds
    margin = max(c_x, c_y, offset)
    limits = [xmin - margin, ymin - margin, xmax + margin, ymax + margin]
    for axis, face in enumerate(faces):
        limits[axis + 2 if face > 0 else axis] = face
    if not all(isfinite(value) for value in limits):
        raise ValueError("Geometry exceeds the finite numeric range.")

    opened = clip_perimeter_to_slab(closed, box(*limits))
    if not isinstance(opened, LineString) or opened.is_empty or opened.is_closed:
        raise ValueError("Could not construct a single open contour.")

    # Each end reaches a different edge; preserve the opposite corner.
    coords = [list(point) for point in opened.coords]
    for axis, edge in enumerate(edges):
        ends = [i for i in (0, -1) if coords[i][axis] == faces[axis]]
        if len(ends) != 1:
            raise ValueError("Could not identify the end facing each edge.")
        point = coords[ends[0]]
        low, high = sorted((edge[0][1 - axis], edge[1][1 - axis]))
        if not low <= point[1 - axis] <= high:
            raise ValueError("Each edge must contain its contour endpoint.")
        point[axis] = edge[0][axis]

    return PunchingPerimeterCandidate(
        perimeter=LineString(coords),
        reached_edges=edges,
    )

def candidate_fits_slab(
    candidate: PunchingPerimeterCandidate,
    slab_outline: Polygon,
) -> bool:
    """Check whether a complete candidate fits the actual slab.

    The contour must remain inside, with no length along a free edge.
    Recorded reached edges must belong to the slab boundary.
    All exterior slab edges are assumed free; openings are separate.
    This is a geometric check, not a design-code admissibility check.
    """
    perimeter = candidate.perimeter
    clipped = clip_perimeter_to_slab(perimeter, slab_outline)

    if perimeter.is_empty or not perimeter.equals(clipped):
        return False

    return all(
        slab_outline.boundary.covers(LineString(edge))
        for edge in candidate.reached_edges
    )

def rectangular_slab_perimeter_candidates(
    c_x: float,
    c_y: float,
    offset: float,
    slab_outline: Polygon,
    *,
    corner_style: Literal["sharp", "round"] = "sharp",
    quad_segs: int = 32,
    max_length: float | None = None,
) -> List[PunchingPerimeterCandidate]:
    """Generate closed, one-edge and adjacent-edge punching candidates.

    Both column and slab are axis-aligned rectangles; the column is
    centered at the origin and fully inside/on the slab. All slab edges
    are free, with no holes. Lengths use one common unit.
    No length limit is inferred. This does not select a governing section.
    Opposite-edge and three/four-edge configurations are not generated.
    An empty list means no supported candidate passed, not a safe design.
    """
    closed = column_offset_perimeter(
        c_x, c_y, offset, corner_style=corner_style, quad_segs=quad_segs
    )
    interior = PunchingPerimeterCandidate(perimeter=closed)
    interior_fits = candidate_fits_slab(interior, slab_outline)

    xmin, ymin, xmax, ymax = slab_outline.bounds
    if not slab_outline.equals(box(xmin, ymin, xmax, ymax)):
        raise ValueError("This generator requires an axis-aligned rectangular slab.")
    if not slab_outline.covers(box(-c_x / 2, -c_y / 2, c_x / 2, c_y / 2)):
        raise ValueError("The slab must contain the complete column footprint.")

    left: Segment = ((xmin, ymin), (xmin, ymax))
    right: Segment = ((xmax, ymin), (xmax, ymax))
    bottom: Segment = ((xmin, ymin), (xmax, ymin))
    top: Segment = ((xmin, ymax), (xmax, ymax))

    combinations: List[Tuple[Segment | None, Segment | None]] = [(None, None)]
    combinations.extend((edge, None) for edge in (left, right))
    combinations.extend((None, edge) for edge in (bottom, top))
    combinations.extend((v, h) for v in (left, right) for h in (bottom, top))

    candidates: List[PunchingPerimeterCandidate] = []
    for vertical, horizontal in combinations:
        gap_x = None if vertical is None else abs(vertical[0][0]) - c_x / 2
        gap_y = None if horizontal is None else abs(horizontal[0][1]) - c_y / 2
        if not candidate_fits_length_limit(
            c_x, c_y, offset,
            max_length=max_length,
            face_gap_x=gap_x,
            face_gap_y=gap_y,
            corner_style=corner_style,
            quad_segs=quad_segs,
        ):
            continue

        # Check expected bounds before building: finite edges must reach the ends.
        bounds = list(closed.bounds)
        for axis, edge in enumerate((vertical, horizontal)):
            if edge is not None:
                coordinate = edge[0][axis]
                bounds[axis + 2 if coordinate > 0 else axis] = coordinate
        if bounds[0] < xmin or bounds[1] < ymin or bounds[2] > xmax or bounds[3] > ymax:
            continue

        if vertical is None and horizontal is None:
            if interior_fits:
                candidates.append(interior)
            continue
        elif vertical is not None and horizontal is not None:
            candidate = corner_perimeter_candidate(
                c_x, c_y, offset, vertical, horizontal,
                corner_style=corner_style, quad_segs=quad_segs,
            )
        else:
            edge = vertical if vertical is not None else horizontal
            assert edge is not None
            candidate = edge_perimeter_candidate(
                c_x, c_y, offset, edge,
                corner_style=corner_style, quad_segs=quad_segs,
            )

        if candidate_fits_slab(candidate, slab_outline):
            candidates.append(candidate)

    return candidates

def rectangular_opening_shadow_angles(
    x: float,
    y: float,
    b: float,
    h: float,
    *,
    rotation: float = 0.0,
) -> Tuple[float, float]:
    """Return the tangent angular interval of a rectangular opening.

    The center (x, y) is measured from the column center in global axes.
    Dimensions b and h follow the opening's local axes. All lengths must
    use the same unit. Rotation is in radians, counterclockwise about
    the opening center; zero keeps its sides aligned with the global axes.

    Return (start, end) in global radians, counterclockwise from +x,
    with 0 < end - start < pi. Angles may lie outside [-pi, pi].
    Containment or contact with the origin, including numerically
    indistinguishable contact, is rejected. This does not check overlap
    with the column footprint or design-code proximity limits.
    """
    if not all(isfinite(value) for value in (x, y, b, h, rotation)):
        raise ValueError("Opening coordinates, dimensions and rotation must be finite.")
    if b <= 0 or h <= 0:
        raise ValueError("Opening dimensions must be positive.")

    # Express the opening center in its local axes, keeping the origin
    # at the column center. The rectangle is axis-aligned in this frame.
    c, s = cos(rotation), sin(rotation)
    local_x = x * c + y * s
    local_y = -x * s + y * c
    xmin, xmax = local_x - b / 2, local_x + b / 2
    ymin, ymax = local_y - h / 2, local_y + h / 2

    if (
        not all(isfinite(value) for value in (xmin, xmax, ymin, ymax))
        or not xmin < xmax
        or not ymin < ymax
    ):
        raise ValueError("Opening bounds cannot be represented at this numeric scale.")

    # Allow for floating-point roundoff when testing contact with the origin.
    tolerance = 8 * ulp(max(abs(x), abs(y), b, h))
    if abs(local_x) <= b / 2 + tolerance and abs(local_y) <= h / 2 + tolerance:
        raise ValueError("The opening contains or touches the column center within numeric precision.")

    local_reference = atan2(local_y, local_x)
    relative_angles = [
        (atan2(py, px) - local_reference + pi) % (2 * pi) - pi
        for px in (xmin, xmax)
        for py in (ymin, ymax)
    ]
    global_reference = atan2(y, x)
    start = global_reference + min(relative_angles)
    end = global_reference + max(relative_angles)
    if not 0 < end - start < pi:
        raise ValueError("Opening tangents cannot be resolved at this numeric scale.")
    return start, end

def circular_opening_shadow_angles(
    x: float,
    y: float,
    radius: float,
) -> Tuple[float, float]:
    """Return the tangent angular interval of a circular opening.

    The center (x, y) is measured from the column center in global axes.
    All lengths must use the same unit. Return (start, end) in radians,
    counterclockwise from +x, with 0 < end - start < pi. Angles may lie
    outside [-pi, pi]. Rotation does not affect a circular opening.

    Containment or contact with the origin, including numerically
    indistinguishable contact, is rejected. This does not check overlap
    with the column footprint or design-code proximity limits.
    """
    if not all(isfinite(value) for value in (x, y, radius)):
        raise ValueError("Opening coordinates and radius must be finite.")
    if radius <= 0:
        raise ValueError("Opening radius must be positive.")

    distance = hypot(x, y)
    if not isfinite(distance):
        raise ValueError("Opening distance cannot be represented at this numeric scale.")

    tolerance = 8 * ulp(max(distance, radius))
    if distance - radius <= tolerance:
        raise ValueError("The opening contains or touches the column center within numeric precision.")

    center_angle = atan2(y, x)
    half_angle = asin(radius / distance)
    start = center_angle - half_angle
    end = center_angle + half_angle
    if not 0 < end - start < pi:
        raise ValueError("Opening tangents cannot be resolved at this numeric scale.")
    return start, end

def subtract_opening_shadow(
    perimeter: LineString | MultiLineString,
    start: float,
    end: float,
) -> LineString | MultiLineString:
    """Remove the angular shadow of an opening from a perimeter.

    The shadow starts at the column center (0, 0). Angles are in radians,
    measured counterclockwise from global +x, with 0 < end - start < pi.
    They may lie outside [-pi, pi]. The input geometry must be valid and 2D.

    Return a new linear geometry, possibly disconnected or empty. Do not
    close the remaining paths or alter the original perimeter. Opening
    proximity limits and design-code applicability are checked elsewhere.
    """
    if not isinstance(perimeter, (LineString, MultiLineString)):
        raise TypeError("Perimeter must be a LineString or MultiLineString.")
    if perimeter.has_z or perimeter.has_m:
        raise ValueError("Perimeter coordinates must be two-dimensional.")
    if not perimeter.is_valid:
        raise ValueError("Perimeter geometry must be valid.")
    if not all(isfinite(angle) for angle in (start, end)):
        raise ValueError("Shadow angles must be finite.")

    width = end - start
    if not 0 < width < pi:
        raise ValueError("The shadow must have an angular width between zero and pi.")
    if perimeter.is_empty:
        return LineString()

    # Enclose the entire perimeter within a circle centered at the origin.
    xmin, ymin, xmax, ymax = perimeter.bounds
    radius = hypot(max(abs(xmin), abs(xmax)), max(abs(ymin), abs(ymax)))
    reach = 2 * radius
    if not isfinite(reach) or reach <= 0:
        raise ValueError("Shadow extent cannot be represented at this numeric scale.")

    # The middle point keeps both outer chords beyond the enclosing circle,
    # even when the angular width approaches pi.
    middle = start + width / 2
    shadow = Polygon([
        (0.0, 0.0),
        (reach * cos(start), reach * sin(start)),
        (reach * cos(middle), reach * sin(middle)),
        (reach * cos(end), reach * sin(end)),
    ])
    if not shadow.is_valid:
        raise ValueError("The shadow polygon could not be constructed.")

    reduced = perimeter.difference(shadow)
    merged = line_merge(reduced)
    return LineString() if merged.is_empty else merged