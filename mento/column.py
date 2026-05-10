from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional, TYPE_CHECKING

from mento.units import cm

if TYPE_CHECKING:
    from pint import Quantity


@dataclass
class Column:
    """
    Column geometry for punching shear analysis.

    Parameters
    ----------
    shape : "rectangular" | "circular"
    position : "interior" | "edge" | "corner"
    b : column width in x-direction (or diameter for circular)
    h : column width in y-direction (ignored for circular)
    edge_distance_x : distance from column centroid to free slab edge in x
        Required for edge and corner columns.
    edge_distance_y : distance from column centroid to free slab edge in y
        Required for corner columns.
    """

    shape: Literal["rectangular", "circular"]
    position: Literal["interior", "edge", "corner"]
    b: Quantity = field(default=0 * cm)
    h: Quantity = field(default=0 * cm)
    edge_distance_x: Optional[Quantity] = field(default=None)
    edge_distance_y: Optional[Quantity] = field(default=None)

    def __post_init__(self) -> None:
        if self.shape not in ("rectangular", "circular"):
            raise ValueError(f"shape must be 'rectangular' or 'circular', got {self.shape!r}")
        if self.position not in ("interior", "edge", "corner"):
            raise ValueError(f"position must be 'interior', 'edge', or 'corner', got {self.position!r}")
        if self.position in ("edge", "corner") and self.edge_distance_x is None:
            raise ValueError(f"edge_distance_x is required for '{self.position}' columns")
        if self.position == "corner" and self.edge_distance_y is None:
            raise ValueError("edge_distance_y is required for corner columns")

    def __repr__(self) -> str:
        shape_str = self.shape.capitalize()
        pos_str = self.position.capitalize()
        b_cm = self.b.to("cm").magnitude
        if self.shape == "rectangular":
            h_cm = self.h.to("cm").magnitude
            dims = f"b={b_cm:.1f} cm, h={h_cm:.1f} cm"
        else:
            dims = f"d={b_cm:.1f} cm"
        return f"Column({shape_str}, {pos_str}, {dims})"
