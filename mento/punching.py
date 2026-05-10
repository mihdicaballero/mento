from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, Optional

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pandas import DataFrame

from mento.column import Column
from mento.material import Concrete, SteelBar
from mento.units import mm, cm, inch

if TYPE_CHECKING:
    from pint import Quantity
    from mento.forces import Forces


@dataclass
class PunchingSlab:
    """
    One-way or two-way slab section for punching shear analysis.

    Parameters
    ----------
    concrete : Concrete
    steel_bar : SteelBar
    h : Quantity — slab thickness
    c_c : Quantity — clear cover to outermost bar
    rho_x : float — flexural reinforcement ratio in x (used for EN 1992 and ACI eccentricity)
    rho_y : float — flexural reinforcement ratio in y
    d_avg is computed as h - c_c - 16 mm (metric) or h - c_c - 5/8 in (imperial).
    Override after construction if needed: slab.d_avg = custom_value
    """

    concrete: Concrete
    steel_bar: SteelBar
    h: Quantity
    c_c: Quantity
    rho_x: float = 0.0
    rho_y: float = 0.0
    d_avg: Quantity = field(init=False)

    def __post_init__(self) -> None:
        if self.concrete.unit_system == "metric":
            bar_estimate = 16 * mm
        else:
            bar_estimate = 0.625 * inch  # ~#5 bar
        self.d_avg = self.h - self.c_c - bar_estimate

    @property
    def unit_system(self) -> str:
        return self.concrete.unit_system


@dataclass
class Opening:
    """
    Slab opening near a column, for punching shear perimeter reduction.

    Parameters
    ----------
    shape : "rectangular" | "circular"
    x : Quantity — x-offset of opening centre from column centroid (+ = right)
    y : Quantity — y-offset of opening centre from column centroid (+ = up)
    b : Quantity — opening width in x (rectangular only)
    h : Quantity — opening height in y (rectangular only)
    diameter : Quantity — opening diameter (circular only)
    """

    shape: Literal["rectangular", "circular"]
    x: Quantity
    y: Quantity
    b: Quantity = field(default=0 * cm)
    h: Quantity = field(default=0 * cm)
    diameter: Quantity = field(default=0 * cm)

    def __post_init__(self) -> None:
        if self.shape not in ("rectangular", "circular"):
            raise ValueError(f"Opening shape must be 'rectangular' or 'circular', got {self.shape!r}")


@dataclass
class Capital:
    """
    Column capital (drop panel) for punching shear.

    Parameters
    ----------
    b : Quantity — capital width in x-direction
    h : Quantity — capital width in y-direction
    thickness : Quantity — capital depth below slab soffit
    """

    b: Quantity
    h: Quantity
    thickness: Quantity

    def __post_init__(self) -> None:
        if self.thickness.to("mm").magnitude <= 0:
            raise ValueError("Capital thickness must be positive")


def _plot_punching_node(node: PunchingNode) -> None:
    """Basic plan-view plot: column, capital, openings, free edges. No critical perimeter yet."""
    d = node.slab.d_avg.to("cm").magnitude
    col = node.column

    b_col = col.b.to("cm").magnitude
    h_col = col.h.to("cm").magnitude if col.shape == "rectangular" else b_col

    extent = max(5 * d, 3 * max(b_col, h_col), 30)

    # Slab boundary (clipped by free edges where present)
    x_slab_max = col.edge_distance_x.to("cm").magnitude if col.edge_distance_x is not None else extent
    y_slab_max = col.edge_distance_y.to("cm").magnitude if col.edge_distance_y is not None else extent
    x_slab_min = -extent
    y_slab_min = -extent

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_facecolor("white")

    # Slab background — exact slab area only, RGB(227, 227, 227)
    slab_bg = mpatches.Rectangle(
        (x_slab_min, y_slab_min),
        x_slab_max - x_slab_min,
        y_slab_max - y_slab_min,
        linewidth=0,
        facecolor="#e3e3e3",
        zorder=1,
    )
    ax.add_patch(slab_bg)

    # Free edge lines — dark blue, drawn only along the slab boundary
    _EDGE_COLOR = "#003399"
    _EDGE_LW = 2.5
    if col.position in ("edge", "corner"):
        ax.plot(
            [x_slab_max, x_slab_max],
            [y_slab_min, y_slab_max],
            color=_EDGE_COLOR,
            linewidth=_EDGE_LW,
            solid_capstyle="butt",
            label="Free edge",
            zorder=6,
        )
    if col.position == "corner":
        ax.plot(
            [x_slab_min, x_slab_max],
            [y_slab_max, y_slab_max],
            color=_EDGE_COLOR,
            linewidth=_EDGE_LW,
            solid_capstyle="butt",
            zorder=6,
        )

    # Capital — clipped to the actual slab boundary
    if node.capital is not None:
        cap_b = node.capital.b.to("cm").magnitude
        cap_h = node.capital.h.to("cm").magnitude
        cx0 = max(-cap_b / 2, x_slab_min)
        cy0 = max(-cap_h / 2, y_slab_min)
        cx1 = min(cap_b / 2, x_slab_max)
        cy1 = min(cap_h / 2, y_slab_max)
        if cx1 > cx0 and cy1 > cy0:
            cap_patch = mpatches.Rectangle(
                (cx0, cy0),
                cx1 - cx0,
                cy1 - cy0,
                linewidth=1.5,
                edgecolor="#555555",
                facecolor="#c0c0c0",
                linestyle="--",
                zorder=3,
                label="Capital",
            )
            ax.add_patch(cap_patch)

    # Column — medium gray (lighter than old #404040, darker than slab)
    if col.shape == "rectangular":
        col_patch = mpatches.Rectangle(
            (-b_col / 2, -h_col / 2),
            b_col,
            h_col,
            linewidth=2,
            edgecolor="black",
            facecolor="#808080",
            zorder=4,
        )
    else:
        col_patch = mpatches.Circle(
            (0, 0),
            b_col / 2,
            linewidth=2,
            edgecolor="black",
            facecolor="#808080",
            zorder=4,
        )
    ax.add_patch(col_patch)

    # Openings — white fill + border + diagonal X (no hatch)
    _OP_COLOR = "#009431"
    for i, opening in enumerate(node.openings):
        ox = opening.x.to("cm").magnitude
        oy = opening.y.to("cm").magnitude
        label = "Opening" if i == 0 else None

        if opening.shape == "rectangular":
            ob = opening.b.to("cm").magnitude
            oh = opening.h.to("cm").magnitude
            ax.add_patch(
                mpatches.Rectangle(
                    (ox - ob / 2, oy - oh / 2),
                    ob,
                    oh,
                    linewidth=1.5,
                    edgecolor=_OP_COLOR,
                    facecolor="white",
                    zorder=5,
                    label=label,
                )
            )
            # X from corner to corner
            ax.plot([ox - ob / 2, ox + ob / 2], [oy - oh / 2, oy + oh / 2], color=_OP_COLOR, linewidth=1.2, zorder=6)
            ax.plot([ox - ob / 2, ox + ob / 2], [oy + oh / 2, oy - oh / 2], color=_OP_COLOR, linewidth=1.2, zorder=6)
        else:
            od = opening.diameter.to("cm").magnitude
            ax.add_patch(
                mpatches.Circle(
                    (ox, oy),
                    od / 2,
                    linewidth=1.5,
                    edgecolor=_OP_COLOR,
                    facecolor="white",
                    zorder=5,
                    label=label,
                )
            )
            # X spanning the bounding box of the circle
            r = od / 2
            ax.plot([ox - r, ox + r], [oy - r, oy + r], color=_OP_COLOR, linewidth=1.2, zorder=6)
            ax.plot([ox - r, ox + r], [oy + r, oy - r], color=_OP_COLOR, linewidth=1.2, zorder=6)

    ax.set_xlim(-extent * 1.1, extent * 1.1)
    ax.set_ylim(-extent * 1.1, extent * 1.1)
    ax.set_aspect("equal")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.set_title(f"Punching Node — {col.position.capitalize()} {col.shape.capitalize()} Column")
    ax.grid(True, alpha=0.3, zorder=0)

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    plt.show()
    plt.close("all")


class PunchingNode:
    """
    Associates a PunchingSlab with a Column and one or more Forces objects.

    Parameters
    ----------
    slab : PunchingSlab
    column : Column
    forces : Forces or list[Forces]
    openings : list[Opening], optional
    capital : Capital, optional
    """

    _last_id: int = 0

    def __init__(
        self,
        slab: PunchingSlab,
        column: Column,
        forces: Forces | list[Forces],
        openings: Optional[list[Opening]] = None,
        capital: Optional[Capital] = None,
    ) -> None:
        PunchingNode._last_id += 1
        self._id = PunchingNode._last_id
        self.slab = slab
        self.column = column
        self.forces = forces if isinstance(forces, list) else [forces]
        self.openings = openings if openings is not None else []
        self.capital = capital

    @property
    def id(self) -> int:
        return self._id

    def check(self) -> DataFrame:
        """Run punching shear check for all forces. (Available from Phase 2.)"""
        raise NotImplementedError("Punching check not yet implemented — coming in Phase 2 (ACI) / Phase 6 (EN 1992).")

    def design(self) -> DataFrame:
        """Design punching shear reinforcement. (Available from Phase 4.)"""
        raise NotImplementedError("Punching design not yet implemented — coming in Phase 4.")

    def plot(self) -> None:
        """Display a plan-view of the punching node geometry."""
        _plot_punching_node(self)

    def __repr__(self) -> str:
        n_forces = len(self.forces)
        n_openings = len(self.openings)
        if self.capital:
            b_cm = self.capital.b.to("cm").magnitude
            h_cm = self.capital.h.to("cm").magnitude
            cap = f", capital={b_cm:.1f}×{h_cm:.1f} cm"
        else:
            cap = ""
        op = f", openings={n_openings}" if n_openings else ""
        return f"PunchingNode(id={self._id}, {self.column!r}, " f"forces={n_forces}{op}{cap})"
