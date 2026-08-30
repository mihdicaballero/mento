"""Plan view of a punching node.

The drawing used to sit next to :class:`~mento.punching.PunchingNode` itself,
which made the punching module part matplotlib. Phase 3 of the architecture
roadmap moves it here; ``PunchingNode.plot()`` stays as a one-line delegation,
so nothing calling it changes.

A module function taking the node, the same shape the design-code modules use.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

if TYPE_CHECKING:
    from mento.punching import PunchingNode


def plot_punching_node(self: "PunchingNode") -> None:
    """Basic plan-view plot: column, capital, openings, free edges. No critical perimeter yet."""
    d = self.slab.d_avg.to("cm").magnitude
    col = self.column

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
    if self.capital is not None:
        cap_b = self.capital.b.to("cm").magnitude
        cap_h = self.capital.h.to("cm").magnitude
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
    col_patch: mpatches.Patch
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
    for i, opening in enumerate(self.openings):
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
