"""Notebook views of a punching shear node.

The punching counterpart to :mod:`mento.reports.views` and
:mod:`mento.reports.walls`: turning what a user typed into something they can
read back. ``PunchingSlab.data`` and ``PunchingNode.data`` delegate here so the
elements never import IPython (ADR-0004).

Only the *inputs* are rendered for now — the check that would fill a results
view lands in Phase 2.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from IPython.display import Markdown, display
from mento.units import Quantity

if TYPE_CHECKING:
    from mento.punching import PunchingNode, PunchingSlab


def _show(markdown: str) -> None:
    """Render Markdown in a notebook; IPython ships no type information."""
    display(Markdown(markdown))  # type: ignore[no-untyped-call]


def _length(value: Quantity, imperial: bool) -> str:
    """A length in the units of its unit system, at drawing precision."""
    return f"{value.to('inch'):~P.2f}" if imperial else f"{value.to('cm'):~P.1f}"


def _diameter(value: Quantity, imperial: bool) -> str:
    """Bare bar size, as a drawing writes it — the ``Ø`` already says what it is."""
    if imperial:
        return f"{value.to('inch').magnitude:.3f}"
    return f"{value.to('mm').magnitude:.0f}"


def _area_per_width(value: Quantity, imperial: bool) -> str:
    return f"{value.to('inch**2/ft'):~P.2f}" if imperial else f"{value.to('cm**2/m'):~P.2f}"


def _rebar_label(slab: "PunchingSlab", direction: str) -> str:
    """The bars of one direction, base mat first: ``Ø12/15 cm ++ Ø16/15 cm``."""
    imperial = slab.unit_system == "imperial"
    parts = []
    for d_b, s_b in slab._bar_sets(direction):
        if d_b.magnitude <= 0:
            continue
        parts.append(f"Ø{_diameter(d_b, imperial)}/{_length(s_b, imperial)}")
    if not parts:
        return "not declared"
    return " ++ ".join(parts)


def _rho_label(rho: Optional[float]) -> str:
    return "—" if rho is None else f"{rho:.4f}"


def slab_data(self: "PunchingSlab") -> None:
    """Slab thickness, cover, derived depths, derived ratios and materials."""
    imperial = self.unit_system == "imperial"
    inner = "y" if self.outer_direction == "x" else "x"

    lines = [
        f"Punching slab, $h$={_length(self.h, imperial)}, "
        f"$c_{{c}}$={_length(self.c_c, imperial)}, "
        f"Concrete {self.concrete.name}, Rebar {self.steel_bar.name}.",
        "",
        f"Top rebar x: {_rebar_label(self, 'x')}, "
        f"$d_{{x}}$={_length(self.d_x, imperial)}, "
        f"$\\rho_{{x}}$={_rho_label(self.rho_x)}",
        "",
        f"Top rebar y: {_rebar_label(self, 'y')}, "
        f"$d_{{y}}$={_length(self.d_y, imperial)}, "
        f"$\\rho_{{y}}$={_rho_label(self.rho_y)}",
        "",
        f"$d_{{avg}}$={_length(self.d_avg, imperial)}, "
        f"$\\rho_{{l}}$={_rho_label(self.rho_l)} "
        f"({self.outer_direction} outside {inner}).",
    ]

    A_s_x, A_s_y = self.A_s_x, self.A_s_y
    if A_s_x is not None or A_s_y is not None:
        as_x = "—" if A_s_x is None else _area_per_width(A_s_x, imperial)
        as_y = "—" if A_s_y is None else _area_per_width(A_s_y, imperial)
        lines += ["", f"$A_{{s,x}}$={as_x}, $A_{{s,y}}$={as_y}."]
    else:
        lines += [
            "",
            "No reinforcement declared — the effective depths assume a two-mat "
            "estimate and ρ is unknown. Declare it with `set_rebar_x()` / `set_rebar_y()`.",
        ]

    markdown_content = "\n".join(lines)
    self._md_data = markdown_content
    _show(markdown_content)
    return None


def _forces_label(self: "PunchingNode") -> str:
    """One line per load combination: the punching load and both moments.

    ``V_z`` is the design punching load — ``V_u`` under ACI 318-19, ``V_Ed``
    under EN 1992-1-1 — and ``M_x`` / ``M_y`` are the unbalanced moments
    transferred to the slab about its two in-plane axes.
    """
    rows = []
    for force in self.forces:
        parts = [f"$V$={force.V_z:~P.1f}"]
        if force.M_x.magnitude != 0:
            parts.append(f"$M_{{x}}$={force.M_x:~P.1f}")
        if force.M_y.magnitude != 0:
            parts.append(f"$M_{{y}}$={force.M_y:~P.1f}")
        rows.append(f"- {force.label or 'unlabelled'}: " + ", ".join(parts))
    return "\n".join(rows) if rows else "- no forces"


def node_data(self: "PunchingNode") -> None:
    """The slab, the column it sits on, and the forces at the node.

    The view a user reaches for at the REPL: a slab on its own does not know
    its column, and checking the geometry you just typed means seeing both.
    """
    imperial = self.slab.unit_system == "imperial"
    column = self.column

    if column.shape == "rectangular":
        dims = f"$b$={_length(column.b, imperial)}, $h$={_length(column.h, imperial)}"
    else:
        dims = f"$D$={_length(column.b, imperial)}"

    header = f"**Punching node {self.id}** — {column.position} {column.shape} column, {dims}."

    edges = []
    if column.edge_distance_x is not None:
        edges.append(f"$x$={_length(column.edge_distance_x, imperial)}")
    if column.edge_distance_y is not None:
        edges.append(f"$y$={_length(column.edge_distance_y, imperial)}")
    if edges:
        header += " Distance to free edge: " + ", ".join(edges) + "."

    if self.capital is not None:
        cap = self.capital
        header += (
            f" Capital {_length(cap.b, imperial)} × {_length(cap.h, imperial)}, "
            f"thickness {_length(cap.thickness, imperial)}."
        )
    if self.openings:
        header += f" {len(self.openings)} opening(s)."

    _show(header)
    slab_data(self.slab)
    _show("Forces:\n\n" + _forces_label(self))
    return None
