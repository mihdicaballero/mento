"""Drawing of a beam cross-section.

The section drawing used to live on ``RectangularBeam`` itself, which made the
element class part matplotlib. Phase 3 of the architecture roadmap moves it
here; ``beam.plot()`` stays as a one-line delegation, so nothing calling it
changes.

These are module functions taking the beam, the same shape the design-code
modules use. They still write ``_fig`` and ``_ax`` back onto it, because that
is what the notebook views and the Word reports pick the figure up from.

A beam is drawn from its public :class:`~mento.section_geometry.SectionGeometry`
-- the legs, stirrups and bars where the checks assume them -- by helpers that
take the axes and the geometry, not the beam. A slab strip, which the geometry
gives no bars, keeps the drawing it always had.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Dict, List, Sequence, Tuple, cast

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle
from mento.units import Quantity

from mento.design_results import GRID, DesignNotRunError, format_transverse_rebar, placed_bars
from mento.results import CUSTOM_COLORS
from mento.section_geometry import BarPosition, Crosstie, SectionGeometry

if TYPE_CHECKING:
    from matplotlib.axes import Axes

    from mento.beam import RectangularBeam
    from mento.settings import BeamSettings


def _axes(beam: "RectangularBeam") -> "Axes":
    """The axes currently being drawn on.

    Typed ``Optional`` on the section because a section that was never plotted
    has none; every helper here runs inside :func:`plot_beam_section`, which
    creates them first.
    """
    return cast("Axes", beam._ax)


def _settings(beam: "RectangularBeam") -> "BeamSettings":
    """The beam's settings, which ``__post_init__`` always fills in."""
    return cast("BeamSettings", beam.settings)


def _plot_rebar_layer(
    self: "RectangularBeam",
    width_cm: float,
    height_cm: float,
    c_c_cm: float,
    stirrup_d_b_cm: float,
    layers_spacing_cm: float,
    n1: float,
    d_b1: Quantity,
    n2: float,
    d_b2: Quantity,
    max_db: Quantity,
    is_bottom: bool = True,
    is_second_layer: bool = False,
) -> None:
    """
    Helper method to plot a single layer of rebars.
    """
    # A slab strip carries width / s bars, which need not be a whole number
    # (mento.slab._bars_at_spacing); what is drawn is the whole bars that
    # cover the strip. A beam's count is whole already.
    n1, n2 = placed_bars(n1), placed_bars(n2)

    # Calculate y-position based on layer and bottom/top
    y_base = c_c_cm + stirrup_d_b_cm if is_bottom else height_cm - c_c_cm - stirrup_d_b_cm

    if is_second_layer:
        y_base += (
            layers_spacing_cm + max_db.to("cm").magnitude
            if is_bottom
            else -layers_spacing_cm - max_db.to("cm").magnitude
        )

    # Plot side bars (position 1 or 3)
    if n1 > 0:
        diameter_cm = d_b1.to("cm").magnitude
        radius_cm = diameter_cm / 2.0

        # nominal vertical center before corner correction
        y_center_nominal = y_base + radius_cm if is_bottom else y_base - radius_cm

        # width available between inner faces of stirrup legs, for this bar diameter
        clear_span = width_cm - 2 * (c_c_cm + stirrup_d_b_cm + radius_cm)

        # corner offset depends on stirrup diameter (controls bend radius)
        corner_offset = 0.43 * stirrup_d_b_cm  # tune factor if needed

        for i in range(n1):
            # even if n1 == 1, just center it
            if n1 == 1:
                x_nominal = width_cm / 2.0
            else:
                x_nominal = c_c_cm + stirrup_d_b_cm + radius_cm + i * (clear_span / (n1 - 1))

            # default: no shift
            x_shift = 0.0
            y_shift = 0.0

            # leftmost bar
            if i == 0:
                x_shift = corner_offset  # push inward (to the right)
                y_shift = corner_offset if is_bottom else -corner_offset
            # rightmost bar
            elif i == n1 - 1:
                x_shift = -corner_offset  # push inward (to the left)
                y_shift = corner_offset if is_bottom else -corner_offset

            x_plot = x_nominal + x_shift
            y_plot = y_center_nominal + y_shift

            circle = Circle(
                (x_plot, y_plot),
                radius_cm,
                color=CUSTOM_COLORS["dark_gray"],
                fill=True,
            )
            _axes(self).add_patch(circle)

    # ---------------------------------
    # Plot intermediate bars (group n2)
    # ---------------------------------
    if n2 > 0:
        diameter_cm = d_b2.to("cm").magnitude
        radius_cm = diameter_cm / 2.0

        y_center_nominal = y_base + radius_cm if is_bottom else y_base - radius_cm

        clear_span = width_cm - 2 * (c_c_cm + stirrup_d_b_cm + radius_cm)

        for i in range(n2):
            x_nominal = c_c_cm + stirrup_d_b_cm + radius_cm + (i + 1) * (clear_span / (n2 + 1))

            # intermediate bars: no special offset
            x_plot = x_nominal
            y_plot = y_center_nominal

            circle = Circle(
                (x_plot, y_plot),
                radius_cm,
                color=CUSTOM_COLORS["dark_gray"],
                fill=True,
            )
            _axes(self).add_patch(circle)


def _format_rebar_layer_text(
    self: "RectangularBeam",
    n1: int,
    d_b1: Quantity,
    n2: int,
    d_b2: Quantity,
) -> str:
    """
    Devuelve un string tipo '2Ø16+3Ø10' a partir de n1, d1, n2, d2.
    Si un grupo tiene n=0, no se incluye.
    Diámetros en mm.

    En slab:
        - Siempre combina en un único grupo: '5Ø12'.
    En beam:
        - Si n1 y n2 tienen el mismo diámetro, combina: '4Ø16'.
        - Si son distintos, deja el formato '2Ø16+3Ø10'.
    """

    mode = getattr(self, "mode", "beam")

    # -------------------------------
    # MODO SLAB: siempre combinar
    # -------------------------------
    if mode == "slab":
        total_bars = n1 + n2
        if total_bars == 0:
            return ""

        # Tomar el diámetro "no nulo"
        if n1 > 0 and d_b1 is not None:
            phi = d_b1.to("mm").magnitude
        elif n2 > 0 and d_b2 is not None:
            phi = d_b2.to("mm").magnitude
        else:
            return ""  # por seguridad

        return f"{total_bars}Ø{phi:.0f}"

    # -------------------------------
    # MODO BEAM
    # -------------------------------
    # Si n1 y n2 tienen el mismo diámetro y ambos > 0 → combinar
    if n1 > 0 and n2 > 0 and d_b1 is not None and d_b2 is not None:
        phi1 = d_b1.to("mm").magnitude
        phi2 = d_b2.to("mm").magnitude

        # Igualdad con una pequeña tolerancia
        if abs(phi1 - phi2) < 1e-6:
            total_bars = n1 + n2
            return f"{total_bars}Ø{phi1:.0f}"

    # Caso general: como lo tenías antes
    parts: list[str] = []

    if n1 > 0 and d_b1 is not None:
        phi1 = d_b1.to("mm").magnitude
        parts.append(f"{n1}Ø{phi1:.0f}")

    if n2 > 0 and d_b2 is not None:
        phi2 = d_b2.to("mm").magnitude
        parts.append(f"{n2}Ø{phi2:.0f}")

    return "+".join(parts) if parts else ""


def _annotate_rebar_layer_text(
    self: "RectangularBeam",
    width_cm: float,
    height_cm: float,
    c_c_cm: float,
    stirrup_d_b_cm: float,
    layers_spacing_cm: float,
    n1: float,
    d_b1: Quantity,
    n2: float,
    d_b2: Quantity,
    max_db: Quantity,
    is_bottom: bool = True,
    is_second_layer: bool = False,
) -> None:
    """
    Escribe a la derecha de la sección la leyenda de armadura para un layer.
    Ejemplo: '2Ø16+3Ø10'.
    """

    # The whole bars the drawing shows; see _plot_rebar_layer.
    text = _format_rebar_layer_text(self, placed_bars(n1), d_b1, placed_bars(n2), d_b2)
    if not text:
        return  # nada que mostrar

    # misma lógica de y_base que en _plot_rebar_layer
    y_base = c_c_cm + stirrup_d_b_cm if is_bottom else height_cm - c_c_cm - stirrup_d_b_cm

    if is_second_layer:
        shift = layers_spacing_cm + max_db.to("cm").magnitude
        y_base = y_base + shift if is_bottom else y_base - shift

    # posición vertical aproximada del centro del layer
    rep_db_cm = max_db.to("cm").magnitude
    y_center = y_base + rep_db_cm / 2.0 if is_bottom else y_base - rep_db_cm / 2.0

    # posición horizontal del texto (a la derecha de la sección)
    x_text = width_cm + 0.1 * width_cm

    _axes(self).text(
        x_text,
        y_center,
        text,
        ha="left",
        va="center",
        color=CUSTOM_COLORS["dark_gray"],
        # fontsize=10,
    )


def _annotate_stirrups_text(
    self: "RectangularBeam",
    width_cm: float,
    height_cm: float,
) -> None:
    """
    Escribe la leyenda de la grilla de una losa a la derecha, a media altura.
    Ejemplo: 'Ø10/8×15'. Una viga escribe la suya con :func:`_annotate_cage_text`.
    """
    if self._stirrup_n == 0:
        return  # nothing to show

    transverse = self.reinforcement.transverse
    # Bare magnitudes, as the drawing has always shown them: mm for the bar,
    # cm for the spacings, no unit suffix.
    text = format_transverse_rebar(
        transverse.layout,
        transverse.n_stirrups,
        f"{self._stirrup_d_b.to('mm').magnitude:.0f}",
        f"{self._stirrup_s_l.to('cm').magnitude:.0f}",
        f"{transverse.s_w.to('cm').magnitude:.0f}",
    )

    x_text = width_cm + 0.1 * width_cm
    y_text = height_cm / 2.0

    _axes(self).text(
        x_text,
        y_text,
        text,
        ha="left",
        va="center",
        color=CUSTOM_COLORS["dark_gray"],
        # fontsize=10,
    )


def _cm(value: Quantity) -> float:
    """A length of the geometry in the centimetres the drawing is made in."""
    return float(value.to("cm").magnitude)


def _add_rounded_stirrup(
    ax: "Axes",
    x0: float,
    y0: float,
    width: float,
    height: float,
    db_cm: float,
    bend_cm: float,
    facecolor: str,
) -> None:
    """
    Add one closed stirrup with rounded corners and thickness db_cm.

    All dimensions in cm. (x0, y0) is the bottom-left of the OUTER stirrup
    line, and ``bend_cm`` the inside diameter of its bends.
    """
    inner_radius = bend_cm / 2
    outer_radius = inner_radius + db_cm

    outer = FancyBboxPatch(
        (x0, y0),
        width,
        height,
        boxstyle=f"Round, pad=0, rounding_size={outer_radius}",
        edgecolor=CUSTOM_COLORS["dark_blue"],
        facecolor="white",
        linewidth=db_cm,  # thickness of the steel
    )
    ax.add_patch(outer)

    inner = FancyBboxPatch(
        (x0 + db_cm, y0 + db_cm),
        width - 2 * db_cm,
        height - 2 * db_cm,
        boxstyle=f"Round, pad=0, rounding_size={inner_radius}",
        edgecolor=CUSTOM_COLORS["dark_blue"],
        facecolor=facecolor,
        linewidth=1,
    )
    ax.add_patch(inner)


def _add_crosstie(ax: "Axes", tie: Crosstie, db_cm: float) -> None:
    """One crosstie: a leg the stirrup's thickness wide, with a short hook stub at each end.

    The leg carries ``gid="crosstie"``; the stubs, at the tie's hook angles, are
    marks of the ends rather than a detail of the bend.
    """
    x, y_bottom, y_top = _cm(tie.x), _cm(tie.y_bottom), _cm(tie.y_top)
    leg = Rectangle(
        (x - db_cm / 2, y_bottom - db_cm / 2),
        db_cm,
        y_top - y_bottom + db_cm,
        facecolor=CUSTOM_COLORS["dark_blue"],
        edgecolor=CUSTOM_COLORS["dark_blue"],
        gid="crosstie",
    )
    ax.add_patch(leg)
    stub = 4 * db_cm
    for (y_end, sign), angle in zip(((y_bottom, 1), (y_top, -1)), tie.hooks):
        radians = math.radians(180 - angle)
        ax.plot(
            [x, x + stub * math.sin(radians)],
            [y_end, y_end + sign * stub * math.cos(radians)],
            color=CUSTOM_COLORS["dark_blue"],
            linewidth=db_cm,
            gid="crosstie_hook",
        )


def _plot_stirrups_in_section(ax: "Axes", geometry: SectionGeometry) -> None:
    """Draw every closed stirrup and crosstie of the cage, at the legs the check assumes.

    Each stirrup is drawn on the centrelines of its legs and branches, so its
    outer line sits half a bar outside them: the perimeter stirrup's outer
    line is the cover. No stirrup is drawn on a section that has none.
    """
    d = _cm(geometry.stirrup_d_b)
    bend = _cm(geometry.stirrup_bend_inner_diameter)
    for stirrup in geometry.stirrups:
        x_left, x_right = _cm(stirrup.x_left), _cm(stirrup.x_right)
        y_bottom, y_top = _cm(stirrup.y_bottom), _cm(stirrup.y_top)
        _add_rounded_stirrup(
            ax,
            x0=x_left - d / 2,
            y0=y_bottom - d / 2,
            width=x_right - x_left + d,
            height=y_top - y_bottom + d,
            db_cm=d,
            bend_cm=bend,
            facecolor=CUSTOM_COLORS["light_gray"],
        )
    for tie in geometry.crossties:
        _add_crosstie(ax, tie, d)


def _plot_bars(ax: "Axes", geometry: SectionGeometry) -> None:
    """One circle per bar, where the check's clear-spacing model puts it."""
    for bar in geometry.bars:
        ax.add_patch(
            Circle(
                (_cm(bar.x), _cm(bar.y)),
                _cm(bar.d_b) / 2.0,
                color=CUSTOM_COLORS["dark_gray"],
                fill=True,
            )
        )


def _layer_text(bars: Tuple[BarPosition, ...]) -> str:
    """``2Ø16+3Ø10`` for the bars of one layer, from their groups; one group when they share a diameter."""
    groups: Dict[int, List[BarPosition]] = {}
    for bar in bars:
        groups.setdefault(bar.group, []).append(bar)
    counts = [(len(members), round(members[0].d_b.to("mm").magnitude)) for _, members in sorted(groups.items())]
    if len(counts) == 2 and counts[0][1] == counts[1][1]:
        counts = [(counts[0][0] + counts[1][0], counts[0][1])]
    return "+".join(f"{n}Ø{d:.0f}" for n, d in counts)


def _annotate_layers(ax: "Axes", geometry: SectionGeometry) -> None:
    """Write each layer's bars to the right of the section, at the height the bars are drawn."""
    x_text = 1.1 * _cm(geometry.width)
    for face in ("bottom", "top"):
        for layer in (1, 2):
            bars = geometry.bars_on(face, layer)
            if not bars:
                continue
            # The middle of the band the layer's bars occupy.
            low = min(_cm(bar.y) - _cm(bar.d_b) / 2 for bar in bars)
            high = max(_cm(bar.y) + _cm(bar.d_b) / 2 for bar in bars)
            ax.text(
                x_text,
                (low + high) / 2,
                _layer_text(bars),
                ha="left",
                va="center",
                color=CUSTOM_COLORS["dark_gray"],
            )


def _annotate_cage_text(ax: "Axes", geometry: SectionGeometry, lines: Sequence[str]) -> None:
    """The stirrup text, one artist per line, to the right of the section at mid-height.

    ``lines`` are the two lines of the notation and the arrangement of the
    cage; they are stacked around the middle of the section with a fixed
    offset in points, so the spacing reads the same at any section size.
    """
    x_text = 1.1 * _cm(geometry.width)
    y_text = _cm(geometry.height) / 2.0
    offsets = [14 * (len(lines) - 1) / 2 - 14 * i for i in range(len(lines))]
    for line, offset in zip(lines, offsets):
        ax.annotate(
            line,
            xy=(x_text, y_text),
            xytext=(0, offset),
            textcoords="offset points",
            ha="left",
            va="center",
            color=CUSTOM_COLORS["dark_gray"],
        )


def _cage_lines(self: "RectangularBeam") -> List[str]:
    """The stirrup text of a beam: its notation on two lines and the cage, in the current language.

    The notation of the last shear check, with the limit the legs were held
    to, when one has run on the section as it is; the configuration otherwise.
    """
    if self._stirrup_n == 0:
        return []
    transverse = self.reinforcement.transverse
    try:
        notation = self.shear_design.notation(separator="\n")
    except DesignNotRunError:
        notation = transverse.notation(separator="\n")
    return [*notation.split("\n"), transverse.arrangement()]


def _dimensions(self: "RectangularBeam", width_cm: float, height_cm: float) -> None:
    """The width and height of the section, with their arrows."""
    # Text and dimension offsets
    dim_offset = 2.5
    text_offset = dim_offset + 2
    # Add width dimension
    _axes(self).annotate(
        "",  # No text here, text is added separately
        xy=(0, -dim_offset),  # Start of arrow (left side)
        xytext=(width_cm, -dim_offset),  # End of arrow (right side)
        arrowprops={
            "arrowstyle": "<->",
            "lw": 1,
            "color": CUSTOM_COLORS["dark_blue"],
        },
    )
    if self.concrete.unit_system == "imperial":
        # Example: format to 2 decimal places, then use pint's compact (~P) format
        width = "{:.0f~P}".format(self.width.to("inch"))
        height = "{:.0f~P}".format(self.height.to("inch"))
    else:
        width = "{:.0f~P}".format(self.width.to("cm"))
        height = "{:.0f~P}".format(self.height.to("cm"))
    # Add width dimension text below the arrow
    _axes(self).text(
        width_cm / 2,  # Center of the arrow
        -text_offset,  # Slightly below the arrow
        width,
        ha="center",
        va="top",
        color=CUSTOM_COLORS["dark_gray"],
    )

    # Add height dimension
    _axes(self).annotate(
        "",  # No text here, text is added separately
        xy=(-dim_offset, 0),  # Start of arrow (bottom)
        xytext=(-dim_offset, height_cm),  # End of arrow (top)
        arrowprops={
            "arrowstyle": "<->",
            "lw": 1,
            "color": CUSTOM_COLORS["dark_blue"],
        },
    )
    # Add height dimension text to the left of the arrow
    _axes(self).text(
        -text_offset,  # Slightly to the left of the arrow
        height_cm / 2,  # Center of the arrow
        height,
        ha="right",
        va="center",
        color=CUSTOM_COLORS["dark_gray"],
        rotation=90,  # Rotate text vertically
    )


def _plot_grid_section(self: "RectangularBeam", width_cm: float, height_cm: float) -> None:
    """A slab strip's bars and grid label, drawn as they always were.

    A strip is detailed by spacing and has no cage, so its geometry publishes
    no bars (see :mod:`mento.section_geometry`); it keeps the whole bars that
    cover the strip and its ``Ø10/8×15`` label.
    """
    c_c_cm: float = self.c_c.to("cm").magnitude
    stirrup_d_b_cm: float = self._stirrup_d_b.to("cm").magnitude
    layers_spacing_cm: float = _settings(self).layers_spacing.to("cm").magnitude
    faces = (
        (
            True,
            (self._n1_b, self._d_b1_b, self._n2_b, self._d_b2_b),
            (self._n3_b, self._d_b3_b, self._n4_b, self._d_b4_b),
        ),
        (
            False,
            (self._n1_t, self._d_b1_t, self._n2_t, self._d_b2_t),
            (self._n3_t, self._d_b3_t, self._n4_t, self._d_b4_t),
        ),
    )
    for is_bottom, first, second in faces:
        for layer, second_layer in ((first, False), (second, True)):
            _plot_rebar_layer(
                self,
                width_cm,
                height_cm,
                c_c_cm,
                stirrup_d_b_cm,
                layers_spacing_cm,
                *layer,
                max_db=first[1],
                is_bottom=is_bottom,
                is_second_layer=second_layer,
            )
    for is_bottom, first, second in faces:
        for layer, second_layer in ((first, False), (second, True)):
            _annotate_rebar_layer_text(
                self,
                width_cm,
                height_cm,
                c_c_cm,
                stirrup_d_b_cm,
                layers_spacing_cm,
                *layer,
                max_db=first[1],
                is_bottom=is_bottom,
                is_second_layer=second_layer,
            )
    _annotate_stirrups_text(self, width_cm, height_cm)


def plot_beam_section(self: "RectangularBeam", show: bool = False) -> Figure:
    """
    Plots the rectangular section with a dark gray border, light gray hatch, and dimensions.

    A beam is drawn from its :attr:`~mento.beam.RectangularBeam.section_geometry`:
    every stirrup of the cage at the legs the shear check assumes, every bar
    where the clear-spacing model puts it, and the stirrup text in three
    lines -- the legs, bar and spacing; the spacing of the legs across the
    width with its maximum; and the arrangement of the cage. A slab strip
    keeps the drawing it always had.
    """

    # Convert dimensions to consistent units (cm)
    width_cm: float = self.width.to("cm").magnitude
    height_cm: float = self.height.to("cm").magnitude

    # Create figure and axis
    fig, self._ax = plt.subplots()
    ax = _axes(self)

    # Create a rectangle patch for the section
    rect = Rectangle(
        (0, 0),
        width_cm,
        height_cm,
        linewidth=1.3,
        edgecolor=CUSTOM_COLORS["dark_gray"],
        facecolor=CUSTOM_COLORS["light_gray"],
    )
    ax.add_patch(rect)

    geometry = self.section_geometry
    if geometry.layout != GRID:
        _plot_stirrups_in_section(ax, geometry)

    # Set plot limits with some padding
    padding = max(width_cm, height_cm) * 0.2
    ax.set_xlim(-padding, width_cm + padding)
    ax.set_ylim(-padding, height_cm + padding)

    _dimensions(self, width_cm, height_cm)

    # Set aspect of the plot to be equal
    ax.set_aspect("equal")
    # Remove axes for better visualization
    ax.axis("off")

    if geometry.layout == GRID:
        _plot_grid_section(self, width_cm, height_cm)
    else:
        _plot_bars(ax, geometry)
        _annotate_layers(ax, geometry)
        _annotate_cage_text(ax, geometry, _cage_lines(self))

    # Store the section figure
    self._fig = fig

    if show:
        plt.show()

    # # Close the figure so notebooks don't auto-display it twice
    plt.close(fig)

    return fig
