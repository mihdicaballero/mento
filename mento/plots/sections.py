"""Drawing of a beam cross-section.

The section drawing used to live on ``RectangularBeam`` itself, which made the
element class part matplotlib. Phase 3 of the architecture roadmap moves it
here; ``beam.plot()`` stays as a one-line delegation, so nothing calling it
changes.

These are module functions taking the beam, the same shape the design-code
modules use. They still write ``_fig`` and ``_ax`` back onto it, because that
is what the notebook views and the Word reports pick the figure up from.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle
from pint import Quantity

from mento.results import CUSTOM_COLORS

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
    n1: int,
    d_b1: Quantity,
    n2: int,
    d_b2: Quantity,
    max_db: Quantity,
    is_bottom: bool = True,
    is_second_layer: bool = False,
) -> None:
    """
    Helper method to plot a single layer of rebars.
    """
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
    n1: int,
    d_b1: Quantity,
    n2: int,
    d_b2: Quantity,
    max_db: Quantity,
    is_bottom: bool = True,
    is_second_layer: bool = False,
) -> None:
    """
    Escribe a la derecha de la sección la leyenda de armadura para un layer.
    Ejemplo: '2Ø16+3Ø10'.
    """

    text = _format_rebar_layer_text(self, n1, d_b1, n2, d_b2)
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
    Escribe la leyenda de estribos a la derecha, a media altura.
    Ejemplo: '1eØ6/20' o '2eØ6/20'.
    """
    if self._stirrup_n == 0:
        return  # nothing to show

    phi_mm = self._stirrup_d_b.to("mm").magnitude
    s_cm = self._stirrup_s_l.to("cm").magnitude

    text = f"{self._stirrup_n:.0f}eØ{phi_mm:.0f}/{s_cm:.0f}"

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


def _add_rounded_stirrup(
    self: "RectangularBeam",
    x0: float,
    y0: float,
    width: float,
    height: float,
    db_cm: float,
    facecolor: str,
) -> None:
    """
    Add one closed stirrup with rounded corners and thickness db_cm.

    All dimensions in cm.
    (x0, y0) is bottom-left of the OUTER stirrup line.
    """
    # corner radii (same logic you already have)
    inner_radius = 4 * db_cm / 2
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
    _axes(self).add_patch(outer)

    inner = FancyBboxPatch(
        (x0 + db_cm, y0 + db_cm),
        width - 2 * db_cm,
        height - 2 * db_cm,
        boxstyle=f"Round, pad=0, rounding_size={inner_radius}",
        edgecolor=CUSTOM_COLORS["dark_blue"],
        facecolor=facecolor,
        linewidth=1,
    )
    _axes(self).add_patch(inner)


def _plot_stirrups_in_section(
    self: "RectangularBeam",
    c_c_cm: float,
    section_width_cm: float,
    section_height_cm: float,
    stirrup_db_cm: float,
    n_stirrups: int,
) -> None:
    """
    Plot 1, 2 or 3 stirrups inside the beam section.

    Rules:
    - 1 stirrup: full width (current behavior).
    - 2 stirrups: outer full width, inner with 1/2 width, centered.
    - 3 stirrups: outer full width, plus 2 inner stirrups whose total
    occupied width is about 1/4 of the section width (each ~1/8), placed
    symmetrically left/right.
    """

    # base geometry (outer stirrup like now)
    stirrup_width = section_width_cm - 2 * c_c_cm
    stirrup_height = section_height_cm - 2 * c_c_cm

    # Always draw the outer stirrup
    _add_rounded_stirrup(
        self,
        x0=c_c_cm,
        y0=c_c_cm,
        width=stirrup_width,
        height=stirrup_height,
        db_cm=stirrup_db_cm,
        facecolor=CUSTOM_COLORS["light_gray"],
    )

    if n_stirrups <= 1:
        return

    # Case 2 stirrups: inner one with 1/2 width, centered
    if n_stirrups == 2:
        inner_width = 0.5 * stirrup_width
        x_center = c_c_cm + stirrup_width / 2
        x0_inner = x_center - inner_width / 2

        _add_rounded_stirrup(
            self,
            x0=x0_inner,
            y0=c_c_cm,
            width=inner_width,
            height=stirrup_height,
            db_cm=stirrup_db_cm,
            facecolor=CUSTOM_COLORS["light_gray"],
        )
        return

    # Case 3+ stirrups: outer + 2 small inner ones.
    # Clamp so it doesn't exceed outer stirrup
    inner_width = min(0.25 * section_width_cm, 0.4 * stirrup_width)

    # Place inner stirrups roughly at quarter points of the outer stirrup
    # -> centers at 1/3 and 2/3 of stirrup span
    x_left_center = c_c_cm + stirrup_width * (1 / 3) - stirrup_db_cm / 2
    x_right_center = c_c_cm + stirrup_width * (2 / 3) + stirrup_db_cm / 2

    x0_left = x_left_center - inner_width / 2 - stirrup_db_cm / 2
    x0_right = x_right_center - inner_width / 2 + stirrup_db_cm / 2

    _add_rounded_stirrup(
        self,
        x0=x0_left,
        y0=c_c_cm,
        width=inner_width,
        height=stirrup_height,
        db_cm=stirrup_db_cm,
        facecolor=CUSTOM_COLORS["light_gray"],
    )
    _add_rounded_stirrup(
        self,
        x0=x0_right,
        y0=c_c_cm,
        width=inner_width,
        height=stirrup_height,
        db_cm=stirrup_db_cm,
        facecolor=CUSTOM_COLORS["light_gray"],
    )


def plot_beam_section(self: "RectangularBeam", show: bool = False) -> Figure:
    """
    Plots the rectangular section with a dark gray border, light gray hatch, and dimensions.
    Also plots the stirrup with rounded corners and thickness.
    """

    # Convert dimensions to consistent units (cm)
    width_cm: float = self.width.to("cm").magnitude
    height_cm: float = self.height.to("cm").magnitude
    c_c_cm: float = self.c_c.to("cm").magnitude
    stirrup_d_b_cm: float = self._stirrup_d_b.to("cm").magnitude
    layers_spacing_cm: float = _settings(self).layers_spacing.to("cm").magnitude

    # Create figure and axis
    fig, self._ax = plt.subplots()

    # Create a rectangle patch for the section
    rect = Rectangle(
        (0, 0),
        width_cm,
        height_cm,
        linewidth=1.3,
        edgecolor=CUSTOM_COLORS["dark_gray"],
        facecolor=CUSTOM_COLORS["light_gray"],
    )
    _axes(self).add_patch(rect)

    if self.mode == "beam":
        db = self._stirrup_d_b.to("cm").magnitude  # bar diameter Ø

        # Cap at 3 for drawing (you can show 1, 2, or 3)
        n_stirrups = max(1, min(3, self._stirrup_n))

        _plot_stirrups_in_section(
            self,
            c_c_cm=c_c_cm,
            section_width_cm=width_cm,
            section_height_cm=height_cm,
            stirrup_db_cm=db,
            n_stirrups=n_stirrups,
        )

    # Set plot limits with some padding
    padding = max(width_cm, height_cm) * 0.2
    self._ax.set_xlim(-padding, width_cm + padding)
    self._ax.set_ylim(-padding, height_cm + padding)

    # Text and dimension offsets
    dim_offset = 2.5
    text_offset = dim_offset + 2
    # Add width dimension
    self._ax.annotate(
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
    self._ax.annotate(
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

    # Set aspect of the plot to be equal
    self._ax.set_aspect("equal")
    # Remove axes for better visualization
    self._ax.axis("off")

    # Calculate rebar positions
    # Bottom rebars
    _plot_rebar_layer(
        self,
        width_cm,
        height_cm,
        c_c_cm,
        stirrup_d_b_cm,
        layers_spacing_cm,
        self._n1_b,
        self._d_b1_b,
        self._n2_b,
        self._d_b2_b,
        max_db=self._d_b1_b,
        is_bottom=True,
    )
    _plot_rebar_layer(
        self,
        width_cm,
        height_cm,
        c_c_cm,
        stirrup_d_b_cm,
        layers_spacing_cm,
        self._n3_b,
        self._d_b3_b,
        self._n4_b,
        self._d_b4_b,
        max_db=self._d_b1_b,
        is_bottom=True,
        is_second_layer=True,
    )

    # Top rebars
    _plot_rebar_layer(
        self,
        width_cm,
        height_cm,
        c_c_cm,
        stirrup_d_b_cm,
        layers_spacing_cm,
        self._n1_t,
        self._d_b1_t,
        self._n2_t,
        self._d_b2_t,
        max_db=self._d_b1_t,
        is_bottom=False,
    )
    _plot_rebar_layer(
        self,
        width_cm,
        height_cm,
        c_c_cm,
        stirrup_d_b_cm,
        layers_spacing_cm,
        self._n3_t,
        self._d_b3_t,
        self._n4_t,
        self._d_b4_t,
        max_db=self._d_b1_t,
        is_bottom=False,
        is_second_layer=True,
    )

    ### REBAR TEXT ANNOTATIONS
    # Bottom, 1st layer
    _annotate_rebar_layer_text(
        self,
        width_cm,
        height_cm,
        c_c_cm,
        stirrup_d_b_cm,
        layers_spacing_cm,
        self._n1_b,
        self._d_b1_b,
        self._n2_b,
        self._d_b2_b,
        max_db=self._d_b1_b,
        is_bottom=True,
        is_second_layer=False,
    )

    # Bottom, 2nd layer
    _annotate_rebar_layer_text(
        self,
        width_cm,
        height_cm,
        c_c_cm,
        stirrup_d_b_cm,
        layers_spacing_cm,
        self._n3_b,
        self._d_b3_b,
        self._n4_b,
        self._d_b4_b,
        max_db=self._d_b1_b,
        is_bottom=True,
        is_second_layer=True,
    )

    # Top, 1st layer
    _annotate_rebar_layer_text(
        self,
        width_cm,
        height_cm,
        c_c_cm,
        stirrup_d_b_cm,
        layers_spacing_cm,
        self._n1_t,
        self._d_b1_t,
        self._n2_t,
        self._d_b2_t,
        max_db=self._d_b1_t,
        is_bottom=False,
        is_second_layer=False,
    )

    # Top, 2nd layer
    _annotate_rebar_layer_text(
        self,
        width_cm,
        height_cm,
        c_c_cm,
        stirrup_d_b_cm,
        layers_spacing_cm,
        self._n3_t,
        self._d_b3_t,
        self._n4_t,
        self._d_b4_t,
        max_db=self._d_b1_t,
        is_bottom=False,
        is_second_layer=True,
    )

    # Stirrups text
    _annotate_stirrups_text(self, width_cm, height_cm)

    # Store the section figure
    self._fig = fig

    if show:
        plt.show()

    # # Close the figure so notebooks don't auto-display it twice
    plt.close(fig)

    return fig
