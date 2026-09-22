from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, Optional, Tuple

from mento.units import Quantity

from mento.codes.registry import design_code
from mento.column import Column
from mento.material import Concrete, SteelBar
from mento.plots import plotting_import_error
from mento.punching_results import PunchingCheck, PunchingCheckNotRunError, envelope_punching
from mento.reports import punching as punching_reports
from mento.units import mm, cm, inch

if TYPE_CHECKING:
    from mento.forces import Forces


#: Bar assumed in each direction when the user declares no reinforcement.
#: ``d_avg = h - c_c - d_b`` falls out of a two-mat layout of this diameter,
#: which is the estimate the class has always used.
_BAR_ESTIMATE_METRIC = 16 * mm
_BAR_ESTIMATE_IMPERIAL = 0.625 * inch  # ~#5 bar


@dataclass
class PunchingSlab:
    """
    Two-way slab at a column, for punching shear analysis.

    The top (tension) reinforcement over the column is declared as bars, and the
    effective depths and reinforcement ratios are *derived* from them::

        slab = PunchingSlab(concrete=conc, steel_bar=steel, h=25*cm, c_c=25*mm)
        slab.set_rebar_x(d_b1=12*mm, s_b1=15*cm,    # base mat, x
                         d_b3=16*mm, s_b3=15*cm)    # extra bars over the column
        slab.set_rebar_y(d_b1=12*mm, s_b1=15*cm)
        slab.d_x, slab.d_y, slab.d_avg
        slab.rho_x, slab.rho_y, slab.rho_l

    Position 1 is the base mat and position 3 a second, interleaved set, the
    same convention as ``OneWaySlab.set_slab_longitudinal_rebar_top``.

    ρ and *d* are not independent — ρ = A_s/(b·d) — so neither can be set on its
    own. ``rho_x`` / ``rho_y`` / ``rho_l`` and ``d_x`` / ``d_y`` / ``d_avg`` are
    read-only; :meth:`set_effective_depth` is the escape hatch for a layout the
    bar declaration cannot express, and it re-derives ρ from the depths it sets.

    Parameters
    ----------
    concrete : Concrete
    steel_bar : SteelBar
    h : Quantity — slab thickness
    c_c : Quantity — clear cover to the outermost bar
    outer_direction : "x" | "y" — which top mat sits closest to the slab face.
        It moves each effective depth by one bar diameter, so it is a real
        detailing input rather than a convention.

    Notes
    -----
    With no reinforcement declared, the effective depths fall back to a two-mat
    Ø16 (metric) / #5 (imperial) layout, which reproduces the historical
    ``d_avg = h - c_c - 16 mm`` estimate. That is enough for an ACI 318-19
    punching check, whose ``v_c`` does not use ρ at all; an EN 1992 check needs
    ρ and has none, so ``rho_x`` / ``rho_y`` are ``None`` rather than zero —
    ρ = 0 would silently collapse ``v_Rd,c`` to ``v_min``.
    """

    concrete: Concrete
    steel_bar: SteelBar
    h: Quantity
    c_c: Quantity
    outer_direction: Literal["x", "y"] = "x"

    # Derived; never passed in.
    _d_x: Quantity = field(init=False, repr=False)
    _d_y: Quantity = field(init=False, repr=False)

    def __post_init__(self) -> None:
        # Slab thickness must be a physical length.
        if not isinstance(self.h, Quantity) or not self.h.check("[length]"):
            raise TypeError("h must be a length Quantity.")

        h_mm = self.h.to("mm").magnitude
        if not math.isfinite(h_mm) or h_mm <= 0:
            raise ValueError("h must be greater than zero.")

        # Cover must be a physical length.
        if not isinstance(self.c_c, Quantity) or not self.c_c.check("[length]"):
            raise TypeError("c_c must be a length Quantity.")

        cover_mm = self.c_c.to("mm").magnitude
        if not math.isfinite(cover_mm) or cover_mm < 0:
            raise ValueError("c_c must be greater than or equal to zero.")

        if self.outer_direction not in ("x", "y"):
            raise ValueError(f"outer_direction must be 'x' or 'y', got {self.outer_direction!r}")

        zero = self._zero_diameter()
        # Two interleaved bar sets per direction: position 1 (base mat) and
        # position 3 (the extra bars over the column).
        self._d_b1_x, self._s_b1_x = zero, 0 * cm
        self._d_b3_x, self._s_b3_x = zero, 0 * cm
        self._d_b1_y, self._s_b1_y = zero, 0 * cm
        self._d_b3_y, self._s_b3_y = zero, 0 * cm

        # Set by set_effective_depth; None means "derive from the bars".
        self._d_x_override: Optional[Quantity] = None
        self._d_y_override: Optional[Quantity] = None

        # Filled by the Markdown view, as on every other element.
        self._md_data: str = ""

        self._update_effective_depths()

    # ------------------------------------------------------------------
    # Reinforcement
    # ------------------------------------------------------------------

    def set_rebar_x(
        self,
        d_b1: Quantity = 0 * mm,
        s_b1: Quantity = 0 * cm,
        d_b3: Quantity = 0 * mm,
        s_b3: Quantity = 0 * cm,
    ) -> None:
        """Declare the top reinforcement spanning in x, and re-derive *d* and ρ.

        Args:
            d_b1: Diameter of the position 1 bars — the base mat.
            s_b1: Spacing of the position 1 bars.
            d_b3: Diameter of the position 3 bars — the extra bars over the column.
            s_b3: Spacing of the position 3 bars.

        The call replaces whatever was declared in x, so omitting position 3
        clears it rather than keeping a previous value.
        """
        self._set_rebar("x", d_b1, s_b1, d_b3, s_b3)

    def set_rebar_y(
        self,
        d_b1: Quantity = 0 * mm,
        s_b1: Quantity = 0 * cm,
        d_b3: Quantity = 0 * mm,
        s_b3: Quantity = 0 * cm,
    ) -> None:
        """Declare the top reinforcement spanning in y, and re-derive *d* and ρ.

        See :meth:`set_rebar_x`; the two directions are independent.
        """
        self._set_rebar("y", d_b1, s_b1, d_b3, s_b3)

    def set_effective_depth(self, d_x: Optional[Quantity] = None, d_y: Optional[Quantity] = None) -> None:
        """Override the derived effective depths, one direction or both.

        For a layout the bar declaration cannot express — a mat at a varying
        depth, a depth read off a drawing. ρ is re-derived from the declared
        bars against the depth given here, so the two stay consistent, which is
        exactly what assigning to ``d_avg`` used to break.

        Passing ``None`` for a direction leaves it derived from the bars.
        """
        if d_x is None and d_y is None:
            raise ValueError("set_effective_depth needs d_x, d_y, or both.")
        for name, value in (("d_x", d_x), ("d_y", d_y)):
            if value is None:
                continue
            if not isinstance(value, Quantity) or not value.check("[length]"):
                raise TypeError(f"{name} must be a length Quantity.")
            if value.to("mm").magnitude <= 0:
                raise ValueError(f"{name} must be greater than zero.")
        if d_x is not None:
            self._d_x_override = d_x
        if d_y is not None:
            self._d_y_override = d_y
        self._update_effective_depths()

    def _set_rebar(
        self,
        direction: str,
        d_b1: Quantity,
        s_b1: Quantity,
        d_b3: Quantity,
        s_b3: Quantity,
    ) -> None:
        for position, d_b, s_b in ((1, d_b1, s_b1), (3, d_b3, s_b3)):
            self._validate_bar_set(direction, position, d_b, s_b)
        setattr(self, f"_d_b1_{direction}", d_b1)
        setattr(self, f"_s_b1_{direction}", s_b1)
        setattr(self, f"_d_b3_{direction}", d_b3)
        setattr(self, f"_s_b3_{direction}", s_b3)
        self._update_effective_depths()

    @staticmethod
    def _validate_bar_set(direction: str, position: int, d_b: Quantity, s_b: Quantity) -> None:
        """A diameter and a spacing are meaningless without each other."""
        label = f"position {position} in {direction}"
        for name, value in ((f"d_b{position}", d_b), (f"s_b{position}", s_b)):
            if not isinstance(value, Quantity) or not value.check("[length]"):
                raise TypeError(f"{name} must be a length Quantity.")
            if value.magnitude < 0:
                raise ValueError(f"{name} must be greater than or equal to zero.")
        has_bar = d_b.magnitude > 0
        has_spacing = s_b.magnitude > 0
        if has_bar != has_spacing:
            missing = "a spacing" if has_bar else "a diameter"
            raise ValueError(f"The bars at {label} need {missing} as well.")

    # ------------------------------------------------------------------
    # Derived geometry
    # ------------------------------------------------------------------

    def _zero_diameter(self) -> Quantity:
        return 0 * mm if self.concrete.unit_system == "metric" else 0 * inch

    def _bar_estimate(self) -> Quantity:
        if self.concrete.unit_system == "metric":
            return _BAR_ESTIMATE_METRIC
        return _BAR_ESTIMATE_IMPERIAL

    def _governing_diameter(self, direction: str) -> Quantity:
        """The bar that sets the depth of a mat: the largest one declared in it.

        Two interleaved sets of different diameters do not sit at one depth. The
        larger one governs, which is the conservative reading of the pair and
        the one a drawing dimensions to.
        """
        d_b1 = getattr(self, f"_d_b1_{direction}")
        d_b3 = getattr(self, f"_d_b3_{direction}")
        largest = max((d_b1, d_b3), key=lambda d_b: d_b.to("mm").magnitude)
        if largest.magnitude <= 0:
            return self._bar_estimate()
        return largest

    def _update_effective_depths(self) -> None:
        """Stack the two mats under the cover and record where each one sits."""
        d_b_x = self._governing_diameter("x")
        d_b_y = self._governing_diameter("y")
        free: Quantity = self.h - self.c_c

        if self.outer_direction == "x":
            d_x: Quantity = free - d_b_x / 2
            d_y: Quantity = free - d_b_x - d_b_y / 2
        else:
            d_y = free - d_b_y / 2
            d_x = free - d_b_y - d_b_x / 2

        self._d_x = self._d_x_override if self._d_x_override is not None else d_x
        self._d_y = self._d_y_override if self._d_y_override is not None else d_y

        if self._d_x.to("mm").magnitude <= 0 or self._d_y.to("mm").magnitude <= 0:
            raise ValueError("d_avg must be greater than zero: h is too small for the given c_c and bar layout.")

    @property
    def d_x(self) -> Quantity:
        """Effective depth of the reinforcement spanning in x."""
        return self._d_x

    @property
    def d_y(self) -> Quantity:
        """Effective depth of the reinforcement spanning in y."""
        return self._d_y

    @property
    def d_avg(self) -> Quantity:
        """Mean effective depth, ``(d_x + d_y)/2`` — EN 1992-1-1 eq. (6.32)."""
        return (self._d_x + self._d_y) / 2

    @d_avg.setter
    def d_avg(self, value: Quantity) -> None:
        raise AttributeError(
            "d_avg is derived from d_x and d_y and cannot be assigned: it would leave "
            "rho_x and rho_y computed against a different depth. Declare the bars with "
            "set_rebar_x() / set_rebar_y(), or set the depths directly with "
            "set_effective_depth(d_x=..., d_y=...), which re-derives rho with them."
        )

    # ------------------------------------------------------------------
    # Derived reinforcement ratios
    # ------------------------------------------------------------------

    def _bar_sets(self, direction: str) -> Tuple[Tuple[Quantity, Quantity], ...]:
        return (
            (getattr(self, f"_d_b1_{direction}"), getattr(self, f"_s_b1_{direction}")),
            (getattr(self, f"_d_b3_{direction}"), getattr(self, f"_s_b3_{direction}")),
        )

    def _A_s(self, direction: str) -> Optional[Quantity]:
        """Top steel per unit width in one direction, base mat plus extra bars."""
        total: Quantity = 0 * mm**2 / mm
        declared = False
        for d_b, s_b in self._bar_sets(direction):
            if d_b.magnitude <= 0:
                continue
            declared = True
            total = total + (math.pi * d_b**2 / 4) / s_b
        if not declared:
            return None
        return total

    @property
    def A_s_x(self) -> Optional[Quantity]:
        """Top steel per unit width spanning in x; ``None`` if none declared."""
        return self._A_s("x")

    @property
    def A_s_y(self) -> Optional[Quantity]:
        """Top steel per unit width spanning in y; ``None`` if none declared."""
        return self._A_s("y")

    def _rho(self, direction: str) -> Optional[float]:
        A_s = self._A_s(direction)
        if A_s is None:
            return None
        d = self._d_x if direction == "x" else self._d_y
        return float((A_s / d).to("").magnitude)

    @property
    def rho_x(self) -> Optional[float]:
        """Top reinforcement ratio spanning in x; ``None`` if no bars declared.

        ``A_s,x / d_x`` per unit width. EN 1992-1-1 §6.4.4(1) averages ρ over
        the column width plus 3·d each side; the extra bars over the column are
        assumed to span that band, so the width cancels out of the ratio.
        The code's ρ ≤ 0.02 cap belongs to the check, not to the geometry, and
        is not applied here.
        """
        return self._rho("x")

    @property
    def rho_y(self) -> Optional[float]:
        """Top reinforcement ratio spanning in y; ``None`` if no bars declared."""
        return self._rho("y")

    @property
    def rho_l(self) -> Optional[float]:
        """``√(ρ_x·ρ_y)`` — EN 1992-1-1 §6.4.4(1). ``None`` unless both are known."""
        rho_x, rho_y = self.rho_x, self.rho_y
        if rho_x is None or rho_y is None:
            return None
        return math.sqrt(rho_x * rho_y)

    @property
    def has_rebar(self) -> bool:
        """Whether enough bars are declared for ρ to exist in both directions."""
        return self.rho_l is not None

    @property
    def unit_system(self) -> str:
        return self.concrete.unit_system

    # ------------------------------------------------------------------
    # Presentation — delegates, as on every other element
    # ------------------------------------------------------------------

    @property
    def data(self) -> None:
        """Show the slab's geometry, reinforcement and materials as Markdown."""
        return punching_reports.slab_data(self)


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
        self._punching_checks: Tuple[PunchingCheck, ...] = ()
        self._punching_checked = False

    @property
    def id(self) -> int:
        return self._id

    def check(self) -> PunchingCheck:
        """Check every load combination and return the governing result.

        Dispatch goes through the registry, so the connection never names a
        design code and a code that has no punching check says so by name.

        The per-combination results stay on :attr:`punching_checks`; this
        returns their envelope, which for punching *is* one of them — a single
        stress against a single resistance on one perimeter, not a mixture.
        """
        if not self.forces:
            raise ValueError("check() requires at least one Forces object.")
        code = design_code(self.slab.concrete)
        checker = code.requires("check_punching")
        self._punching_checks = tuple(checker(self, force) for force in self.forces)
        self._punching_checked = True
        return envelope_punching(self._punching_checks)

    @property
    def punching_checks(self) -> Tuple[PunchingCheck, ...]:
        """One result per load combination, in the order they were given."""
        if not self._punching_checked:
            raise PunchingCheckNotRunError("No punching results yet: call check() first.")
        return self._punching_checks

    def design(self) -> None:
        """Size the punching shear reinforcement. (Phase 4.)"""
        design_code(self.slab.concrete).requires("design_punching")(self)

    def plot(self) -> None:
        """Display a plan-view of the punching node geometry.

        The drawing itself lives in :mod:`mento.plots.punching`.
        """
        try:
            from mento.plots.punching import plot_punching_node
        except ImportError as error:
            raise plotting_import_error(error) from error

        plot_punching_node(self)

    @property
    def data(self) -> None:
        """Show the node's slab, column and forces as Markdown."""
        return punching_reports.node_data(self)

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
        return f"PunchingNode(id={self._id}, {self.column!r}, forces={n_forces}{op}{cap})"
