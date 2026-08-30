from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, Optional

from pandas import DataFrame
from pint import Quantity

from mento.column import Column
from mento.material import Concrete, SteelBar
from mento.plots.punching import plot_punching_node
from mento.units import mm, cm, inch

if TYPE_CHECKING:
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

        if self.concrete.unit_system == "metric":
            bar_estimate = 16 * mm
        else:
            bar_estimate = 0.625 * inch  # ~#5 bar
        self.d_avg = self.h - self.c_c - bar_estimate

        # The cover and bar allowance must leave a usable effective depth.
        if self.d_avg.to("mm").magnitude <= 0:
            raise ValueError("d_avg must be greater than zero: h is too small for the given c_c.")

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
        """Display a plan-view of the punching node geometry.

        The drawing itself lives in :mod:`mento.plots.punching`.
        """
        plot_punching_node(self)

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
