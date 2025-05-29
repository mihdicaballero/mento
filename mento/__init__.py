"""Mento structural analysis package.

This package provides tools for structural analysis and design of concrete elements.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

# Expose units and Quantity directly to the user
from .units import (
    GPa,
    MPa,
    kPa,
    kN,
    kNm,
    kg,
    kip,
    ksi,
    lb,
    m,
    mm,
    cm,
    psi,
    sec,
    ureg,
    deg,
    ft,
    inch,
)

# Re-export Quantity for user convenience
__all__ = [
    "ureg",
    "m",
    "cm",
    "mm",
    "kN",
    "kNm",
    "kPa",
    "MPa",
    "GPa",
    "kg",
    "sec",
    "psi",
    "lb",
    "kip",
    "ksi",
    "inch",
    "ft",
    "deg",
    "Node",
    "Forces",
    "Concrete_ACI_318_19",
    "SteelBar",
    "Concrete_CIRSOC_201_25",
    "Concrete_EN_1992_2004",
    "RectangularBeam",
    "Formatter",
    "TablePrinter",
    "DocumentBuilder",
    "EN_1992_2004_beam",
    "ACI_318_19_beam",
    "BeamSummary",
]

if TYPE_CHECKING:
    from .beam import RectangularBeam
    from .codes import ACI_318_19_beam, EN_1992_2004_beam
    from .forces import Forces
    from .material import (
        Concrete_ACI_318_19,
        Concrete_CIRSOC_201_25,
        Concrete_EN_1992_2004,
        SteelBar,
    )
    from .node import Node
    from .results import DocumentBuilder, Formatter, TablePrinter
    from .summary import BeamSummary

# Lazy imports to avoid circular dependencies
def __getattr__(name: str) -> object:
    if name in {
        "Node",
        "Forces",
        "Concrete_ACI_318_19",
        "SteelBar",
        "Concrete_CIRSOC_201_25",
        "Concrete_EN_1992_2004",
        "RectangularBeam",
        "Formatter",
        "TablePrinter",
        "DocumentBuilder",
        "EN_1992_2004_beam",
        "ACI_318_19_beam",
        "BeamSummary",
    }:
        import importlib

        module = importlib.import_module(f".{name.lower()}", __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")