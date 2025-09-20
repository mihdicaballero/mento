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
    kgf,
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
    "kgf",
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
    "OneWaySlab",
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
    "BeamSettings",
    "BeamSummary",
]

if TYPE_CHECKING:
    from mento.beam import RectangularBeam
    from mento.codes import ACI_318_19_beam, EN_1992_2004_beam
    from mento.forces import Forces
    from mento.material import (
        Concrete_ACI_318_19,
        Concrete_CIRSOC_201_25,
        Concrete_EN_1992_2004,
        SteelBar,
    )
    from mento.settings import BeamSettings
    from mento.node import Node
    from mento.results import DocumentBuilder, Formatter, TablePrinter
    from mento.summary import BeamSummary


def __getattr__(name: str) -> object:
    # Map class names to their actual module files
    module_mapping = {
        "RectangularBeam": "beam",
        "Node": "node",
        "BeamSettings": "settings",
        "Forces": "forces",
        "Concrete_ACI_318_19": "material",
        "SteelBar": "material",
        "Concrete_CIRSOC_201_25": "material",
        "Concrete_EN_1992_2004": "material",
        "Formatter": "results",
        "TablePrinter": "results",
        "DocumentBuilder": "results",
        "EN_1992_2004_beam": "codes",
        "ACI_318_19_beam": "codes",
        "BeamSummary": "summary",
    }

    if name in module_mapping:
        import importlib

        module = importlib.import_module(f".{module_mapping[name]}", __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
