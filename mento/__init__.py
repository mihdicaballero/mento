"""Mento structural analysis package.

This package provides tools for structural analysis and design of concrete elements.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ._version import __version__

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
    "__version__",
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
    "Footing",
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
    "Column",
    "PunchingSlab",
    "Opening",
    "Capital",
    "PunchingNode",
    "ShearWall",
    "ShearWallSummary",
    "set_language",
    "get_language",
    "available_languages",
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
    from mento.slab import Footing, OneWaySlab
    from mento.settings import BeamSettings
    from mento.node import Node
    from mento.results import DocumentBuilder, Formatter, TablePrinter
    from mento.beam_summary import BeamSummary
    from mento.column import Column
    from mento.punching import Capital, Opening, PunchingNode, PunchingSlab
    from mento.shear_wall import ShearWall
    from mento.shear_wall_summary import ShearWallSummary
    from mento.i18n import available_languages, get_language, set_language


def __getattr__(name: str) -> object:
    # Map class names to their actual module files
    module_mapping = {
        "RectangularBeam": "beam",
        "OneWaySlab": "slab",
        "Footing": "slab",
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
        "BeamSummary": "beam_summary",
        "Column": "column",
        "PunchingSlab": "punching",
        "Opening": "punching",
        "Capital": "punching",
        "PunchingNode": "punching",
        "ShearWall": "shear_wall",
        "ShearWallSummary": "shear_wall_summary",
        "set_language": "i18n",
        "get_language": "i18n",
        "available_languages": "i18n",
    }

    if name in module_mapping:
        import importlib

        module = importlib.import_module(f".{module_mapping[name]}", __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
