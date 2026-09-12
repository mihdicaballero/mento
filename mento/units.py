from typing import TYPE_CHECKING

from pint import UnitRegistry
import pint

# Initialize UnitRegistry
ureg = UnitRegistry(system="mks")
pint.set_application_registry(ureg)
ureg.formatter.default_format = ".2f~P"  # Standardize output formatting

# Metric system units
m = ureg.meter
cm = ureg.centimeter
mm = ureg.millimeter
N = ureg.newton
kN = ureg.kilonewton
kgf = ureg.kilogram_force
kNm = kN * m
Pa = ureg.pascal
kPa = ureg.kilopascal
MPa = ureg.megapascal
GPa = ureg.gigapascal
kg = ureg.kilogram
sec = ureg.second
deg = ureg.degree
dimensionless = ureg.Quantity(1, "")

# Imperial system units
psi = ureg.psi
lb = ureg.pound
lbf = ureg.pound_force
kip = ureg.kip
ksi = ureg.ksi
inch = ureg.inch
ft = ureg.foot

# The type mento annotates its quantities with.
#
# Since pint 0.26 the stubs type every arithmetic result -- ``3 * cm``,
# ``h - c_c``, ``q.to("mm")`` -- as ``PlainQuantity``, the base class, while the
# ``pint.Quantity`` the registry hands out is a subclass of it. Annotating with
# ``pint.Quantity`` therefore rejects the result of any calculation. Under the
# type checker ``Quantity`` is the base class, which accepts both; at runtime it
# is the registry class, so ``isinstance(x, Quantity)`` keeps working and users
# can still construct one with ``Quantity(1, "m")``.
if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity as Quantity
else:
    from pint import Quantity

__all__ = [
    "Quantity",
    "ureg",
    "dimensionless",
    "m",
    "cm",
    "mm",
    "N",
    "kN",
    "kgf",
    "kNm",
    "Pa",
    "kPa",
    "MPa",
    "GPa",
    "kg",
    "sec",
    "deg",
    "psi",
    "lb",
    "lbf",
    "kip",
    "ksi",
    "inch",
    "ft",
]
