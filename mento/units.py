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
