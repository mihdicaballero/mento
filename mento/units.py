from pint import UnitRegistry

# Initialize UnitRegistry
ureg = UnitRegistry(system='mks')
ureg.formatter.default_format = '.2f~P'  # Standardize output formatting

# Metric system units
m = ureg.meter
cm = ureg.centimeter
mm = ureg.millimeter
kN = ureg.kilonewton
kNm = kN * m
Pa = ureg.pascal
MPa = ureg.megapascal
GPa = ureg.gigapascal
kg = ureg.kilogram
sec = ureg.second

# Imperial system units
psi = ureg.psi
lb = ureg.pound
lbf = ureg.pound_force
kip = ureg.kip
ksi = ureg.ksi
inch = ureg.inch
ft = ureg.foot