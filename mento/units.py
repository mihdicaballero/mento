from pint import UnitRegistry
import math

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

def main() -> None:
    print(45*deg)
    a= 45*ureg.degree
    print(a, math.cos(a))
    b = 0.000545*dimensionless
    print(b.magnitude)

if __name__ == "__main__":
    main()