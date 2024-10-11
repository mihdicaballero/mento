from devtools import debug
from pint import UnitRegistry

# Initialize UnitRegistry
ureg = UnitRegistry(system='mks')
ureg.formatter.default_format = '.2f~P'  # Standardize output formatting

# Metric system units
m = ureg.meter
cm = ureg.centimeter
mm = ureg.millimeter
kN = 1e3 * ureg.newton
kNm = kN * m
MPa = ureg.megapascal
GPa = ureg.gigapascal
kg = ureg.kilogram
sec = ureg.second

# Imperial system units
psi = ureg.psi
lb = ureg.pound_force
kip = 1e3 * lb
ksi = ureg.ksi
inch = ureg.inch
ft = ureg.foot

def main() -> None:
    debug(2*cm, 3*MPa, 4*kg, 1*mm, 1*m, 3*kN, 2*kNm)
    debug(1*psi, 1*lb, 1*kip, 1*psi, 1*ksi, 1*inch, 1*ft)
    a = 3*ureg.meter
    N = 3*kN
    A = 3*m**2
    s = N/A
    print(a, N.to_compact(), N.to('kip'),A, s)
    # Get only the value
    debug(N.magnitude, N.to('kN').magnitude)
    debug(a, N, A, s)
    wavelength = 1550 * ureg.nm
    frequency = (ureg.speed_of_light / wavelength).to('Hz')
    print(frequency)
    print(frequency.to_compact())

    # How to format specific output
    debug('Length is {:~P}'.format(a))
    # Print all available units
    # print(dir(ureg.sys.mks))


if __name__ == "__main__":
    main()