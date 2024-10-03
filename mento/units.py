from devtools import debug
from pint import UnitRegistry

# Initialize the environment
ureg = UnitRegistry(system='mks')
ureg.formatter.default_format = '.2f~P'

# Metric system units
m = ureg.meter
cm = 1e-2 * ureg.meter
mm = 1e-3 * ureg.meter
kN = 1e3 * ureg.newton
kNm = 1e3 * ureg.newton * ureg.meter
MPa = 1e6 * ureg.pascal
kg = 1e3 * ureg.gram
sec = ureg.second

# Imperial system units
psi = ureg.psi
lb = ureg.pound
kip = 1e3 * ureg.pound
ksi = 1e3 * ureg.psi
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
    debug(N.magnitude, N.to('kN'))
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