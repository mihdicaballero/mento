import forallpeople as si
si.environment('structural')

# Metric system units
m = si.m
cm = 1e-2 * si.m
mm = 1e-3 * si.m
kN = 1e3 * si.N
MPa = 1e6 * si.Pa
kg = si.kg
sec = si.s

# Imperial system units
psi = si.psi
lb = si.lb
kip = 1e3 * si.lb
psi = si.psi
ksi = 1e3 * si.psi
inch = si.inch
ft = si.ft

def main() -> None:
    print(2*cm, 3*MPa, 4*kg, 1*cm, 1*m, 3*kN)
    print(1*psi, 1*lb, 1*kip, 1*psi, 1*ksi, 1*inch, 1*ft)

if __name__ == "__main__":
    main()