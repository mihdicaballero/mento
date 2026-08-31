"""Pure ACI 318-19 equations.

Every function here is a plain function of plain floats: no element objects, no
pint, no state. Each one cites the clause or table it comes from and is tested
against values taken from the code itself.

**Units.** Arguments and results are floats in one of two consistent systems,
selected by the ``is_imperial`` keyword — pass ``concrete.is_imperial``:

===========  ==========================  =======================
quantity     SI (``is_imperial=False``)  US customary (``True``)
===========  ==========================  =======================
force        N                           lb
length       mm                          in
stress       MPa                         psi
===========  ==========================  =======================

ACI publishes separate, individually rounded coefficients for the two systems
— ``0.17*sqrt(f_c[MPa])`` against ``2*sqrt(f_c[psi])``, which differ by about
2 % rather than being exact conversions of each other. Both are reproduced as
printed, so a result computed in one system will not exactly equal the other
converted; that difference belongs to the code, not to mento.
"""
