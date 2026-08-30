"""Pure EN 1992-2004 equations.

Every function here is a plain function of plain floats: no element objects, no
pint, no state. Each one cites the clause it comes from and is tested against
values taken from the code itself.

**Units.** EN 1992-1-1 is written in a single, coherent SI system and prints no
US customary variant of any coefficient, so — unlike the ACI package — nothing
here takes an ``is_imperial`` keyword. Arguments and results are floats in:

===========  ==============
quantity     unit
===========  ==============
force        N
length       mm
area         mm²
stress       MPa (= N/mm²)
moment       N·mm
===========  ==============

Reinforcement "areas per unit length" (A_sw/s) are therefore mm²/mm, which is
dimensionally a length — the same convention the ACI package uses.

Where EN publishes a *recommended* value that a National Annex may change
(C_Rd,c = 0.18/gamma_c, k_1 = 0.15, v_min from Eq. (6.3N), rho_max = 0.04),
mento uses the recommended value; the clause citation names it so a reader can
see what would have to change for another Annex.
"""
