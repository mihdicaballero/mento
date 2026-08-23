Beam — ACI 318-19
=================

Provisions implemented for ``RectangularBeam`` when the section is built with
``Concrete_ACI_318_19`` or ``Concrete_CIRSOC_201_25``.

.. admonition:: CIRSOC 201-25
   :class: note

   ``Concrete_CIRSOC_201_25`` subclasses ``Concrete_ACI_318_19``. Every equation on
   this page applies unchanged. The differences are that CIRSOC is metric-only and
   ships a different reinforcing bar catalogue (it allows Ø6 mm stirrups, which the
   ACI catalogue does not). Results tables report ``design_code = "CIRSOC 201-25"``.

Scope
-----

**Implemented**

- Flexure: nominal moment of singly and doubly reinforced rectangular sections,
  required steel area, minimum and maximum reinforcement limits, ductility check.
- One-way shear: concrete contribution :math:`V_c` per Table 22.5.5.1, steel
  contribution, upper bound on total shear, minimum stirrups and spacing limits.

**Not implemented**: torsion (Ch. 22.7), deflection and crack control (Ch. 24),
development and splice lengths (Ch. 25), seismic detailing (Ch. 18), fatigue,
two-way action in the beam module (see ``PunchingSlab``), and any strut-and-tie
treatment of deep beams — the shear provisions here assume sectional behaviour.

Symbols
-------

.. list-table::
   :header-rows: 1
   :widths: 14 46 40

   * - Symbol
     - Meaning
     - Source
   * - :math:`f'_c`
     - Specified compressive strength of concrete
     - ``concrete.f_c``
   * - :math:`f_y`
     - Specified yield strength of longitudinal reinforcement
     - ``steel_bar.f_y``
   * - :math:`f_{yt}`
     - Yield strength of transverse reinforcement, capped (§20.2.2.4)
     - ``_calculate_f_yt_aci``
   * - :math:`\beta_1`
     - Stress-block depth factor (Table 22.2.2.4.3)
     - ``concrete.beta_1``
   * - :math:`\varepsilon_{cu}`
     - Ultimate concrete strain, 0.003 (§22.2.2.1)
     - ``concrete._epsilon_c``
   * - :math:`\varepsilon_y`
     - Reinforcement yield strain, :math:`f_y/E_s`
     - ``steel_bar.epsilon_y``
   * - :math:`\phi_t,\ \phi_v`
     - Strength reduction factors, 0.90 and 0.75 (§21.2)
     - ``concrete.phi_t``, ``concrete.phi_v``
   * - :math:`\lambda`
     - Lightweight concrete factor
     - ``concrete.lambda_factor``
   * - :math:`d,\ d'`
     - Effective depth to tension steel / cover to compression steel
     - ``_d_bot``, ``_c_mec_top`` (or swapped)
   * - :math:`a,\ c`
     - Stress block depth and neutral axis depth, :math:`a = \beta_1 c`
     - —
   * - :math:`\rho_w`
     - Longitudinal tension reinforcement ratio :math:`A_s/(b_w d)`
     - ``_rho_w``

Flexure
-------

Stress block
^^^^^^^^^^^^

Section strength uses the equivalent rectangular stress distribution of §22.2.2.4:
a uniform stress of :math:`0.85 f'_c` over a depth :math:`a = \beta_1 c`, with
:math:`\beta_1` from Table 22.2.2.4.3:

.. math::

   \beta_1 = \begin{cases}
   0.85 & f'_c \le 28\ \text{MPa} \\[2pt]
   0.85 - \dfrac{0.05}{7}\,(f'_c - 28) & 28 < f'_c \le 55\ \text{MPa} \\[2pt]
   0.65 & f'_c > 55\ \text{MPa}
   \end{cases}

Nominal moment — singly reinforced
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From force equilibrium :math:`0.85 f'_c\, a\, b = A_s f_y`:

.. math::

   a = \frac{A_s f_y}{0.85 f'_c b}
   \qquad
   M_n = A_s f_y \left(d - \frac{a}{2}\right)

Implemented in ``_determine_nominal_moment_simple_reinf_ACI_318_19``. Used when
:math:`A_s \le A_{s,max}`; the tension steel is assumed to yield, which
:math:`A_{s,max}` guarantees.

Nominal moment — doubly reinforced
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When :math:`A_s > A_{s,max}` the compression steel enters the resisting couple.
Equilibrium is written with a **displaced-concrete correction**: the compression bar
occupies a volume already accounted for by the :math:`0.85 f'_c a b` block, so its
effective contribution is :math:`(f_y - 0.85 f'_c)`, not :math:`f_y`.

Assuming the compression steel yields:

.. math::

   c = \frac{A_s f_y - A'_s\,(f_y - 0.85 f'_c)}{0.85 f'_c\, b\, \beta_1}
   \qquad
   \varepsilon'_s = \frac{c - d'}{c}\,\varepsilon_{cu}

If :math:`\varepsilon'_s \ge \varepsilon_y` the assumption holds and

.. math::

   M_n = 0.85 f'_c\, a\, b \left(d - \frac{a}{2}\right)
       + A'_s\,(f_y - 0.85 f'_c)(d - d')

Otherwise the compression steel does not yield and :math:`c` is recovered from the
quadratic that results from substituting :math:`f'_s = \varepsilon_{cu} E_s (c-d')/c`
into equilibrium:

.. math::

   \underbrace{0.85 f'_c b \beta_1}_{A}\,c^2
   + \underbrace{\left[A'_s\,(\varepsilon_{cu} E_s - 0.85 f'_c) - A_s f_y\right]}_{B}\,c
   + \underbrace{\left(-d' A'_s \varepsilon_{cu} E_s\right)}_{C} = 0

with :math:`f'_{s,net} = f'_s - 0.85 f'_c` used for the couple. Implemented in
``_determine_nominal_moment_double_reinf_ACI_318_19``.

Required reinforcement
^^^^^^^^^^^^^^^^^^^^^^

Inverting the singly-reinforced expression for a demand :math:`M_u`:

.. math::

   R_n = \frac{M_u}{\phi_t\, b\, d^2}
   \qquad
   A_{s,calc} = \frac{0.85 f'_c b d}{f_y}
                \left(1 - \sqrt{1 - \frac{2 R_n}{0.85 f'_c}}\right)

A negative radicand means the section cannot develop :math:`M_u` as singly
reinforced; mento then sets :math:`A_{s,calc} = A_{s,max}` so the calculation
completes and the downstream check reports ``DCR > 1`` instead of raising.

Minimum reinforcement
^^^^^^^^^^^^^^^^^^^^^

§9.6.1.2 gives the mechanical minimum:

.. math::

   \rho_{min} = \max\left(\frac{0.25\sqrt{f'_c}}{f_y},\ \frac{1.4}{f_y}\right)
   \quad\text{(MPa)}
   \qquad
   \max\left(\frac{3\sqrt{f'_c}}{f_y},\ \frac{200}{f_y}\right)
   \quad\text{(psi)}

§9.6.1.3 allows :math:`A_{s,min}` to be waived if the steel provided is at least
one third greater than required, i.e. :math:`\tfrac{4}{3}A_{s,calc}`. mento applies
this only when it actually saves steel and the result still clears the geometric
minimum below — see :ref:`aci-decisions`.

Maximum reinforcement and ductility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The upper limit keeps the section tension-controlled (CRSI Guide, Beam Theory p. 6-3,
consistent with §21.2.2):

.. math::

   \rho_{max} = 0.85\,\beta_1\,\frac{f'_c}{f_y}\,
                \frac{\varepsilon_{cu}}{\varepsilon_y + 2\varepsilon_{cu}}

The neutral axis at the tension-controlled boundary follows from strain
compatibility with :math:`\varepsilon_t = \varepsilon_y + \varepsilon_{cu}`:

.. math::

   c_t = \frac{0.003\, d}{\varepsilon_y + 0.006}

Sections with :math:`c < c_t` are ductile. The ratio :math:`c/d` is reported in the
results tables.

Compression steel at the ductility limit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When compression reinforcement is required, its stress is evaluated at :math:`c_t`
and corrected for displaced concrete:

.. math::

   f'_s = \min\!\left(0.003\,E_s\left(1 - \frac{d'}{c_t}\right),\ f_y\right)
   \qquad
   f'_{s,net} = f'_s - 0.85 f'_c

.. math::

   A'_s = \frac{M_u/\phi_t - M_{n,t}}{f'_{s,net}\,(d - d')}
   \qquad
   M_{n,t} = \rho f_y \left(d - \frac{a_{max}}{2}\right) b\, d,
   \quad a_{max} = \beta_1 c_t

Shear
-----

Concrete contribution
^^^^^^^^^^^^^^^^^^^^^

Table 22.5.5.1, in the two-branch form that depends on whether minimum stirrups are
present. With the axial term :math:`\sigma_{Nu} = \min\!\big(N_u/(6 A_g),\,0.05 f'_c\big)`:

When :math:`A_v \ge A_{v,min}`:

.. math::

   k_c = \max\left(0.17\lambda\sqrt{f'_c},\
                   0.66\,\lambda\,\rho_w^{1/3}\sqrt{f'_c}\right) + \sigma_{Nu}

When :math:`A_v < A_{v,min}`, the size effect factor :math:`\lambda_s` applies
(§22.5.5.1.3):

.. math::

   k_c = 0.66\,\lambda_s\,\lambda\,\rho_w^{1/3}\sqrt{f'_c} + \sigma_{Nu}
   \qquad
   \lambda_s = \min\!\left(\sqrt{\frac{2}{1 + 0.004\,d}},\ 1.0\right)

(:math:`d` in mm; the imperial form is :math:`\sqrt{2/(1 + d/(10\,\text{in}))}`.)

.. note::

   The :math:`\le 1.0` cap is part of the code equation and is easy to drop. The
   size effect factor exists to **reduce** :math:`V_c` in deep members; left
   uncapped it exceeds 1 for :math:`d < 250` mm and inflates :math:`V_c` instead.
   The combination that reaches this branch — no minimum stirrups **and** a shallow
   section — is in practice a slab, which is exactly where it matters most.

The result is bounded above (§22.5.5.1.1):

.. math::

   V_c = \min\Big(0.42\,\lambda\sqrt{f'_c}\,A_{cv},\ \max(0,\ k_c A_{cv})\Big)
   \qquad A_{cv} = b_w\, d

.. important::

   :math:`\rho_w` uses the reinforcement on the **tension** face for the combination
   being checked — bottom steel for a positive moment, top steel for a negative one.
   A combination with no moment therefore uses the bottom steel. Because
   :math:`V_c \propto \rho_w^{1/3}`, a section with no longitudinal steel would
   compute :math:`V_c = 0`; see :ref:`aci-decisions`.

Steel contribution and total strength
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   V_s = A_v f_{yt} d
   \qquad
   \phi V_n = \phi_v \left(V_c + A_v f_{yt} d\right)

with :math:`f_{yt} \le 420` MPa (60 ksi) per §20.2.2.4. The section is also capped
by §22.5.1.2:

.. math::

   \phi V_{max} = \phi_v\left(V_c + 0.66\,\lambda\sqrt{f'_c}\,A_{cv}\right)

and the reported DCR uses :math:`\min(\phi V_n,\ \phi V_{max})` as the capacity, so
exceeding the upper bound is surfaced rather than silently allowed.

Minimum shear reinforcement
^^^^^^^^^^^^^^^^^^^^^^^^^^^

§9.6.3.1 — no minimum stirrups are required while

.. math::

   V_u < 0.083\,\phi_v\,\lambda\sqrt{f'_c}\,A_{cv}
   \quad\text{(MPa)}

Above that threshold, Table 9.6.3.4 gives

.. math::

   A_{v,min} = \max\left(0.062\sqrt{f'_c}\,\frac{b_w}{f_{yt}},\
                         0.35\,\frac{b_w}{f_{yt}}\right)
   \quad\text{(MPa)}

and the required area is

.. math::

   A_{v,req} = \max\left(\frac{V_u - \phi V_c}{\phi_v\, f_{yt}\, d},\ A_{v,min}\right)

Spacing limits
^^^^^^^^^^^^^^

Table 9.7.6.2.2, keyed on :math:`V_s` demand:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Condition
     - :math:`s_{max}` along length
     - :math:`s_{max}` across width
   * - :math:`V_s \le 0.083\lambda\sqrt{f'_c}A_{cv}`
     - :math:`\min(d/2,\ 600\ \text{mm})`
     - :math:`\min(d,\ 600\ \text{mm})`
   * - :math:`V_s > 0.083\lambda\sqrt{f'_c}A_{cv}`
     - :math:`\min(d/4,\ 300\ \text{mm})`
     - :math:`\min(d/2,\ 300\ \text{mm})`

.. _aci-decisions:

Implementation decisions
------------------------

Where the code leaves room for interpretation, this is what mento does.

Geometric minimum of 1.8‰
^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the §9.6.1.2 mechanical minimum, mento enforces a **geometric**
minimum of :math:`1.8\text{‰}\,b\,h` on the tension face. This is not an ACI
provision; it is common detailing practice and it prevents two degenerate outcomes:

1. When §9.6.1.3's 4/3 rule would drop the steel below any sensible detailing
   amount, the geometric minimum takes over.
2. When :math:`M_u = 0` — a shear-only load combination — §9.6.1.2 requires no
   flexural steel at all, so the design would return :math:`A_s = 0`. That is not a
   buildable layout, and it makes :math:`\rho_w = 0`, which collapses :math:`V_c` to
   zero in Table 22.5.5.1. mento adopts the geometric minimum instead.

Order of the minimum-reinforcement rules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When :math:`A_{s,calc} < A_{s,min}`, mento evaluates :math:`\tfrac{4}{3}A_{s,calc}`
and adopts it **only if** it is both smaller than :math:`A_{s,min}` (so it actually
saves steel) and not smaller than the geometric minimum. Otherwise the governing
value is :math:`A_{s,min}`, or the geometric minimum, whichever applies. The results
tables flag whether the 4/3 rule was used.

Displaced-concrete correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compression reinforcement contributes :math:`(f_y - 0.85 f'_c)` rather than
:math:`f_y`, both in equilibrium and in :math:`M_n`. Many textbook treatments omit
this and slightly overestimate capacity. mento applies it consistently.

Exact :math:`a_{max}` instead of the 0.59 shortcut
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Whitney block depth at the ductility limit is taken as
:math:`a_{max} = \beta_1 c_t` exactly, rather than the common shortcut
:math:`d - 0.59\rho f_y d/f'_c`, which relies on :math:`0.59 \approx 1/1.7` and loses
roughly 0.3% of precision.

Effective depth from the real bar layout
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:math:`d` is recomputed from the centroid of the bars actually placed, and the design
iterates on it. Software that assumes a fixed
:math:`d` — ETABS, for instance — will report slightly different required areas for
the same demand. The difference is geometric, not a disagreement about the code.

For shear, :math:`d_{shear} = \min(d_{bot},\ d_{top})`, the conservative choice when
both faces carry steel.

Validation
----------

Every row is pinned by a test in ``tests/test_beam.py`` and was verified against the
external source named.

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Check
     - Test
     - Verified against
   * - :math:`M_n`, singly reinforced
     - ``test_determine_nominal_moment_simple_reinf_ACI_318_19_Test_Etabs_03``
     - ETABS + Excel *BEAM-01-Flexure-Rectangle ACI 318-19-v6*
   * - :math:`M_n`, doubly reinforced (steel not yielding)
     - ``test_determine_nominal_moment_double_reinf_ACI_318_19_Test_Etabs_01``
     - ETABS + Excel, same workbook
   * - :math:`M_n`, doubly reinforced (steel yielding)
     - ``test_determine_nominal_moment_double_reinf_ACI_318_19_Test_Etabs_23_yielding``
     - ETABS + Excel, same workbook
   * - :math:`A_{s,req}`, :math:`A_{s,min}` governing
     - ``test_calculate_flexural_reinforcement_ACI_318_19_Test_Etabs_05``
     - ETABS, :math:`A_{s,req} = 1.7041\ \text{in}^2`
   * - :math:`A_{s,req}`, doubly reinforced
     - ``test_calculate_flexural_reinforcement_ACI_318_19_doubly_reinforced_Test_Etabs_01``
     - ETABS
   * - :math:`\rho_{max}`
     - ``test_maximum_flexural_reinforcement_ratio_ACI_318_19_Test_Etabs_05``
     - CRSI Guide, Beam Theory p. 6-3
   * - :math:`\rho_{min}`
     - ``test_minimum_flexural_reinforcement_ratio_ACI_318_19_Test_Etabs_05``
     - §9.6.1.2
   * - Flexure check, full section
     - ``test_check_flexure_ACI_318_19_1`` … ``_3``
     - Excel workbook
   * - Flexure design, doubly reinforced
     - ``test_design_flexure_ACI_318_19_Test_Etabs_01``
     - ETABS, verified through :math:`\phi M_n \ge M_u`
   * - Shear check, with stirrups
     - ``test_shear_check_ACI_318_19_1``, ``_2``
     - Calcpad *ACI 318-19 Beam Shear 01 — Imperial*
   * - Shear check, no stirrups
     - ``test_shear_check_ACI_318_19_no_rebar_1``, ``_2``
     - Calcpad, same sheet
   * - Shear design
     - ``test_shear_design_ACI_318_19``
     - Calcpad, same sheet
   * - Shear design, CIRSOC catalogue
     - ``test_shear_design_CIRSOC_201_2025``
     - CIRSOC 201-25 bar catalogue
   * - Geometric minimum at :math:`M_u = 0`
     - ``test_design_flexure_ACI_318_19_zero_moment_adopts_geometric_minimum``
     - Internal consistency: :math:`\rho_w > 0`, :math:`V_c > 0`
   * - Never raises on an unfittable section
     - ``test_design_flexure_rebar_infeasible_does_not_crash``
     - Public contract
   * - :math:`\lambda_s \le 1.0` cap
     - ``test_lambda_s_is_capped_at_one_imperial`` / ``_metric``
       (in ``tests/test_slab.py``)
     - §22.5.5.1.3
