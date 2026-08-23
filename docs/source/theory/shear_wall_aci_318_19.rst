Shear Wall — ACI 318-19
=======================

Provisions implemented for ``ShearWall``, covering **in-plane shear** per
ACI 318-19 Chapter 11. Also valid for ``Concrete_CIRSOC_201_25``, which differs only
in the reinforcing bar catalogue used for design.

Scope
-----

**Implemented**: in-plane shear check and design — concrete contribution with the
aspect-ratio factor, distributed mesh contribution, upper bound on nominal shear,
minimum reinforcement ratios in both directions, and spacing limits.

.. warning::

   **Flexure is not implemented for walls.** ``ShearWall`` inherits flexural methods
   from ``RectangularBeam``, but those assume beam geometry and a lumped tension
   chord — they are not validated for a wall section with distributed vertical steel
   and must not be used. Also absent: boundary element design (§18.10.6),
   out-of-plane bending, sliding shear at construction joints (§11.5.4.4), coupling
   beams, and any seismic detailing of Chapter 18.

How it differs from a beam
--------------------------

Three things change relative to the beam shear model:

1. The concrete contribution uses an **aspect-ratio factor** :math:`\alpha_c` in
   place of the longitudinal reinforcement term. A squat wall mobilises far more
   concrete shear than a slender one, and :math:`\rho_w` is not the right variable
   to express that.
2. Reinforcement is a **distributed mesh in two orthogonal directions**, placed on
   both faces, not stirrups. The horizontal mesh :math:`\rho_t` resists shear; the
   vertical mesh :math:`\rho_l` has its own minimum tied to :math:`\rho_t`.
3. The effective area is the **gross in-plane area** :math:`A_{cv} = l_w t`, not
   :math:`b_w d`.

Symbols
-------

.. list-table::
   :header-rows: 1
   :widths: 16 44 40

   * - Symbol
     - Meaning
     - Source
   * - :math:`t`
     - Wall thickness (out-of-plane)
     - ``thickness``
   * - :math:`l_w`
     - Wall length (in-plane, resists shear)
     - ``length``
   * - :math:`h_w`
     - Wall height, for the aspect ratio
     - ``height``
   * - :math:`A_{cv}`
     - Gross area resisting shear, :math:`l_w t`
     - ``_Acv``
   * - :math:`\alpha_c`
     - Aspect-ratio factor (§11.5.4.6)
     - ``_alpha_c``
   * - :math:`\rho_t,\ \rho_l`
     - Horizontal and vertical distributed reinforcement ratios
     - ``_rho_t``, ``_rho_l``
   * - :math:`f_{yt}`
     - Yield strength of the mesh, capped
     - ``_f_yt_wall``

Shear strength
--------------

Aspect-ratio factor
^^^^^^^^^^^^^^^^^^^

§11.5.4.6 — linear interpolation on :math:`h_w/l_w`:

.. math::

   \alpha_c = \begin{cases}
   0.25 & h_w/l_w \le 1.5 \\[2pt]
   \text{linear} & 1.5 < h_w/l_w < 2.0 \\[2pt]
   0.17 & h_w/l_w \ge 2.0
   \end{cases}
   \quad \text{(MPa)}

with 3.0 and 2.0 as the corresponding imperial values (psi). Squat walls
(:math:`h_w/l_w \le 1.5`) get the higher factor.

Nominal strength
^^^^^^^^^^^^^^^^

.. math::

   V_c = \alpha_c\,\lambda\sqrt{f'_c}\,A_{cv}
   \qquad
   V_s = \rho_t\, f_{yt}\, A_{cv}
   \qquad
   V_n = V_c + V_s

bounded above by §11.5.4.3:

.. math::

   V_{n,max} = 0.66\,\lambda\sqrt{f'_c}\,A_{cv} \quad\text{(MPa)}
   \qquad
   8\,\lambda\sqrt{f'_c}\,A_{cv} \quad\text{(psi)}

The reported design strength is :math:`\phi_v \min(V_n,\ V_{n,max})` with
:math:`\phi_v = 0.75`, so exceeding the upper bound is surfaced in the DCR rather
than silently allowed.

Required horizontal reinforcement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Inverting :math:`\phi_v V_n \ge V_u` for the mesh ratio:

.. math::

   \rho_{t,req} = \max\left(
       \frac{V_u/\phi_v - \alpha_c \lambda \sqrt{f'_c}\, A_{cv}}
            {f_{yt}\, A_{cv}},
       \ \rho_{t,min}\right)

Minimum reinforcement
---------------------

Horizontal, §11.6.1 — a flat minimum:

.. math::

   \rho_{t,min} = 0.0025

Vertical, §11.6.2 — the minimum depends on how much horizontal steel the strength
check demanded, and on the aspect ratio:

.. math::

   \rho_{l,min} = \max\Big(0.0025,\
       0.0025 + 0.5\,(2.5 - h_w/l_w)\,(\rho_{t,req} - 0.0025)\Big)

with :math:`h_w/l_w` clamped to :math:`[0.5,\ 2.5]`. The two ends of that clamp are
the physically meaningful cases:

- :math:`h_w/l_w \ge 2.5` — slender wall, the interpolation term vanishes and only
  the flat minimum applies.
- :math:`h_w/l_w \le 0.5` — very squat wall, and :math:`\rho_{l,min}` rises to match
  :math:`\rho_{t,req}`. A squat wall carries shear through a diagonal strut that
  needs vertical steel to anchor it.

Per §11.5.4.3, :math:`\rho_{l,req}` need not exceed the :math:`\rho_t` required for
strength.

Spacing limits
--------------

§11.7.3, with the absolute cap being 450 mm (18 in):

.. math::

   s_{h,max} = \min\left(\frac{l_w}{5},\ 3t,\ 450\ \text{mm}\right)
   \qquad
   s_{v,max} = \min\left(\frac{l_w}{3},\ 3t,\ 450\ \text{mm}\right)

Implementation decisions
------------------------

Reinforcement on both faces
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mesh is assumed placed on **both faces** (E.F.), so the ratio is computed as
:math:`\rho = n_{curtains}\, A_b/(s\,t)` with ``_n_curtains`` fixed at 2. mento does
not currently model single-curtain walls; ACI §11.7.2.3 requires two curtains above a
shear demand threshold in any case.

CIRSOC bar catalogue
^^^^^^^^^^^^^^^^^^^^

The provisions are identical under CIRSOC 201-25. The design routine draws from a
different bar list: Ø6 mm is available for the horizontal mesh, and the vertical
mesh starts at Ø10 mm.

:math:`f_{yt}` cap
^^^^^^^^^^^^^^^^^^

As for beams, the yield strength used for the mesh is capped per §20.2.2.4, so a
high-strength mesh does not translate into unlimited shear strength.

Validation
----------

Tests live in ``tests/test_shear_wall.py``, organised in classes by topic.

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Check
     - Test
     - Verified against
   * - :math:`\alpha_c`, squat / slender / interpolated
     - ``test_alpha_c_low_hw_lw``, ``test_alpha_c_high_hw_lw``,
       ``test_alpha_c_interpolated``
     - §11.5.4.6
   * - :math:`V_c`, :math:`\phi V_c`, :math:`\phi V_n`, :math:`\phi V_{n,max}`
     - ``test_Vc``, ``test_phi_Vc``, ``test_phi_Vn``, ``test_phi_Vn_max``
     - §11.5.4.6, §11.5.4.3
   * - DCR and capacity flags
     - ``test_DCR``, ``test_Vu_le_phi_Vn``, ``test_Vu_le_phi_Vn_max``
     - Internal consistency
   * - :math:`\rho_{t,min}`
     - ``test_rho_t_min_is_0025``
     - §11.6.1
   * - :math:`\rho_{l,min}`, slender and interpolated
     - ``test_rho_l_min_is_0025_when_hw_lw_gt2``,
       ``test_rho_l_min_interpolated_per_11_6_2``
     - §11.6.2
   * - Spacing limits
     - ``test_s_h_max_metric``, ``test_s_v_max_metric``,
       ``test_spacing_ok_flag_when_within_limit``,
       ``test_spacing_fail_flag_when_exceeds_limit``
     - §11.7.3
   * - Mesh ratios from bar layout
     - ``test_set_horizontal_rebar_updates_rho_t``,
       ``test_set_vertical_rebar_updates_rho_l``
     - Internal consistency
