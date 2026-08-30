# ADR-0005: Equation signatures — floats inside, pint at the boundary

Status: accepted (2026-08-29)

## Context

The Phase 1 example in the roadmap originally showed equation functions typed
`-> Quantity`, i.e. pint end to end. That signature was assumed, not decided —
and it is the decision that most affects mako, because it fixes the cost of every
call in the station loop. Profiling (2026-08-29, rame-env / Python 3.12, commit
6299b78) settled it:

- `design_flexure` costs **375 ms/call**; `check_shear` **8 ms**; `check_flexure`
  **7 ms**. Roughly 65 % of profiler self-time is inside pint (~967k pint-internal
  calls vs ~2.8k mento calls per design).
- Of the design cost, **97 % is the rebar selection search in `rebar.py`**, whose
  inner loops do pint arithmetic. The flexure equations themselves cost ~6 ms.
- A representative flexure chain (ρ_min → A_s → M_n → DCR) evaluated 2000 times:
  pint with string units **1071 µs/call**, pint with prebuilt unit objects
  **313 µs** (3.4×), plain floats **0.9 µs** (~1200×).

The intermediate fix — replacing the ~384 `.to("cm")` string conversions with
`.to(cm)` unit objects — was considered as a way to keep pint in the engine. At
313 µs per chain it is still ~350× off the float floor, which matters at mako
scale (20,000 checks per structure). concrete-properties, the closest production
reference, reaches the same conclusion structurally: its engine is float-only and
units exist purely for display (`UnitDisplay`). Adopting its display pattern while
keeping pint in the engine would be an undeclared divergence on exactly the axis
that costs the most.

## Decision

**`codes/<code>/equations/` functions take and return plain floats in a
consistent unit system: N, mm, MPa (= N/mm²), N·mm.** Each docstring states the
expected unit per argument alongside the clause citation.

*Amended 2026-08-29, on extracting the first module.* The system cannot be a
single fixed one. ACI 318-19 publishes separately rounded coefficients for SI
and US customary — `0.17*sqrt(f_c[MPa])` against `2*sqrt(f_c[psi])`, which
differ by 2.3 % rather than being conversions of each other — and mento
reproduces both because its validation suite checks against worked examples in
both. So an equation whose coefficients differ takes an `imperial: bool`
keyword and works in either (N, mm, MPa) or (lb, in, psi); the *unit system* is
a parameter, the *units within it* are fixed and documented per argument. This
is a property of the design codes, not of mento: an engineer comparing the two
forms against the printed code should find both, unconverted.

**pint stays everywhere users touch mento.** Materials, forces, sections,
settings, and the result dataclasses of ADR-0001 keep Quantities. The checker
(ADR-0002) is the boundary: it converts inputs once on entry
(`f_c.to(MPa).magnitude`, using prebuilt unit objects, never strings), calls the
float equations, and wraps outputs back into Quantities when assembling the
result. Imperial inputs need no special handling — pint normalizes psi/kip to the
base system in that same single conversion.

Code that does not migrate into `equations/` (presentation, rebar orchestration)
still applies the mechanical `.to("x")` → `.to(x)` fix; it is free and worth 3–6×
on those paths.

## Consequences

- Equation evaluation drops ~1000×; the mako loop's per-station cost becomes
  boundary conversion + dataclass assembly, not arithmetic. Target: 20,000 checks
  (1000 stations × 2 faces × 10 combos) in **< 5 s serial** after Phase 2b.
- Equations become dependency-free `float → float` functions — the same shape
  structuralcodes uses, and the easiest possible unit for an engineer to verify
  against the printed code.
- Risk: unit mistakes inside equations are no longer caught by pint. Mitigation
  is the existing contract: every equation ships parametrized tests against the
  code's own tables (which are stated in known units), plus the per-argument unit
  in the docstring. A wrong unit cannot pass a table test.
- The checker gains a second responsibility besides factors and orchestration:
  unit boundary. The ADR-0002 rule extends to: **the checker does not compute,
  and the equations do not convert.**
