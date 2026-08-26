# ADR-0002: Hybrid code layer — pure equations + per-code checker class

Status: accepted (2026-08-26)

## Context

The per-code modules must scale to more codes (CIRSOC as a real implementation,
future codes) and more verifications (torsion, punching), contributed by
structural engineers with limited Python experience. Today's pattern — functions
taking `self: RectangularBeam` and mutating it — has the highest possible
contribution threshold: touching one verification requires understanding the
whole beam's private state.

Alternatives considered, from the two reference projects:

- **A. Pure functions only** (structuralcodes): per-code subpackages of
  `float → float` functions citing clauses. Maximally testable and auditable, but
  orchestration, safety factors, and the rich result assembly mento's Word
  reports need have no natural home — structuralcodes can afford this because it
  produces no reports.
- **B. One class per code** (concrete-properties `DesignCode`): factors and
  orchestration live together, but equations become private methods — less
  individually testable, and the class trends toward the same 1500-line god the
  current modules are.

## Decision

Both, with a strict boundary. Per code, a subpackage:

- `codes/<code>/equations/*.py` — pure functions, one clause per docstring,
  parametrized tests against the code's own tables. No element objects, no state.
- `codes/<code>/checker.py` — a light class holding safety factors (φ, γ),
  calling equation functions in order, and assembling the result dataclasses of
  ADR-0001.

Boundary rule: **the checker does not compute.** If algebra appears in a checker
method, it moves to `equations/`.

## Consequences

- Contributors touch one small file per equation; an engineer can verify a
  formula against the printed code without reading the system.
- Reports and mako get one obvious entry point per code
  (`ACI318_19().check_shear(...)`).
- The `if design_code == ...` chains in `beam.py` collapse into a checker
  dispatch, and later into the Phase-4 registry.
- Cost: two concepts instead of one, and the boundary needs review discipline.
