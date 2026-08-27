# ADR-0001: Checks return immutable result dataclasses

Status: accepted (2026-08-26)

## Context

Today, check/design functions in `codes/` mutate ~164 private attributes on the
element instance (`self._phi_V_c = ...`) and return `None` or a display-oriented
DataFrame. Each combination overwrites the previous one's state, which forced two
workarounds already in the codebase: the envelope accumulators in `beam.py`
(whose own docstring documents the aliasing problem) and `design_results.py`, a
frozen-dataclass facade that *reads* the private state after the fact and needs a
`DesignNotRunError` to detect unmutated beams.

mako will run mento over ~1000 stations of a continuous beam per load combination.
A loop like that cannot be built on checks that share mutable state on one object.

## Decision

Every check and design function builds and **returns** a frozen dataclass
(`ShearCheckResult`, `FlexureCheckResult`, `ShearDesign`, `FlexureDesign`, ...),
defined in one module (`design_results.py`, extended). Result objects carry
demands, capacities, ratios, governing clauses, and the intermediate values that
reports need. Elements are not mutated by checks.

During migration (Phase 2 of the roadmap), the documented private attributes are
kept synchronized as a deprecated compatibility layer, then removed in a minor
release per ADR-0003.

## Consequences

- Envelopes become plain code over a list of results (`max(r.DCR for r in rs)`);
  the stateful accumulators and `DesignNotRunError` are deleted by construction.
- Checks become stateless and parallelizable — the mako loop is a list
  comprehension.
- Presentation (Word, plots, tables) reads result objects instead of private
  state, enabling ADR-0004.
- Both reference projects (concrete-properties `results.py`, structuralcodes
  `core/_section_results.py`) use exactly this pattern.
- Cost: the compatibility sync layer must be maintained until removal, and tests
  that assert against DataFrame cells migrate to asserting against result fields.
