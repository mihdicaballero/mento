# ADR-0004: Presentation is a separate layer, not a deleted feature

Status: accepted (2026-08-26)

## Context

Word reports, section/rebar plots, styled tables, and pint-formatted output are
core differentiators of mento for production office use — this is where mento
deliberately goes beyond the academic reference projects (structuralcodes ships
no plotting at all). But today that presentation code is smeared through the
calculation layers: ~30 % of `beam.py` is matplotlib and ~360 lines are docx
assembly; 29–38 % of the per-code modules are `_initialize_dicts_*` report-table
builders; a code module even calls a string formatter defined on the beam.

The consequences are concrete: importing an element drags in matplotlib, code
modules cannot be reviewed as "just the code equations", and the mako loop would
load display machinery it never uses.

## Decision

Keep every presentation feature; move all of it into presentation modules
(`plots/`, `reports/`, plus the already-clean `results.py`) that consume the
result dataclasses of ADR-0001 and nothing else.

- `_initialize_dicts_*` table builders move out of `codes/` into the report layer.
- Plotting moves out of `beam.py` into a plotting module, with a
  `plotting_context`-style helper (concrete-properties pattern) centralizing
  matplotlib boilerplate.
- Calculations run in pint quantities; **display rounding and display-unit choice
  belong to the presentation layer** (concrete-properties' `UnitDisplay` idea).

Dependency rule: presentation imports results; calculation layers never import
presentation. matplotlib/docx are imported only inside presentation modules.

## Consequences

- `codes/` shrinks to equations + checkers and becomes auditable clause by clause.
- The mako loop runs without display dependencies on the import path.
- Result dataclasses must carry everything reports need — a forcing function that
  keeps them complete (demands, capacities, ratios, intermediate values).
- Cost: methods users may call today on elements (`beam.plot()`, `*_doc()`)
  become thin delegating wrappers or move, following the deprecation policy of
  ADR-0003.
