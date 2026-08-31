"""Drawing — the matplotlib half of the presentation layer.

One module per element family:

- :mod:`~mento.plots.sections` — beam and slab cross-sections.
- :mod:`~mento.plots.walls` — shear wall elevations.
- :mod:`~mento.plots.punching` — punching perimeters.

Kept apart from :mod:`mento.reports` so that a calculation-only import path
never pulls matplotlib in.
"""
