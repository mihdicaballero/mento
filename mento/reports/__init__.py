"""Report assembly — the presentation layer of mento.

Three output media, one per module, all of them consumers of results rather
than producers of them:

- :mod:`~mento.reports.tables` — the labelled, rounded, unit-annotated tables a
  check produces.
- :mod:`~mento.reports.views` — the Markdown and console output for a notebook.
- :mod:`~mento.reports.documents` — the Word reports for a single element.
- :mod:`~mento.reports.summaries` — the Word reports for a summary of many.

Nothing here is imported by ``codes/`` or by the element classes' calculation
paths; the dependency only points this way.
"""
