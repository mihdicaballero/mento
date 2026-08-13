"""Deprecated alias for :mod:`mento.beam_summary`.

The module was renamed to ``mento.beam_summary`` so that it matches
``mento.shear_wall_summary``. Importing from here still works, and
``from mento import BeamSummary`` is unaffected, but this shim will be
removed in a future release.
"""

import warnings

from mento.beam_summary import BeamSummary

warnings.warn(
    "mento.summary has been renamed to mento.beam_summary. Import BeamSummary from "
    "mento.beam_summary, or simply use `from mento import BeamSummary`.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["BeamSummary"]
