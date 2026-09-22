"""Drawing — the matplotlib half of the presentation layer.

One module per element family:

- :mod:`~mento.plots.sections` — beam and slab cross-sections.
- :mod:`~mento.plots.walls` — shear wall elevations.
- :mod:`~mento.plots.punching` — punching perimeters.

Kept apart from :mod:`mento.reports` so that a calculation-only import path
never pulls matplotlib in. That is enforced rather than hoped for: the elements
import these modules inside their ``plot()`` methods, and this ``__init__``
imports nothing, so ``import mento`` and a whole ``design()`` leave matplotlib
unloaded (``tests/test_lazy_imports.py``).
"""

from __future__ import annotations

#: The libraries a drawing needs and a calculation does not.
_PLOTTING_LIBRARIES = {"matplotlib", "seaborn"}

PLOTTING_IMPORT_ERROR = (
    "Plotting needs {name}, which is not installed. Install it with `pip install matplotlib seaborn`."
)


def plotting_import_error(error: ImportError) -> ImportError:
    """``error``, reworded if it is a plotting library that is missing.

    An environment without matplotlib is a legitimate one now that nothing but
    ``plot()`` needs it, and a bare ``No module named 'matplotlib'`` from three
    frames down does not say which call wanted it or what to install. Any other
    import failure is returned as it came.
    """
    name = (error.name or "").split(".")[0]
    if name in _PLOTTING_LIBRARIES:
        return ImportError(PLOTTING_IMPORT_ERROR.format(name=name), name=error.name)
    return error
