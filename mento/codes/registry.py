"""What a design code has to provide, and which ones provide it.

Every element used to dispatch by comparing ``concrete.design_code`` against a
string literal -- in :mod:`mento.beam`, :mod:`mento.shear_wall`,
:mod:`mento.rebar` and the summaries -- so adding a code meant finding and
editing each of those chains, and forgetting one failed silently at runtime
rather than at import.

A code declares itself here instead. :class:`DesignCode` is the whole contract:
if a new code fills one in and registers it, every element picks it up, and
``tests/test_architecture_boundaries.py`` fails the build if an element ever
grows a string comparison again.

The hooks are plain functions rather than methods on a class, matching how
``codes/`` is already written (ADR-0002): the code modules hold free functions
typed as ``self: RectangularBeam``, and a registry of them adds no indirection
that has to be stepped through when reading a check.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Callable, Dict, Tuple, Type

if TYPE_CHECKING:
    from mento.material import Concrete


@dataclass(frozen=True)
class DesignCode:
    """One design code's metadata and the hooks the elements call.

    ``materials`` lists the concrete classes the code accepts; it is what makes
    a mismatch between a section's concrete and its code detectable rather than
    an obscure failure deep in a check.
    """

    title: str
    year: int
    materials: Tuple[Type["Concrete"], ...]

    #: Verification. Each returns the check state for one combination.
    check_shear: Callable[..., Any]
    check_flexure: Callable[..., Any]
    #: Copy a state back onto the element, for the report tables (ADR-0001).
    apply_shear_state: Callable[..., None]
    apply_flexure_state: Callable[..., None]
    #: Sizing. Each mutates the element, which is the point of a design.
    design_shear: Callable[..., None]
    design_flexure: Callable[..., None]
    #: Bar selection: the catalogue and the detailing rules differ per code.
    transverse_rebar: Callable[..., Any]
    longitudinal_rebar: Callable[..., Any]
    #: The zeroed result attributes this code's report tables expect to read.
    initialize_attributes: Callable[..., None]
    #: Shear on a wall. Not every code has one; see :func:`requires`.
    check_shear_wall: Callable[..., Any] | None = None
    apply_wall_shear_state: Callable[..., None] | None = None
    design_shear_wall: Callable[..., None] | None = None

    #: Presentation the code owns: the symbols it names its quantities with and
    #: the unit rows of its summary tables. Not a calculation, but it is
    #: per-code, so a new code must be able to supply it without editing a
    #: report module either.
    flexure_symbols: Dict[str, str] = field(default_factory=dict)
    units_row_shear: Dict[str, str] = field(default_factory=dict)
    units_row_flexure: Dict[str, str] = field(default_factory=dict)
    #: What this code calls the demands and capacities in a multi-beam summary.
    #: ACI says Mu / ØMn / ØVn where EN says MEd / MRd / VRd, and the summaries
    #: have to name their columns one way or the other.
    #: Keys: moment_demand, shear_demand, shear_demand_source, axial_demand,
    #: moment_capacity_top, moment_capacity_bot, shear_capacity.
    summary_columns: Dict[str, str] = field(default_factory=dict)
    #: The moment capacities off a checked beam, already keyed by column name.
    capacity_columns: Callable[..., Dict[str, float]] | None = None
    #: Shear symbols for the one-line markdown summary, and the row its capacity
    #: sits on in this code's detail table.
    shear_symbols: Dict[str, Any] = field(default_factory=dict)
    #: Summary columns this code does not report, dropped from the Word tables.
    summary_drop_columns: Tuple[str, ...] = ()
    #: Smallest stirrup bar this code's detailing rules allow. Only the codes
    #: whose report tables quote it need to supply one.
    min_stirrup_diameter: Callable[..., Any] | None = None

    def requires(self, hook: str) -> Callable[..., Any]:
        """The hook, or a clear error naming the code that lacks it."""
        value = getattr(self, hook)
        if value is None:
            raise NotImplementedError(f"{hook.replace('_', ' ')} is not implemented for design code: {self.title}")
        return value


_REGISTRY: Dict[str, DesignCode] = {}
_DISCOVERED = False


def _discover() -> None:
    """Import every ``mento/codes/<code>/code.py``, which registers itself.

    Discovery rather than a hand-maintained list, so adding a code is adding a
    subpackage and nothing else -- the exit criterion of the roadmap's Phase 4.

    Lazy, and deliberately so: the code modules import ``mento.rebar``, which
    reaches back into ``codes/`` for its equations. Importing them from this
    module at import time would rebuild the cycle that emptying
    ``codes/__init__.py`` removed. By the time anything asks for a code, both
    ends of that chain are already imported.
    """
    global _DISCOVERED
    if _DISCOVERED:
        return
    _DISCOVERED = True  # set first: a code module importing this one must not recurse
    import importlib
    import pkgutil

    import mento.codes

    for module in pkgutil.iter_modules(mento.codes.__path__):
        if module.ispkg:
            importlib.import_module(f"mento.codes.{module.name}.code")


def register(code: DesignCode) -> DesignCode:
    """Add a code to the registry, keyed by the title elements dispatch on."""
    if code.title in _REGISTRY:
        raise ValueError(f"design code already registered: {code.title}")
    _REGISTRY[code.title] = code
    return code


def design_code(concrete: "Concrete") -> DesignCode:
    """The registered code for this concrete."""
    _discover()
    try:
        return _REGISTRY[concrete.design_code]
    except KeyError:
        raise NotImplementedError(
            f"unknown design code: {concrete.design_code}. Registered: {', '.join(sorted(_REGISTRY)) or 'none'}"
        ) from None


def registered_codes() -> Tuple[str, ...]:
    """Every registered code's title, for diagnostics and tests."""
    _discover()
    return tuple(sorted(_REGISTRY))
