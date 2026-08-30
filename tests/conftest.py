"""Shared fixtures and test-suite configuration.

A fixture belongs here when more than one test module needs it and defines it
identically. Fixtures tied to a specific validated example — a beam from a
Calcpad case, a wall from a code table — stay in the module that asserts against
them, next to the expected values they were chosen to reproduce.
"""

import matplotlib
import pytest

from mento.material import Concrete, SteelBar
from mento.units import MPa

# Force a non-interactive backend for the whole suite, before any test module
# imports pyplot. Several modules exercise the plotting code, and without this
# they depend on whatever backend the machine happens to default to.
matplotlib.use("Agg")


@pytest.fixture()
def concrete_c25() -> Concrete:
    """Generic C25 concrete, with no design code attached."""
    return Concrete("C25")


@pytest.fixture()
def steel_b500s() -> SteelBar:
    """B500S reinforcement — the grade used by the EN examples."""
    return SteelBar(name="B500S", f_y=500 * MPa)


@pytest.fixture()
def steel() -> SteelBar:
    """ADN 420 reinforcement — the grade used by the ACI/CIRSOC examples."""
    return SteelBar(name="ADN 420", f_y=420 * MPa)
