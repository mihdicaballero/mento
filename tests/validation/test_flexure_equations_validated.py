"""Flexure clause functions checked against a case validated outside mento.

Marked ``published_example``, with a ``Source:`` paragraph in each docstring; see
``tests/architecture/test_published_examples.py``.
"""

import math
import pytest
from mento.codes.aci_318_19.equations import flexure as eq

#: The f_y cap of ACI 318-19 Table 20.2.2.4(a) on the minimum ratio, in psi.
ACI_FY_CAP_US = 80_000.0


@pytest.mark.published_example
def test_min_reinforcement_ratio_us_matches_the_validated_case():
    """Test_Etabs_05 of the beam suite: f_c = 6000 psi, f_y = 60 ksi.

    3*sqrt(6000)/60000 = 0.003873, which governs over 200/60000 = 0.003333.

    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion,
    row 31 (Test_Etabs_05), column BC: A_s = A_s,min = 1.7041 in² = 0.003873 * 16 in * 27.5 in.
    """
    got = eq.min_reinforcement_ratio(6000.0, 60_000.0, ACI_FY_CAP_US, is_imperial=True)
    assert got == pytest.approx(0.003873, rel=1e-4)
    assert got == pytest.approx(3 * math.sqrt(6000.0) / 60_000.0, rel=1e-12)
