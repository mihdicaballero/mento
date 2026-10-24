"""Column shear checked against cases published outside mento.

Marked ``published_example``, with a ``Source:`` paragraph in each docstring; see
``tests/architecture/test_published_examples.py``.
"""

import pytest

from mento.codes.en_1992_2004.equations import shear as eq


@pytest.mark.published_example
def test_the_share_of_a_circular_link_is_the_published_0_85() -> None:
    """λ1 of a closed circular link with z = 0.8·D, r_sv = 0.45·D and the chords
    equidistant from the centre (z_0 = z/2): Eq. (6) integrates to 0.8450,
    and the paper gives 0.85, 0.6 % above it -- rounded, as the paper writes
    every factor to two figures. The test holds the two to 1 %.

    Source: Orr, Darby, Ibell, Denton and Shave, "Shear design of circular concrete
    sections using the Eurocode 2 truss model", The Structural Engineer 88(23-24), 2010,
    §3.2, Eq. (7) and the paragraph before it: "an effectiveness factor of λ1 = 0.85
    was obtained", after Turmo, Ramos and Aparicio, Engineering Structures 31(2), 2009.
    """
    D = 1000.0
    assert eq.circular_link_efficiency(0.4 * D, 0.8 * D, 0.45 * D) == pytest.approx(0.85, rel=0.01)


@pytest.mark.published_example
def test_a_spiral_of_no_pitch_is_a_closed_link() -> None:
    """λ2 = ((p/(2π·r_sv))² + 1)^-0.5 is 1 when the pitch is zero.

    Source: Orr et al., The Structural Engineer 88(23-24), 2010, §3.2.1, Eq. (9) and
    the sentence after it: "When the pitch, p, is set to zero (as for closed links),
    λ2 is equal to unity".
    """
    assert eq.spiral_efficiency(0.0, 200.0) == 1.0
