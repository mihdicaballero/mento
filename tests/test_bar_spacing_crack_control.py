"""The crack-control cap on bar spacing: ACI 318-19 / CIRSOC 201-25 §24.3.2.

§7.7.2.2 (one-way slabs) and §9.7.2.2 (beams) send the bars closest to the
tension face to Table 24.3.2: s <= min(380*(280/f_s) - 2.5*c_c, 300*(280/f_s)),
with f_s = (2/3)*f_y permitted by §24.3.2.1. mento used to apply the 3h and
450 mm of §7.7.2.3 alone to slabs, and nothing to beams. Every number below is
worked by hand in its docstring.
"""

import pytest

from mento import (
    Concrete_ACI_318_19,
    Concrete_CIRSOC_201_25,
    Concrete_EN_1992_2004,
    Forces,
    Node,
    OneWaySlab,
    SteelBar,
)
from mento.slab import Footing
from mento.units import MPa, cm, kNm, mm

ADN_420 = SteelBar(name="ADN 420", f_y=420 * MPa)


def _slab(concrete, height, steel=ADN_420, c_c=25 * mm):  # type: ignore[no-untyped-def]
    return OneWaySlab(label="L1", concrete=concrete, steel_bar=steel, width=100 * cm, height=height, c_c=c_c)


# ---------------------------------------------------------------------------
# Slabs: the cap enters the design limit beside §7.7.2.3
# ---------------------------------------------------------------------------


def test_an_aci_slab_is_held_to_300_mm_by_table_24_3_2() -> None:
    """A 12 cm ACI slab under 7.5 kN·m, f'c 25, f_y 420, c_c 25 mm.

    §7.7.2.3 allows min(3h, 450) = 360 mm. §7.7.2.2 sends the bars to Table
    24.3.2 with f_s = (2/3)*420 = 280 MPa and c_c = 25 mm (no stirrup):
    min(380*1 - 62.5, 300*1) = min(317.5, 300) = 300 mm. The three Ø10 the
    search picks for the 2.26 cm² required would sit at floor(100/3) = 33 cm,
    inside the 36 cm of §7.7.2.3 and past the 30 cm of §24.3.2, so the strip
    is detailed Ø10/30: 3.33 bars, 2.62 cm².
    """
    slab = _slab(Concrete_ACI_318_19(name="H25", f_c=25 * MPa), 12 * cm)
    Node(section=slab, forces=Forces(label="C1", M_y=7.5 * kNm)).design()

    assert slab._max_bar_spacing().to("mm").magnitude == pytest.approx(300.0)
    layer = slab.reinforcement.bottom.layers[0]
    assert (layer.d_b, layer.s.to("cm").magnitude) == (10 * mm, 30)
    assert slab.reinforcement.bottom.A_s.to("cm**2").magnitude == pytest.approx(2.62, abs=5e-3)
    assert slab.reinforcement.bottom.A_s >= slab.flexure_design.bottom.A_s_req
    assert slab.warnings == ()


@pytest.mark.parametrize(
    ("f_y", "c_c", "expected_mm"),
    [
        # 380*(280/280) - 2.5*40 = 280 against 300: the cover term governs.
        (420 * MPa, 40 * mm, 280.0),
        # f_s = (2/3)*500 = 333.3: 380*0.84 - 2.5*25 = 256.7 against 300*0.84 = 252.
        (500 * MPa, 25 * mm, 252.0),
    ],
    ids=["deep_cover", "stronger_steel"],
)
def test_the_slab_cap_follows_the_cover_and_the_steel_grade(f_y, c_c, expected_mm) -> None:  # type: ignore[no-untyped-def]
    """A 25 cm ACI slab: 3h = 750 and 450 mm never bind; §24.3.2 does."""
    slab = _slab(Concrete_ACI_318_19(name="H25", f_c=25 * MPa), 25 * cm, SteelBar(name="S", f_y=f_y), c_c)

    assert slab._max_bar_spacing().to("mm").magnitude == pytest.approx(expected_mm, abs=0.05)


def test_a_cirsoc_slab_reads_the_same_300_mm_twice() -> None:
    """CIRSOC 201-25 art. 7.7.2.3 already prints 300 mm; Tabla 24.3.2 agrees for
    ADN 420 with 25 mm of cover, and takes over below it: with 40 mm the
    cover term gives 380 - 100 = 280 mm."""
    concrete = Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa)
    assert _slab(concrete, 25 * cm)._max_bar_spacing().to("mm").magnitude == pytest.approx(300.0)
    assert _slab(concrete, 25 * cm, c_c=40 * mm)._max_bar_spacing().to("mm").magnitude == pytest.approx(280.0)


def test_an_en_slab_keeps_its_own_limit() -> None:
    """EN 1992-1-1 controls cracking through §7.3.3, not through Table 24.3.2:
    a 25 cm slab keeps the 400 mm of §9.3.1.1(3)."""
    slab = _slab(Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), 25 * cm, SteelBar(name="B500S", f_y=500 * MPa))

    assert slab._max_bar_spacing().to("mm").magnitude == pytest.approx(400.0)


def test_a_slab_spread_past_table_24_3_2_is_warned() -> None:
    """A 20 cm ACI slab with Ø16 every 40 cm carries 5.03 cm²/m, and the slab
    between the bars is bare: 400 mm is inside min(3h, 450) = 450 mm and past
    the 300 mm of §24.3.2, so the check says so and the report marks it."""
    slab = _slab(Concrete_ACI_318_19(name="H25", f_c=25 * MPa), 20 * cm)
    slab.set_slab_longitudinal_rebar_bot(d_b1=16 * mm, s_b1=40 * cm)
    Node(section=slab, forces=Forces(label="C1", M_y=20 * kNm)).check_flexure()

    found = {w.code: w for w in slab.warnings}
    assert set(found) == {"bar_spacing_exceeds_max"}
    assert found["bar_spacing_exceeds_max"].face == "bottom"
    assert found["bar_spacing_exceeds_max"].values["s_max"].to("mm").magnitude == pytest.approx(300.0)
    rows = slab._data_min_max_flexure
    assert rows["Check"][3] == "Bar spacing bottom"
    assert rows["Max."][3] == pytest.approx(300.0)
    assert rows["Ok?"][3] == "❌"
    # A slab's spacing row carries the cap already: no rows are added.
    assert len(rows["Check"]) == 4


def test_a_footing_takes_the_cap_with_its_own_cover() -> None:
    """§13.3.2.1 sends a one-way footing to Chapter 7, and §7.7.2.2 to Table
    24.3.2 with the footing's cover: 50 mm gives 380 - 125 = 255 mm, tighter
    than the 300 mm practice put on a mat. EN 1992-1-1 keeps the 300 mm."""
    aci = Footing(
        label="Z1",
        concrete=Concrete_ACI_318_19(name="H25", f_c=25 * MPa),
        steel_bar=ADN_420,
        width=100 * cm,
        height=60 * cm,
        c_c=50 * mm,
    )
    en = Footing(
        label="Z1",
        concrete=Concrete_EN_1992_2004(name="C25", f_c=25 * MPa),
        steel_bar=SteelBar(name="B500S", f_y=500 * MPa),
        width=100 * cm,
        height=60 * cm,
        c_c=50 * mm,
    )

    assert aci._max_bar_spacing().to("mm").magnitude == pytest.approx(255.0)
    assert en._max_bar_spacing().to("mm").magnitude == pytest.approx(300.0)
