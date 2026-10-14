"""Design or check a grid of flexure sections with whichever mento ``PYTHONPATH`` points at.

Not part of the test suite. It measures what a change to the flexure model
moves: run it against two checkouts (``git worktree add``), then read the two
outputs with ``compare.py``. Used for the EN compression steel change
(CHANGELOG, [Unreleased]).

Usage::

    PYTHONPATH=<checkout> python sweep.py OUT.json design EN [shard nshards]
    PYTHONPATH=<checkout> python sweep.py OUT.json check EN
    PYTHONPATH=<checkout> python sweep.py OUT.json recheck:DESIGNS.json EN

``design`` designs every case of :func:`design_cases` under a positive moment
and checks the result. ``check`` checks the fixed layouts of
:func:`check_cases`. ``recheck:DESIGNS.json`` checks the layouts a ``design``
run left (another checkout's, typically) with this checkout's resistance, which
tells a design that got worse from a layout the new check rejects. The code is
``EN``, ``ACI`` or ``CIRSOC``. Shards split a run across processes, by case
index; each writes its own file.
"""

import json
import sys
import time
import traceback
import warnings

warnings.simplefilter("ignore")

from mento import (  # noqa: E402
    Concrete_ACI_318_19,
    Concrete_CIRSOC_201_25,
    Concrete_EN_1992_2004,
    Footing,
    Forces,
    Node,
    OneWaySlab,
    RectangularBeam,
    SteelBar,
)
from mento.settings import BeamSettings  # noqa: E402
from mento.units import MPa, cm, kN, kNm, mm  # noqa: E402

#: The private attributes a layout is read from and written back to. Enough to
#: rebuild the bars in another checkout; not an API.
LAYOUT = {
    "beam": [f"_{a}{i}_{f}" for f in "bt" for i in (1, 2, 3, 4) for a in ("n", "d_b")],
    "slab": [f"_{a}{i}_{f}" for f in "bt" for i in (1, 3) for a in ("d_b", "s_b")],
}


def materials(code):
    if code == "EN":
        return Concrete_EN_1992_2004(name="C25", f_c=25 * MPa), SteelBar(name="B500S", f_y=500 * MPa)
    if code == "ACI":
        return Concrete_ACI_318_19(name="H25", f_c=25 * MPa), SteelBar(name="ADN 420", f_y=420 * MPa)
    return Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa), SteelBar(name="ADN 420", f_y=420 * MPa)


def make(code, kind, b, h, c_c):
    concrete, steel = materials(code)
    if kind == "beam":
        return RectangularBeam(label="B", concrete=concrete, steel_bar=steel, width=b * cm, height=h * cm, c_c=c_c * mm)
    cls = OneWaySlab if kind == "slab" else Footing
    return cls(
        label="S",
        concrete=concrete,
        steel_bar=steel,
        width=b * cm,
        height=h * cm,
        c_c=c_c * mm,
        settings=BeamSettings(unit_system="metric"),
    )


def design_cases():
    """Beams 20-30 x 30-60 cm, and one-way slabs and footings 100 x 15-30 cm."""
    cases = []
    for b in (20, 25, 30):
        for h in (30, 40, 50, 60):
            for c_c in (25, 40):
                for M in (40, 80, 120, 160, 200, 250, 300, 400):
                    cases.append(("beam", b, h, c_c, M))
    for kind in ("slab", "footing"):
        for h in (15, 20, 25, 30):
            for c_c in (20, 30, 50, 70):
                for M in (20, 40, 60, 80, 100, 125, 150, 200, 300):
                    cases.append((kind, 100, h, c_c, M))
    return cases


def check_cases():
    """Fixed layouts, with and without top bars, on the same kinds of section."""
    cases = []
    for kind in ("slab", "footing"):
        for h in (15, 20, 25, 30):
            for c_c in (20, 30, 50, 70):
                for d_b in (12, 16, 20, 25):
                    for s in (10, 15):
                        for d_t, s_t in ((0, 0), (10, 15), (12, 10), (16, 10), (20, 10), (25, 10)):
                            cases.append((kind, 100, h, c_c, ("slab", d_b, s, d_t, s_t)))
    for b in (20, 30):
        for h in (30, 40, 50, 60):
            for c_c in (25, 40):
                for n_b, d_b in ((3, 16), (3, 20), (4, 20), (3, 25), (4, 25), (6, 20)):
                    for n_t, d_t in ((0, 0), (2, 10), (2, 16), (2, 20), (3, 20), (2, 25)):
                        cases.append(("beam", b, h, c_c, ("beam", n_b, d_b, n_t, d_t)))
    return cases


def raw_layout(section):
    keys = LAYOUT["slab" if isinstance(section, OneWaySlab) else "beam"]
    layout = {}
    for key in keys:
        value = getattr(section, key)
        layout[key] = value if isinstance(value, int) else float(value.to("mm").magnitude)
    return layout


def apply_layout(section, layout):
    for key, value in layout.items():
        setattr(section, key, value if key.startswith("_n") else value * mm)
    if isinstance(section, OneWaySlab):
        section._calculate_longitudinal_rebars()
    section._update_longitudinal_rebar_attributes()
    section._drop_results()


def set_layout(section, layout):
    kind, a, b, c, d = layout
    if kind == "slab":
        section.set_slab_longitudinal_rebar_bot(d_b1=a * mm, s_b1=b * cm)
        if c:
            section.set_slab_longitudinal_rebar_top(d_b1=c * mm, s_b1=d * cm)
    else:
        section.set_longitudinal_rebar_bot(a, b * mm)
        if c:
            section.set_longitudinal_rebar_top(c, d * mm)


def face(check):
    M_capacity = check.M_capacity
    return {
        "DCR": round(float(check.DCR), 4),
        "M_cap": round(float(M_capacity.to("kN*m").magnitude), 3) if M_capacity is not None else None,
    }


def summary(section, node):
    check = section.flexure_checks[0]
    return {
        "layout": raw_layout(section),
        "bot": str(section.reinforcement.bottom),
        "top": str(section.reinforcement.top),
        "check_bot": face(check.bottom),
        "check_top": face(check.top),
        "warnings": sorted({warning.code for warning in node.warnings}),
    }


def run_case(code, mode, case):
    kind, b, h, c_c, extra = case
    record = {"case": [kind, b, h, c_c, list(extra) if isinstance(extra, tuple) else extra]}
    started = time.time()
    try:
        section = make(code, kind, b, h, c_c)
        if mode == "design":
            node = Node(section=section, forces=[Forces(label="M", M_y=extra * kNm, V_z=0 * kN)])
            node.design_flexure()
        else:
            set_layout(section, extra)
            # M_Rd does not depend on M, only on its sign.
            node = Node(section=section, forces=[Forces(label="M", M_y=1 * kNm, V_z=0 * kN)])
        node.check_flexure()
        record.update(summary(section, node))
    except Exception as exc:  # noqa: BLE001 -- a crash is a result here
        record["error"] = f"{type(exc).__name__}: {exc}"[:200]
        record["traceback"] = traceback.format_exc()[-600:]
    record["t"] = round(time.time() - started, 2)
    return record


def recheck(code, source):
    """Check the layouts a design run left, under this checkout's resistance."""
    records = []
    for previous in json.load(open(source, encoding="utf-8")):
        if "layout" not in previous:
            continue
        kind, b, h, c_c, M = previous["case"]
        record = {"case": previous["case"]}
        try:
            section = make(code, kind, b, h, c_c)
            apply_layout(section, previous["layout"])
            node = Node(section=section, forces=[Forces(label="M", M_y=M * kNm, V_z=0 * kN)])
            node.check_flexure()
            record.update(summary(section, node))
            record["bot"], record["top"] = previous["bot"], previous["top"]
        except Exception as exc:  # noqa: BLE001
            record["error"] = f"{type(exc).__name__}: {exc}"[:200]
        records.append(record)
    return records


def main(argv):
    out, mode, code = argv[1], argv[2], argv[3]
    shard, nshards = (int(argv[4]), int(argv[5])) if len(argv) > 5 else (0, 1)
    if mode.startswith("recheck:"):
        records = recheck(code, mode.split(":", 1)[1])
    else:
        cases = design_cases() if mode == "design" else check_cases()
        records = [run_case(code, mode, case) for i, case in enumerate(cases) if i % nshards == shard]
    with open(out, "w", encoding="utf-8") as f:
        json.dump(records, f)
    print(out, len(records))


if __name__ == "__main__":
    main(sys.argv)
