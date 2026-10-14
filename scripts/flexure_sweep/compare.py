"""Compare two runs of ``sweep.py``: what passed, what fails now, and why.

Usage::

    python compare.py OUT_DIR BEFORE AFTER MODE CODE [--recheck BEFORE_ON_AFTER]

Reads ``OUT_DIR/{BEFORE,AFTER}_{MODE}_{CODE}_*.json``, the shard files of two
sweeps, and compares them case by case. For designs: pass/fail transitions,
whether the layout or the warnings moved, and the failures no shortfall
warning explains. With ``--recheck``, a third run (``sweep.py``
``recheck:`` of BEFORE's designs on AFTER's checkout) tells, for each design
that passed and fails now, whether the old layout still passes the new check:
if it does, the design got worse, not the section. For checks: how M_Rd moved.
"""

import argparse
import collections
import glob
import json
import statistics

EXPLAINED = {"As_below_required", "section_too_small_for_moment", "As_above_max", "not_tension_controlled"}


def load(out_dir, tree, mode, code):
    records = {}
    for path in sorted(glob.glob(f"{out_dir}/{tree}_{mode}_{code}_*.json")):
        for record in json.load(open(path, encoding="utf-8")):
            records[json.dumps(record["case"])] = record
    return records


def passes(record):
    return record is not None and "error" not in record and record["check_bot"]["DCR"] <= 1 + 1e-9


def state(record):
    return "error" if "error" in record else "pass" if passes(record) else "fail"


def silent(records):
    """Designs past DCR 1 with no warning that says so."""
    return [
        key
        for key, record in records.items()
        if "error" not in record and record["check_bot"]["DCR"] > 1 + 1e-9 and not EXPLAINED & set(record["warnings"])
    ]


def compare_designs(before, after, rechecked):
    transitions = collections.Counter()
    moved = collections.Counter()
    for key in before.keys() & after.keys():
        a, b = before[key], after[key]
        transitions[(state(a), state(b))] += 1
        if "error" in a or "error" in b:
            continue
        if (a["bot"], a["top"]) != (b["bot"], b["top"]):
            moved["layout"] += 1
        elif a["check_bot"]["DCR"] != b["check_bot"]["DCR"]:
            moved["DCR, same layout"] += 1
        elif a["warnings"] != b["warnings"]:
            moved["warnings only"] += 1
    print("transitions:", {f"{x} -> {y}": n for (x, y), n in sorted(transitions.items())})
    print("changes:", dict(moved))
    print("silent failures before:", len(silent(before)), "after:", len(silent(after)))
    for name, records in (("before", before), ("after", after)):
        for key in silent(records):
            print(f"  silent ({name}): {key}")

    regressions = sorted(key for key in before.keys() & after.keys() if passes(before[key]) and not passes(after[key]))
    print(f"\npass -> fail: {len(regressions)}")
    print("| case | M | before | DCR | before, rechecked | after | DCR |")
    print("|---|---|---|---|---|---|---|")
    misses = 0
    for key in regressions:
        a, b, c = before[key], after[key], rechecked.get(key) if rechecked else None
        kind, width, h, c_c, M = json.loads(key)
        miss = passes(c)
        misses += miss
        rechecked_DCR = f"{c['check_bot']['DCR']:.3f}" if c and "error" not in c else "-"
        print(
            f"| {kind} {width}x{h} c{c_c} | {M} | {a['bot']} / {a['top']} | {a['check_bot']['DCR']:.3f} "
            f"| {rechecked_DCR} | {b['bot']} / {b['top']} | {b['check_bot']['DCR']:.3f} |"
            + (" design got worse" if miss else "")
        )
    if rechecked:
        print(f"\nold layouts that still pass the new check: {misses} of {len(regressions)}")


def compare_checks(before, after):
    raised, lowered, same = [], [], 0
    for key in before.keys() & after.keys():
        a, b = before[key]["check_bot"]["M_cap"], after[key]["check_bot"]["M_cap"]
        if abs(a - b) < 1e-6:
            same += 1
        elif b > a:
            raised.append((key, a, b))
        else:
            lowered.append((b / a, key, a, b))
    print(f"M_Rd unchanged: {same}, higher: {len(raised)}, lower: {len(lowered)}")
    if lowered:
        ratios = sorted(ratio for ratio, *_ in lowered)
        print("lower, new/old: min %.3f, median %.3f, max %.3f" % (ratios[0], statistics.median(ratios), ratios[-1]))
    for key, a, b in sorted(raised, key=lambda item: item[1])[:20]:
        print(f"  higher: {key} {a} -> {b}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out_dir")
    parser.add_argument("before")
    parser.add_argument("after")
    parser.add_argument("mode", choices=["design", "check"])
    parser.add_argument("code", choices=["EN", "ACI", "CIRSOC"])
    parser.add_argument("--recheck", help="tree name of BEFORE's designs rechecked on AFTER")
    args = parser.parse_args()
    before = load(args.out_dir, args.before, args.mode, args.code)
    after = load(args.out_dir, args.after, args.mode, args.code)
    print(f"== {args.before} -> {args.after}, {args.mode} {args.code}: {len(before.keys() & after.keys())} cases")
    if args.mode == "check":
        compare_checks(before, after)
    else:
        rechecked = load(args.out_dir, args.recheck, "design", args.code) if args.recheck else {}
        compare_designs(before, after, rechecked)


if __name__ == "__main__":
    main()
