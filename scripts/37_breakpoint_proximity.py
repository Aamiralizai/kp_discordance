#!/usr/bin/env python3
"""37_breakpoint_proximity.py — how much of the discordance sits on the breakpoint?

The headline estimate counts a carbapenemase-positive isolate as discordant when
its minimum inhibitory concentration is at or below the susceptible breakpoint.
Broth microdilution is reproducible to about one doubling dilution, so an
isolate reported at exactly the breakpoint is one dilution from the next
category. If most discordant isolates sit there, part of the headline is assay
resolution rather than biology, and a reader has to be told.

This reports the distribution of the reported upper bound among discordant
isolates, and repeats the headline under a definition one doubling dilution
stricter, in which susceptibility requires the whole interval to lie at or below
half the breakpoint. The stricter definition is not proposed as a better
breakpoint. It is a bound on how much of the estimate survives when
breakpoint-adjacent isolates are removed.

The same is done per enzyme family, because the question that matters is not
only how large the headline is but whether the family ordering depends on the
breakpoint-adjacent isolates.

USAGE
    python3 scripts/37_breakpoint_proximity.py [--drug meropenem] [--tag _clsi_bmd]
"""
import argparse
import csv
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
SUP = os.path.join(ROOT, "results", "supplementary")


def wilson(k, n):
    if n == 0:
        return 0.0, 0.0, 0.0
    p = k / n
    z = 1.959963984540054
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * p, 100 * max(0.0, c - h), 100 * min(1.0, c + h)


def bound(v, default):
    if v in ("", None):
        return default
    if v == "inf":
        return float("inf")
    try:
        return float(v)
    except ValueError:
        return default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="meropenem")
    ap.add_argument("--tag", default="_clsi_bmd")
    ap.add_argument("--breakpoint", type=float, default=1.0,
                    help="susceptible breakpoint in mg/L, CLSI meropenem is 1")
    a = ap.parse_args()

    src = os.path.join(TAB, f"recordlevel_discordance_{a.drug}{a.tag}.tsv")
    with open(src, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["call"] in ("S", "R")]
    pos = [r for r in rows if r["carb_pos"] == "True"]
    sus = [r for r in pos if r["call"] == "S"]
    strict_cut = a.breakpoint / 2.0

    # distribution of the reported upper bound among the discordant isolates.
    # An interval whose two bounds are equal is an exact reading at that
    # concentration, which is a stronger statement than an upper bound there,
    # so the two are counted separately.
    dist = {}
    exact = {}
    for r in sus:
        hi = bound(r["hi"], float("inf"))
        lo = bound(r["lo"], 0.0)
        dist[hi] = dist.get(hi, 0) + 1
        if lo == hi and r["lo_inc"] == "True" and r["hi_inc"] == "True":
            exact[hi] = exact.get(hi, 0) + 1
    rows_out = []
    cum = 0
    for hi in sorted(dist):
        cum += dist[hi]
        rows_out.append({
            "drug": a.drug, "stratum": a.tag.lstrip("_"),
            "mic_upper_mgL": ("inf" if hi == float("inf") else f"{hi:g}"),
            "n_isolates": dist[hi],
            "pct_of_discordant": round(100 * dist[hi] / len(sus), 1),
            "cumulative_pct": round(100 * cum / len(sus), 1),
            "n_exact_reading": exact.get(hi, 0),
            "pct_exact_reading": round(100 * exact.get(hi, 0) / len(sus), 1),
            "at_or_above_breakpoint": hi >= a.breakpoint,
        })
    out1 = os.path.join(TAB, f"breakpoint_proximity_{a.drug}{a.tag}.tsv")
    with open(out1, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    at_bp = sum(n for hi, n in dist.items() if hi >= a.breakpoint)
    print(f"  {a.drug}{a.tag}: {len(sus)} discordant of {len(pos)} carbapenemase-positive")
    print(f"  upper bound at or above the breakpoint of {a.breakpoint:g} mg/L: "
          f"{at_bp}/{len(sus)} = {100 * at_bp / len(sus):.1f} percent")
    ex_bp = sum(n for hi, n in exact.items() if hi >= a.breakpoint)
    print(f"  exact reading at or above that breakpoint: "
          f"{ex_bp}/{len(sus)} = {100 * ex_bp / len(sus):.1f} percent")

    # the headline, and the headline one dilution stricter, overall and by family
    fams = ["KPC", "NDM", "OXA-48-like", "VIM"]
    sens = []
    def block(label, subset):
        n = len(subset)
        k1 = sum(1 for r in subset if r["call"] == "S")
        k2 = sum(1 for r in subset
                 if r["call"] == "S" and bound(r["hi"], float("inf")) <= strict_cut)
        p1, l1, h1 = wilson(k1, n)
        p2, l2, h2 = wilson(k2, n)
        sens.append({
            "drug": a.drug, "stratum": a.tag.lstrip("_"), "group": label, "n": n,
            "susceptible_as_published": k1,
            "pct_as_published": round(p1, 1),
            "ci_lo_as_published": round(l1, 1), "ci_hi_as_published": round(h1, 1),
            "susceptible_strict": k2,
            "pct_strict": round(p2, 1),
            "ci_lo_strict": round(l2, 1), "ci_hi_strict": round(h2, 1),
            "strict_definition": f"MIC <= {strict_cut:g} mg/L",
        })
    block("all carbapenemase-positive", pos)
    for f in fams:
        g = [r for r in pos if r["carb_family"] == f]
        if g:
            block(f, g)
    out2 = os.path.join(TAB, f"breakpoint_sensitivity_{a.drug}{a.tag}.tsv")
    with open(out2, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sens[0]), delimiter="\t")
        w.writeheader()
        w.writerows(sens)

    print(f"\n  {'group':28s}{'as published':>22s}{'one dilution stricter':>26s}")
    for s in sens:
        print(f"  {s['group']:28s}"
              f"{s['pct_as_published']:>8.1f}% ({s['ci_lo_as_published']:.1f}-{s['ci_hi_as_published']:.1f})"
              f"{s['pct_strict']:>12.1f}% ({s['ci_lo_strict']:.1f}-{s['ci_hi_strict']:.1f})")
    print(f"\n[out] {out1}")
    print(f"[out] {out2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
