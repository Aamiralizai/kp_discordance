#!/usr/bin/env python3
"""38_provenance_stratified_family.py — does the family ordering survive provenance?

The paper concedes that discordance varies by reporting year and leaves the
concession there. That is the wrong place to stop, because the families are
very unevenly distributed across reporting years. Most OXA-48-like and VIM
isolates come from one year and most KPC isolates from another, so a reader is
entitled to suspect that the family ordering is a provenance artefact.

Two things are established here. First, reporting year and instrument platform
are not two confounders but one: they are perfectly nested in this collection,
so a model containing both is not identifiable and the two heterogeneity tests
the paper reports are the same test twice.

Second, the family ordering is compared within the single provenance stratum
that contains all four families at usable size. A within-stratum comparison
needs no model and is not disturbed by the imbalance, so it answers the
objection directly.

USAGE
    python3 scripts/38_provenance_stratified_family.py [--drug meropenem] [--tag _clsi_bmd]
"""
import argparse
import collections
import csv
import math
import os
import sys

from scipy.stats import fisher_exact

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
FAMS = ["KPC", "NDM", "OXA-48-like", "VIM"]
MIN_CELL = 10


def wilson(k, n):
    if n == 0:
        return 0.0, 0.0, 0.0
    p = k / n
    z = 1.959963984540054
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * p, 100 * max(0.0, c - h), 100 * min(1.0, c + h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="meropenem")
    ap.add_argument("--tag", default="_clsi_bmd")
    a = ap.parse_args()

    src = os.path.join(TAB, f"recordlevel_discordance_{a.drug}{a.tag}.tsv")
    with open(src, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["call"] in ("S", "R")]
    pos = [r for r in rows if r["carb_pos"] == "True"]
    year = lambda r: (r["year"] or "<not reported>")          # noqa: E731
    plat = lambda r: (r["platform"] or "<not reported>")      # noqa: E731

    # 1. are year and platform the same variable?
    pairs = collections.Counter((year(r), plat(r)) for r in pos)
    reported = {y for (y, _) in pairs if y != "<not reported>"}
    nested_reported = all(
        len({p for (y, p) in pairs if y == yy}) == 1 for yy in reported)
    nested = all(
        len({p for (y, p) in pairs if y == yy}) == 1
        for yy in {y for (y, _) in pairs})
    nest_rows = [{"drug": a.drug, "stratum": a.tag.lstrip("_"),
                  "reporting_year": y, "platform": p, "n": n}
                 for (y, p), n in sorted(pairs.items())]
    nest_rows.append({"drug": a.drug, "stratum": a.tag.lstrip("_"),
                      "reporting_year": "ALL", "platform": "ALL",
                      "n": f"platform nested within every reported year: "
                           f"{nested_reported}; nested within all year groups "
                           f"including not-reported: {nested}"})
    out0 = os.path.join(TAB, f"provenance_nesting_{a.drug}{a.tag}.tsv")
    with open(out0, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["drug", "stratum", "reporting_year",
                                           "platform", "n"], delimiter="\t")
        w.writeheader()
        w.writerows(nest_rows)
    print(f"  platform nested within every REPORTED year: {nested_reported}")
    print(f"  nested across all year groups including not-reported: {nested}")
    for (y, pl), n in sorted(pairs.items()):
        print(f"    {y:>16s}  {pl:<16s} {n:5d}")

    # 2. family composition by provenance stratum
    comp = []
    for f in FAMS:
        g = [r for r in pos if r["carb_family"] == f]
        c = collections.Counter(year(r) for r in g)
        top, topn = (c.most_common(1)[0] if c else ("", 0))
        comp.append({"drug": a.drug, "family": f, "n": len(g),
                     "largest_stratum": top,
                     "n_in_largest": topn,
                     "pct_in_largest": round(100 * topn / len(g), 1) if g else ""})
    out1 = os.path.join(TAB, f"family_by_provenance_{a.drug}{a.tag}.tsv")
    with open(out1, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(comp[0]), delimiter="\t")
        w.writeheader()
        w.writerows(comp)
    print("\n  family concentration in a single reporting year")
    for c in comp:
        print(f"    {c['family']:14s} n={c['n']:4d}  "
              f"{c['pct_in_largest']:5.1f}% in {c['largest_stratum']}")

    # 3. the ordering within every stratum that supports it
    strata = sorted({year(r) for r in pos})
    out = []
    for s in strata:
        sub = [r for r in pos if year(r) == s]
        usable = [f for f in FAMS
                  if len([r for r in sub if r["carb_family"] == f]) >= MIN_CELL]
        if len(usable) < 2:
            continue
        ref = "KPC" if "KPC" in usable else usable[0]
        gref = [r for r in sub if r["carb_family"] == ref]
        kref = sum(1 for r in gref if r["call"] == "S")
        for f in usable:
            g = [r for r in sub if r["carb_family"] == f]
            k = sum(1 for r in g if r["call"] == "S")
            p, lo, hi = wilson(k, len(g))
            if f == ref:
                orr, pv = ("reference", "reference")
            else:
                o, pp = fisher_exact([[k, len(g) - k], [kref, len(gref) - kref]])
                # An odds ratio is not estimable when a cell of the 2x2 table
                # is empty. Reporting the floating-point artefact of the
                # attempt, nan for an empty cell in both arms and inf for one
                # in the reference arm only, states a result where there is
                # none, so both become the same explicit label.
                if math.isnan(o) or math.isinf(o):
                    orr = "not estimable"
                else:
                    orr = round(o, 2)
                pv = "not estimable" if math.isnan(pp) else f"{pp:.2g}"
            out.append({"drug": a.drug, "stratum": a.tag.lstrip("_"),
                        "reporting_year": s, "family": f, "n": len(g),
                        "susceptible": k, "pct": round(p, 1),
                        "ci_lo": round(lo, 1), "ci_hi": round(hi, 1),
                        "reference_family": ref,
                        "odds_ratio_vs_reference": orr,
                        "fisher_p_vs_reference": pv,
                        "families_comparable_here": len(usable)})
    out2 = os.path.join(TAB, f"family_within_provenance_{a.drug}{a.tag}.tsv")
    with open(out2, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    print("\n  susceptibility by family within each provenance stratum")
    cur = None
    for r in out:
        if r["reporting_year"] != cur:
            cur = r["reporting_year"]
            print(f"\n    reporting year {cur}  "
                  f"({r['families_comparable_here']} families comparable)")
        print(f"      {r['family']:14s} {r['susceptible']:3d}/{r['n']:<4d} "
              f"{r['pct']:5.1f}% ({r['ci_lo']:.1f}-{r['ci_hi']:.1f})   "
              f"OR vs {r['reference_family']} {r['odds_ratio_vs_reference']}  "
              f"p {r['fisher_p_vs_reference']}")
    print(f"\n[out] {out0}\n[out] {out1}\n[out] {out2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
