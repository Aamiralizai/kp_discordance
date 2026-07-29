#!/usr/bin/env python
"""14_results_manifest.py — single source of truth for every reported number.

Writes results/manuscript_results.json containing every headline value, read
from the generated record-level tables. Reconcile the abstract, Results,
Discussion, tables and figures against THIS FILE, not against remembered
values or earlier drafts.

Usage:
    python scripts/14_results_manifest.py
"""
import json
import math
import sys
from pathlib import Path

import pandas as pd

TAB = Path("results/tables")
OUT = Path("results/manuscript_results.json")


def wilson(k, n):
    if n == 0:
        return (None, None, None)
    p, z = k / n, 1.959964
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (round(100 * p, 1), round(100 * max(0, c - h), 1),
            round(100 * min(1, c + h), 1))


def summarise(drug, tag, label):
    f = TAB / f"recordlevel_discordance_{drug}{tag}.tsv"
    if not f.exists():
        return None
    d = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
    d["carb_pos"] = d["carb_pos"].astype(str).str.lower().isin({"true", "1"})
    b = d[d["call"].isin(["S", "R"])]
    pos, neg = b[b["carb_pos"]], b[~b["carb_pos"]]
    fs = int((pos["call"] == "S").sum()); fr = int((neg["call"] == "R").sum())
    tp = int((pos["call"] == "R").sum()); tn = int((neg["call"] == "S").sum())
    out = {
        "analysis": label, "drug": drug, "tag": tag,
        "linked_genomes": len(d),
        "susceptible": int((d["call"] == "S").sum()),
        "intermediate": int((d["call"] == "I").sum()),
        "resistant": int((d["call"] == "R").sum()),
        "unresolvable": int((d["call"] == "ambiguous").sum()),
        "binary_set": len(b),
        "carbapenemase_positive_n": len(pos),
        "carbapenemase_positive_susceptible_n": fs,
        "carbapenemase_positive_susceptible_pct_ci": wilson(fs, len(pos)),
        "carbapenemase_negative_n": len(neg),
        "carbapenemase_negative_resistant_n": fr,
        "carbapenemase_negative_resistant_pct_ci": wilson(fr, len(neg)),
        "PPV": round(tp / (tp + fs), 3) if tp + fs else None,
        "NPV": round(tn / (tn + fr), 3) if tn + fr else None,
        "sensitivity": round(tp / (tp + fr), 3) if tp + fr else None,
        "specificity": round(tn / (tn + fs), 3) if tn + fs else None,
        "by_family": {},
    }
    for fam, g in pos.groupby("carb_family"):
        if len(g) < 10:
            continue
        k = int((g["call"] == "S").sum())
        out["by_family"][fam] = {"n": len(g), "susceptible": k,
                                 "pct_ci": wilson(k, len(g))}
    ff = TAB / f"recordlevel_fisher_{drug}{tag}.tsv"
    if ff.exists():
        fisher = pd.read_csv(ff, sep="\t")
        out["fisher_intact_vs_ompK36"] = fisher.to_dict("records")
    return out


def main():
    res = {"note": ("Every value here is computed from the record-level tables. "
                    "Reconcile all manuscript text against this file."),
           "analyses": []}
    for drug, tag, label in (
            ("meropenem", "_clsi_bmd", "primary (CLSI, broth dilution)"),
            ("imipenem", "_clsi_bmd", "imipenem (CLSI, broth dilution)"),
            ("meropenem", "_clsi_bmd_eucast", "meropenem (EUCAST breakpoints)"),
            ("imipenem", "_clsi_bmd_eucast", "imipenem (EUCAST breakpoints)"),
            ("meropenem", "_clsi_bmd_dedup", "meropenem deduplicated"),
            ("meropenem", "_unstratified", "meropenem unstratified")):
        s = summarise(drug, tag, label)
        if s:
            res["analyses"].append(s)
        else:
            res.setdefault("missing", []).append(f"{drug}{tag}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=2))
    print(f"[out] {OUT}  ({len(res['analyses'])} analyses)")
    for a in res["analyses"]:
        p = a["carbapenemase_positive_susceptible_pct_ci"]
        print(f"  {a['analysis']:<38} binary {a['binary_set']:>6,}  "
              f"gene+ susceptible {p[0]}% [{p[1]}-{p[2]}]")
    if res.get("missing"):
        print("\n  NOT YET GENERATED: " + ", ".join(res["missing"]))
        print("  (manuscript must not cite values for these)")


if __name__ == "__main__":
    main()
