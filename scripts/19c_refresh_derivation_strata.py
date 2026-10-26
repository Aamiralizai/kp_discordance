#!/usr/bin/env python3
"""19c_refresh_derivation_strata.py — recompute the derivation-set columns of
external_stratum_comparison.tsv from the current record-level table.

WHY THIS EXISTS
19_external_validation.py needs the NCBI Pathogen Detection snapshot, which is
181 MB and fetched separately. When the derivation genotypes change but the
external cohort and the fitted rules do not, rerunning the whole external step
is unnecessary, but the derivation-side columns of the comparison table would
otherwise keep stale values.

This script refreshes only those columns, and it refuses to run if any decision
rule changed its predicted call, because in that case the external predictions
would change too and the full external step must be rerun.

    python3 scripts/19c_refresh_derivation_strata.py [--prev results_round5_prev]
"""
import argparse
import csv
import math
from pathlib import Path

REC = Path("results/tables/recordlevel_discordance_meropenem_clsi_bmd.tsv")
CMP = Path("results/tables/external_stratum_comparison.tsv")
RULES = Path("results/tables/rules_decision_list.tsv")


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return 100 * p, (c - m) / d * 100, (c + m) / d * 100


def rows(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def porin_class(v):
    if not v or v == "none":
        return "intact"
    if v == "ompK35_loss":
        return "ompK35_only"
    return "ompK36_disrupted"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prev", help="previous results tree, to verify no rule flipped")
    a = ap.parse_args()

    if a.prev:
        old = {(r["tier"], r["rule"]): r["then_predict"]
               for r in rows(Path(a.prev) / "tables" / "rules_decision_list.tsv")}
        new = {(r["tier"], r["rule"]): r["then_predict"] for r in rows(RULES)}
        flipped = [k for k in new if k in old and new[k] != old[k]]
        if flipped or set(new) != set(old):
            raise SystemExit(
                "Decision rules changed, so the external predictions change too. "
                "Rerun scripts/19_external_validation.py in full.")

    rec = [r for r in rows(REC) if r["call"] in ("S", "R")]
    strata = {}
    for r in rec:
        fam = r["carb_family"] or "none"
        if r["carb_pos"].lower() not in ("true", "1"):
            fam = "none"
        key = (fam, porin_class(r["porin"]))
        n, k = strata.get(key, (0, 0))
        strata[key] = (n + 1, k + (r["call"] == "R"))

    out, refreshed = [], 0
    for row in rows(CMP):
        key = (row["family"], row["porin"])
        if key in strata:
            n, k = strata[key]
            p, lo, hi = wilson(k, n)
            row["n_derivation"] = f"{float(n):.1f}"
            row["n_resistant_derivation"] = f"{float(k):.1f}"
            row["resistance_pct_derivation"] = f"{p:.1f}"
            row["ci_lower_derivation"] = f"{lo:.1f}"
            row["ci_upper_derivation"] = f"{hi:.1f}"
            if row["resistance_pct_external"]:
                row["delta_pct"] = f"{float(row['resistance_pct_external']) - p:.1f}"
                row["ci_overlap"] = str(
                    not (hi < float(row["ci_lower_external"])
                         or float(row["ci_upper_external"]) < lo))
            refreshed += 1
        out.append(row)

    with open(CMP, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print(f"refreshed the derivation columns of {refreshed} strata in {CMP}")


if __name__ == "__main__":
    main()
