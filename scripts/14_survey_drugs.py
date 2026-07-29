#!/usr/bin/env python
"""14_survey_drugs.py — what phenotype data exists, before committing to expand.

Run this BEFORE widening the analysis to more antibiotics. It reports, per drug:
record counts, censoring composition, testing-method mix and how many genomes
would survive the CLSI broth-microdilution restriction — so you can judge which
drugs are worth the breakpoint-curation effort.

Requires only the raw phenotype table; no breakpoints needed.

Usage:
    python scripts/14_survey_drugs.py
    python scripts/14_survey_drugs.py --min-genomes 200
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from kp_defs import BREAKPOINTS, UNVERIFIED_BREAKPOINTS

RAW = Path("data/raw/phenotypes/kp_amr_anchor.tsv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-genomes", type=int, default=100)
    a = ap.parse_args()
    if not RAW.exists():
        sys.exit(f"{RAW} not found")

    df = pd.read_csv(RAW, sep="\t", dtype=str, low_memory=False).fillna("")
    cols = {c.lower().strip(): c for c in df.columns}
    abx, gid = cols.get("antibiotic"), cols.get("genome_id")
    val, sgn = cols.get("measurement_value"), cols.get("measurement_sign")
    meth, std = cols.get("laboratory_typing_method"), cols.get("testing_standard")
    unit = cols.get("measurement_unit")

    rows = []
    for drug, g in df.groupby(df[abx].str.lower().str.strip()):
        n_rec = len(g)
        g2 = g[~g[unit].str.strip().str.lower().eq("mm")] if unit else g
        g2 = g2[pd.to_numeric(g2[val], errors="coerce").notna()]
        strat = g2
        if std:
            strat = strat[strat[std].str.strip().str.upper() == "CLSI"]
        if meth:
            strat = strat[strat[meth].str.strip().str.lower() == "broth dilution"]
        s = g2[sgn].str.strip() if sgn else pd.Series("", index=g2.index)
        verified = any(drug in BREAKPOINTS.get(k, {}) for k in BREAKPOINTS)
        listed = any(d == drug for d, _ in UNVERIFIED_BREAKPOINTS)
        rows.append(dict(
            drug=drug[:34], records=n_rec, numeric=len(g2),
            genomes_numeric=g2[gid].nunique(),
            clsi_bmd_genomes=strat[gid].nunique(),
            pct_left_cens=round(100 * s.isin(["<=", "<"]).mean(), 1) if len(g2) else 0,
            pct_right_cens=round(100 * s.isin([">=", ">"]).mean(), 1) if len(g2) else 0,
            breakpoints=("VERIFIED" if verified else
                         "listed, unverified" if listed else "ABSENT")))
    t = pd.DataFrame(rows).sort_values("clsi_bmd_genomes", ascending=False)
    pd.set_option("display.width", 200)
    print(t.to_string(index=False))

    ok = t[t["clsi_bmd_genomes"] >= a.min_genomes]
    print(f"\ndrugs with >= {a.min_genomes} CLSI broth-microdilution genomes: {len(ok)}")
    if len(ok):
        print("  " + ", ".join(ok["drug"]))
    need = ok[ok["breakpoints"] != "VERIFIED"]
    if len(need):
        print(f"\n{len(need)} of these lack VERIFIED breakpoints and cannot be")
        print("analysed until config/breakpoints.csv is completed:")
        for _, r in need.iterrows():
            print(f"  {r['drug']:<34} ({r['breakpoints']})")
    Path("results/tables").mkdir(parents=True, exist_ok=True)
    t.to_csv("results/tables/drug_survey.tsv", sep="\t", index=False)
    print("\n[out] results/tables/drug_survey.tsv")


if __name__ == "__main__":
    main()
