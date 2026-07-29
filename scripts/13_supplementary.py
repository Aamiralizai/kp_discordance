#!/usr/bin/env python
"""13_supplementary.py — build supplementary tables from analysis outputs.

Produces:
  Supplementary Table S1  genome-level dataset (accession, ST, genotype, MIC, call)
  Supplementary Table S2  discordance by enzyme family, stratified and unstratified
  Supplementary Table S3  porin-class breakdown within enzyme family
  Supplementary Table S4  sensitivity analyses summary
  Supplementary Table S5  curated carbapenemase kinetic constants
  Supplementary Table S6  provenance composition of the dataset

Outputs both a multi-sheet .xlsx (for submission) and .tsv files.

Usage:
    python scripts/13_supplementary.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TAB = Path("results/tables")
SUPP = Path("results/supplementary")
KLEB = Path("data/processed/kleborate_klebsiella_pneumo_complex_output.tsv")
RAW = Path("data/raw/phenotypes/kp_amr_anchor.tsv")
KIN = Path("data/raw/enzyme_kinetics/kinetic_constants_curated.csv")


def load_disc(drug, tag=""):
    f = TAB / f"recordlevel_discordance_{drug}{tag}.tsv"
    if not f.exists():
        return None   # legacy discordance_*.tsv is INVALID and never used
    d = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
    for c in ("mic", "lo", "hi"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d.attrs["n_linked"] = len(d)
    for c in ("carb_pos", "lc", "rc"):
        if c in d.columns:
            d[c] = d[c].astype(str).str.lower().isin({"true", "1"})
    return d


def wilson(k, n):
    import math
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n; z = 1.959964
    dd = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / dd
    h = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / dd
    return (100 * p, 100 * max(0, c - h), 100 * min(1, c + h))


def main():
    SUPP.mkdir(parents=True, exist_ok=True)
    sheets = {}

    # ---------- S1: genome-level dataset ----------
    mero = load_disc("meropenem", "_clsi_bmd")
    imi = load_disc("imipenem", "_clsi_bmd")
    if mero is None:
        sys.exit("Run 08_discordance_analysis.py --tag _clsi_bmd first")

    s1 = mero.rename(columns={
        "genome_id": "BVBRC_genome_id", "ST": "sequence_type",
        "carb_family": "carbapenemase_family", "carb_allele": "carbapenemase_allele",
        "carb_pos": "carbapenemase_present", "porin": "porin_class",
        "mic": "meropenem_MIC_mgL", "lo": "meropenem_MIC_lower_mgL",
        "hi": "meropenem_MIC_upper_mgL", "n_records": "n_phenotype_records",
        "method": "testing_method", "standard": "interpretive_standard",
        "platform": "testing_platform", "year": "standard_year",
        "call": "meropenem_interpretation"})
    keep = [c for c in ["BVBRC_genome_id", "sequence_type", "carbapenemase_present",
                        "carbapenemase_family", "carbapenemase_allele", "porin_class",
                        "meropenem_MIC_lower_mgL", "meropenem_MIC_upper_mgL",
                        "n_phenotype_records", "testing_method",
                        "interpretive_standard", "testing_platform", "standard_year",
                        "meropenem_interpretation"] if c in s1.columns]
    s1 = s1[keep]
    if imi is not None:
        icols = ["genome_id", "call"] + [c for c in ("lo", "hi") if c in imi.columns]
        add = imi[icols].rename(columns={
            "genome_id": "BVBRC_genome_id", "lo": "imipenem_MIC_lower_mgL",
            "hi": "imipenem_MIC_upper_mgL", "call": "imipenem_interpretation"})
        s1 = s1.merge(add, on="BVBRC_genome_id", how="left")
    sheets["S1_genome_level_data"] = s1.sort_values("BVBRC_genome_id")

    # ---------- S2: discordance by family, all strata ----------
    def binary(df):
        """Rates are computed on S/R only; I and unresolved are not denominators."""
        return df[df["call"].isin(["S", "R"])].copy()

    rows = []
    for drug in ("meropenem", "imipenem"):
        for tag, label in (("", "unstratified"),
                           ("_clsi_bmd", "CLSI broth microdilution"),
                           ("_clsi_bmd_dedup", "CLSI BMD, deduplicated")):
            d = load_disc(drug, tag)
            if d is None:
                continue
            d = binary(d)
            pos = d[d["carb_pos"]]
            for fam, g in pos.groupby("carb_family"):
                if len(g) < 5:
                    continue
                k = int((g["call"] == "S").sum())
                p, lo, hi = wilson(k, len(g))
                rows.append(dict(drug=drug, stratum=label, carbapenemase_family=fam,
                                 n=len(g), n_susceptible=k,
                                 pct_susceptible=round(p, 1),
                                 ci_lower=round(lo, 1), ci_upper=round(hi, 1)))
    sheets["S2_discordance_by_family"] = pd.DataFrame(rows)

    # ---------- S3: porin breakdown ----------
    rows = []
    for drug in ("meropenem", "imipenem"):
        d = load_disc(drug, "_clsi_bmd")
        if d is None:
            continue
        d = d[d["call"].isin(["S", "R"])].copy()
        pos = d[d["carb_pos"]]
        for (fam, pr), g in pos.groupby(["carb_family", "porin"]):
            if len(g) < 5:
                continue
            k = int((g["call"] == "S").sum())
            p, lo, hi = wilson(k, len(g))
            rows.append(dict(drug=drug, carbapenemase_family=fam, porin_class=pr,
                             n=len(g), n_susceptible=k,
                             pct_susceptible=round(p, 1),
                             ci_lower=round(lo, 1), ci_upper=round(hi, 1)))
    sheets["S3_porin_breakdown"] = pd.DataFrame(rows)

    # ---------- S4: sensitivity summary ----------
    s4 = []
    for f, label in ((TAB / "sensitivity_breakpoints_meropenem.tsv", "breakpoint standard (meropenem)"),
                     (TAB / "sensitivity_breakpoints_imipenem.tsv", "breakpoint standard (imipenem)"),
                     (TAB / "sensitivity_fdr_meropenem.tsv", "FDR correction (meropenem)"),
                     (TAB / "sensitivity_fdr_imipenem.tsv", "FDR correction (imipenem)")):
        if f.exists():
            df = pd.read_csv(f, sep="\t")
            df.insert(0, "analysis", label)
            s4.append(df)
    if s4:
        sheets["S4_sensitivity"] = pd.concat(s4, ignore_index=True)

    # ---------- S5: kinetic constants ----------
    if KIN.exists():
        sheets["S5_kinetic_constants"] = pd.read_csv(KIN)

    # ---------- S6: provenance composition ----------
    if RAW.exists():
        raw = pd.read_csv(RAW, sep="\t", dtype=str, low_memory=False).fillna("")
        abc = next((c for c in raw.columns if c.lower().strip() == "antibiotic"), None)
        cols = [c for c in raw.columns if any(
            k in c.lower() for k in ["laboratory_typing", "testing_standard"])]
        rows = []
        if abc and cols:
            for drug in ("meropenem", "imipenem"):
                sub = raw[raw[abc].str.lower().str.strip() == drug]
                for c in cols:
                    vc = sub[c].replace("", "<not reported>").value_counts()
                    for lvl, n in vc.items():
                        rows.append(dict(drug=drug, field=c, level=lvl, n_records=int(n)))
        sheets["S6_provenance"] = pd.DataFrame(rows)

    # ---------- S7-S10: generated rather than static ----------
    ext = Path("results/tables")
    if (ext / "extB_incremental_value.tsv").exists():
        sheets["S7_classification_rules"] = pd.read_csv(
            ext / "extB_incremental_value.tsv", sep="\t")
    if (ext / "extC_ppv_projection.tsv").exists():
        sheets["S8_ppv_projection"] = pd.read_csv(
            ext / "extC_ppv_projection.tsv", sep="\t")
    if (ext / "extA_false_resistant.tsv").exists():
        sheets["S9_carbnegative_resistance"] = pd.read_csv(
            ext / "extA_false_resistant.tsv", sep="\t")
    reg = []
    for drug in ("meropenem", "imipenem"):
        f = ext / f"sensitivity_regression_{drug}_clsi_bmd.tsv"
        if f.exists():
            reg.append(pd.read_csv(f, sep="\t"))
    if reg:
        sheets["S10_regression_estimates"] = pd.concat(reg, ignore_index=True)
    loo = []
    for drug in ("meropenem", "imipenem"):
        f = ext / f"sensitivity_leaveoneout_{drug}_clsi_bmd.tsv"
        if f.exists():
            loo.append(pd.read_csv(f, sep="\t"))
    if loo:
        sheets["S11_leave_one_ST_out"] = pd.concat(loo, ignore_index=True)

    # ---------- write ----------
    for name, df in sheets.items():
        df.to_csv(SUPP / f"{name}.tsv", sep="\t", index=False)
        print(f"  {name:<28} {len(df):>7,} rows")

    try:
        with pd.ExcelWriter(SUPP / "Supplementary_Tables.xlsx", engine="openpyxl") as xl:
            for name, df in sheets.items():
                df.to_excel(xl, sheet_name=name[:31], index=False)
        print(f"\n[out] {SUPP}/Supplementary_Tables.xlsx")
    except Exception as e:
        print(f"\n(xlsx skipped: {e} — install openpyxl for Excel output)")
    print(f"[out] {SUPP}/*.tsv")


if __name__ == "__main__":
    main()
