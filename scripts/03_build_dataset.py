#!/usr/bin/env python
"""03_build_dataset.py — reconcile BV-BRC phenotypes with assemblies on disk.

Drops phenotype rows whose genome failed to download, normalizes MIC values,
and reports what you actually have to model with.

Run from project root:
    python scripts/03_build_dataset.py

Outputs:
    data/processed/phenotypes_clean.tsv   modelling-ready phenotype table
    data/processed/dataset_summary.txt    per-drug counts / MIC coverage
"""
import math
import re
import sys
from pathlib import Path

import pandas as pd

PHEN = Path("data/raw/phenotypes/kp_amr_anchor.tsv")
GEN = Path("data/raw/genomes")
OUT = Path("data/processed")
OUT.mkdir(parents=True, exist_ok=True)


def parse_mic(val: str):
    """Fallback parser for a sign embedded in the value string."""
    if val is None:
        return (None, "")
    s = str(val).strip()
    if not s:
        return (None, "")
    m = re.match(r"^\s*(>=|<=|>|<|=)?\s*([0-9]*\.?[0-9]+)\s*$", s)
    if not m:
        return (None, "")
    return (float(m.group(2)), m.group(1) or "")


def censor_flags(sign: str):
    """BV-BRC stores the operator in `measurement_sign`.
      <= or <  -> LEFT censored  (true MIC at or below the value)
      >= or >  -> RIGHT censored (true MIC at or above the value)
    Returns (left, right)."""
    s = (sign or "").strip()
    return (s in {"<=", "<"}, s in {">=", ">"})


def main():
    if not PHEN.exists():
        sys.exit(f"ERROR: {PHEN} not found. Run 01_download_bvbrc.sh first.")

    df = pd.read_csv(PHEN, sep="\t", dtype=str, low_memory=False).fillna("")
    cols = {c.lower().strip(): c for c in df.columns}
    gid = cols.get("genome_id"); ab = cols.get("antibiotic")
    rp = cols.get("resistant_phenotype"); mv = cols.get("measurement_value")
    if not gid or not ab:
        sys.exit("ERROR: genome_id / antibiotic columns missing.")

    n0 = len(df)
    have = {p.stem for p in GEN.glob("*.fna") if p.stat().st_size > 1000}
    df["_has_genome"] = df[gid].astype(str).isin(have)
    kept = df[df["_has_genome"]].drop(columns="_has_genome").copy()

    # exclude non-MIC units (disk-diffusion zone diameters are mm, inverted scale)
    unit = cols.get("measurement_unit")
    if unit:
        bad = kept[unit].str.strip().str.lower().isin({"mm"})
        if bad.any():
            print(f"  dropping {int(bad.sum())} disk-diffusion (mm) records")
            kept = kept[~bad].copy()

    # normalize MIC
    if mv:
        sign_col = cols.get("measurement_sign")
        parsed = kept[mv].apply(parse_mic)
        kept["mic_value"] = [p[0] for p in parsed]
        if sign_col:
            signs = kept[sign_col].fillna("").astype(str)
        else:
            signs = pd.Series([p[1] for p in parsed], index=kept.index)
        flags = signs.apply(censor_flags)
        kept["left_censored"] = [f[0] for f in flags]
        kept["right_censored"] = [f[1] for f in flags]
        kept["mic_censor"] = signs.str.strip()
        kept["log2_mic"] = kept["mic_value"].apply(
            lambda v: None if v is None or pd.isna(v) or v <= 0 else round(math.log2(v), 4)
        )

    kept.to_csv(OUT / "phenotypes_clean.tsv", sep="\t", index=False)

    # ---- summary ----
    lines = []
    A = lines.append
    A("DATASET SUMMARY")
    A("=" * 60)
    A(f"assemblies on disk (>1KB) : {len(have)}")
    A(f"phenotype rows (raw)      : {n0}")
    A(f"phenotype rows (w/ genome): {len(kept)}   [dropped {n0 - len(kept)}]")
    A(f"unique genomes w/ pheno   : {kept[gid].nunique()}")
    A("")
    A("PER-DRUG BREAKDOWN")
    A("-" * 60)
    A(f"{'antibiotic':<28}{'rows':>7}{'genomes':>9}{'with MIC':>10}")
    for drug, g in kept.groupby(kept[ab].str.lower().str.strip()):
        nmic = g["mic_value"].notna().sum() if "mic_value" in g else 0
        A(f"{drug:<28}{len(g):>7}{g[gid].nunique():>9}{nmic:>10}")
    A("")
    if "left_censored" in kept:
        A("CENSORING")
        A("-" * 60)
        A(f"  left-censored  (<=, <) : {int(kept['left_censored'].sum())}")
        A(f"  right-censored (>=, >) : {int(kept['right_censored'].sum())}")
        A(f"  exact                  : {int((~kept['left_censored'] & ~kept['right_censored']).sum())}")
        A("  NOTE: >50% censored means plain medians/OLS are biased —")
        A("        use interval regression (scripts/07).")
        A("")
    if rp:
        A("RESISTANCE PHENOTYPE DISTRIBUTION")
        A("-" * 60)
        for k, v in kept[rp].str.strip().value_counts().items():
            if k:
                A(f"  {k:<26}{v:>7}")
    txt = "\n".join(lines)
    (OUT / "dataset_summary.txt").write_text(txt)
    print(txt)
    print(f"\nWrote {OUT/'phenotypes_clean.tsv'} and {OUT/'dataset_summary.txt'}")


if __name__ == "__main__":
    main()
