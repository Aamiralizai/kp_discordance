#!/usr/bin/env python
"""08_discordance_analysis.py — quantify carbapenemase genotype-phenotype discordance.

THE QUESTION: what fraction of carbapenemase-POSITIVE K. pneumoniae are
phenotypically SUSCEPTIBLE to carbapenems, and what genomic features predict it?

WHY CENSORING DOESN'T BREAK THIS (unlike the interaction analysis):
Discordance is a BREAKPOINT question, not an effect-size question.
  '<=1'  -> unambiguously susceptible (true MIC is at or below 1)
  '>8'   -> unambiguously resistant
Censored records still resolve which side of the breakpoint an isolate sits on,
even though they cannot resolve exact MIC values. Records are excluded only when
the censoring bound straddles the breakpoint (e.g. '<=8' when breakpoint is 2 —
true value could be either side).

Breakpoints (CLSI/EUCAST, Enterobacterales):
  meropenem : S <= 1,  R >= 4   mg/L
  imipenem  : S <= 1,  R >= 4   mg/L

Usage:
    python scripts/08_discordance_analysis.py
    python scripts/08_discordance_analysis.py --drug imipenem
"""
import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

KLEB = Path("data/processed/kleborate_klebsiella_pneumo_complex_output.tsv")
PHEN = Path("data/processed/phenotypes_clean.tsv")
OUTDIR = Path("results/tables"); FIGDIR = Path("results/figures")

C_STRAIN = "strain"
C_SPECIES = "enterobacterales__species__species"
C_ST = "klebsiella_pneumo_complex__mlst__ST"
C_CARB = "klebsiella_pneumo_complex__amr__Bla_Carb_acquired"
C_OMP = "klebsiella_pneumo_complex__amr__Omp_mutations"
C_ESBL = "klebsiella_pneumo_complex__amr__Bla_ESBL_acquired"

# CLSI Enterobacterales carbapenem breakpoints (mg/L)
BREAKPOINTS = {"meropenem": dict(S=1.0, R=4.0),
               "imipenem": dict(S=1.0, R=4.0),
               "ertapenem": dict(S=0.5, R=2.0)}


def has_carb(v):
    return isinstance(v, str) and v.strip() not in {"", "-", "nan", "NA", "None"}


def carb_family(v):
    if not has_carb(v):
        return "none"
    s = v.upper(); f = []
    if "KPC" in s: f.append("KPC")
    if "NDM" in s: f.append("NDM")
    if any(x in s for x in ("OXA-48", "OXA-181", "OXA-232", "OXA-244", "OXA-162")):
        f.append("OXA-48-like")
    if "VIM" in s: f.append("VIM")
    if "IMP" in s: f.append("IMP")
    return "+".join(sorted(set(f))) if f else "other"


def carb_allele(v):
    """Specific allele string, cleaned of Kleborate confidence flags."""
    if not has_carb(v):
        return "none"
    toks = [re.sub(r"[\*\?\^\$]", "", t).strip() for t in re.split(r"[;,]", v)]
    toks = [t for t in toks if t]
    return ";".join(sorted(set(toks))) if toks else "none"


def porin_status(v):
    if not isinstance(v, str) or v.strip() in {"", "-", "nan", "NA", "None"}:
        return "none"
    k36_loss = k36_loop3 = k35_loss = False
    for tok in re.split(r"[;,]", v):
        m = re.match(r"(OmpK3[56])\s*:\s*(.+)", tok.strip(), flags=re.I)
        if not m:
            continue
        gene, chg = m.group(1).upper(), m.group(2)
        trunc = bool(re.search(r"(fs|del|dup|\*|Ter)", chg, flags=re.I)) or \
            bool(re.search(r"c\.25C>T", chg, flags=re.I))
        ins = bool(re.search(r"\d+_\d+ins", chg, flags=re.I))
        if gene == "OMPK36":
            if trunc: k36_loss = True
            elif ins: k36_loop3 = True
        elif gene == "OMPK35" and trunc:
            k35_loss = True
    if k36_loss: return "ompK36_loss"
    if k36_loop3: return "ompK36_loop3"
    if k35_loss: return "ompK35_loss"
    return "none"


def classify(mic, lc, rc, bp_s, bp_r):
    """Resolve S/I/R accounting for censoring. Returns S / R / I / ambiguous.

    lc = left-censored  ('<=x' -> true value <= x)
    rc = right-censored ('>x'  -> true value >  x)
    """
    if mic is None or (isinstance(mic, float) and math.isnan(mic)):
        return "ambiguous"
    if lc:
        # true value <= mic. Susceptible only if the BOUND is already <= S.
        return "S" if mic <= bp_s else "ambiguous"
    if rc:
        # true value > mic. Resistant only if the BOUND is already >= R.
        return "R" if mic >= bp_r else "ambiguous"
    if mic <= bp_s:
        return "S"
    if mic >= bp_r:
        return "R"
    return "I"


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    z = 1.959964
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (p, max(0, c - h), min(1, c + h))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="meropenem")
    ap.add_argument("--min-n", type=int, default=10)
    ap.add_argument("--standard", default=None,
                    help="restrict to a testing_standard, e.g. CLSI")
    ap.add_argument("--method", default=None,
                    help="restrict to a laboratory_typing_method, e.g. 'Broth dilution'")
    ap.add_argument("--dedup", action="store_true",
                    help="keep one isolate per ST x carb_allele x porin (clonality control)")
    ap.add_argument("--tag", default="",
                    help="suffix for output filenames")
    a = ap.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True); FIGDIR.mkdir(parents=True, exist_ok=True)
    bp = BREAKPOINTS.get(a.drug.lower())
    if bp is None:
        sys.exit(f"No breakpoints defined for {a.drug}")
    if not KLEB.exists(): sys.exit(f"ERROR: {KLEB} not found")
    if not PHEN.exists(): sys.exit(f"ERROR: {PHEN} not found")

    kb = pd.read_csv(KLEB, sep="\t", dtype=str, low_memory=False)
    kb = kb[kb[C_SPECIES].astype(str).str.contains("pneumoniae", case=False, na=False)]
    kb["genome_id"] = kb[C_STRAIN].astype(str).str.strip()
    kb["carb_family"] = kb[C_CARB].apply(carb_family)
    kb["carb_allele"] = kb[C_CARB].apply(carb_allele)
    kb["carb_pos"] = kb[C_CARB].apply(has_carb)
    kb["porin"] = kb[C_OMP].apply(porin_status)
    kb["ST"] = kb[C_ST].astype(str)
    kb["esbl"] = kb[C_ESBL].apply(has_carb) if C_ESBL in kb.columns else False

    ph = pd.read_csv(PHEN, sep="\t", dtype=str, low_memory=False)
    abc = next(c for c in ph.columns if c.lower().strip() == "antibiotic")
    gid = next(c for c in ph.columns if c.lower().strip() == "genome_id")
    ph = ph[ph[abc].str.lower().str.strip() == a.drug.lower()].copy()
    ph["genome_id"] = ph[gid].astype(str).str.strip()
    ph["mic"] = pd.to_numeric(ph.get("mic_value"), errors="coerce")
    if "left_censored" in ph.columns:
        ph["lc"] = ph["left_censored"].astype(str).str.lower().isin({"true", "1"})
        ph["rc"] = ph["right_censored"].astype(str).str.lower().isin({"true", "1"})
    else:
        s = ph.get("mic_censor", pd.Series("", index=ph.index)).fillna("").astype(str)
        ph["lc"] = s.str.strip().isin(["<=", "<"]); ph["rc"] = s.str.strip().isin([">=", ">"])
    ph = ph.dropna(subset=["mic"])
    ph = ph.groupby("genome_id", as_index=False).agg(
        mic=("mic", "median"), lc=("lc", "max"), rc=("rc", "max"))

    df = kb.merge(ph, on="genome_id", how="inner")

    # ---- methodological stratification (removes the dominant confounder) ----
    RAW = Path("data/raw/phenotypes/kp_amr_anchor.tsv")
    if (a.standard or a.method) and RAW.exists():
        raw = pd.read_csv(RAW, sep="\t", dtype=str, low_memory=False).fillna("")
        rabc = next((c for c in raw.columns if c.lower().strip() == "antibiotic"), None)
        rgid = next((c for c in raw.columns if c.lower().strip() == "genome_id"), None)
        raw = raw[raw[rabc].str.lower().str.strip() == a.drug.lower()]
        keep = [rgid] + [c for c in raw.columns
                         if c.lower() in ("testing_standard", "laboratory_typing_method")]
        raw = raw[keep].drop_duplicates(subset=[rgid]).rename(columns={rgid: "genome_id"})
        n_before = len(df)
        df = df.merge(raw, on="genome_id", how="left")
        if a.standard and "testing_standard" in df.columns:
            df = df[df["testing_standard"].astype(str).str.strip().str.upper()
                    == a.standard.upper()]
        if a.method and "laboratory_typing_method" in df.columns:
            df = df[df["laboratory_typing_method"].astype(str).str.strip().str.lower()
                    == a.method.strip().lower()]
        print(f"[stratified] standard={a.standard} method={a.method}: "
              f"{n_before} -> {len(df)} isolates")

    df["call"] = [classify(m, l, r, bp["S"], bp["R"])
                  for m, l, r in zip(df["mic"], df["lc"], df["rc"])]

    if a.dedup:
        n_before = len(df)
        df = df.drop_duplicates(subset=["ST", "carb_allele", "porin"])
        print(f"[dedup] one per ST x allele x porin: {n_before} -> {len(df)} isolates")

    print("=" * 76)
    print(f"CARBAPENEMASE GENOTYPE-PHENOTYPE DISCORDANCE — {a.drug.upper()}")
    print(f"breakpoints: S <= {bp['S']} mg/L, R >= {bp['R']} mg/L")
    print("=" * 76)
    print(f"genomes with genotype + phenotype : {len(df)}")
    print(f"  resolvable despite censoring    : {(df['call'] != 'ambiguous').sum()}"
          f"  ({(df['call'] != 'ambiguous').mean():.1%})")
    print(f"  ambiguous (censor bound straddles breakpoint): "
          f"{(df['call'] == 'ambiguous').sum()}")
    print("\n  NOTE: censored records ARE usable here — '<=1' resolves as susceptible,")
    print("  '>8' as resistant. Only bounds straddling the breakpoint are dropped.")

    d = df[df["call"] != "ambiguous"].copy()
    d["carb_pos"] = d["carb_pos"].astype(bool)

    # ---------- headline: discordance rates ----------
    print("\n" + "=" * 76)
    print("HEADLINE: discordance in both directions")
    print("=" * 76)
    pos = d[d["carb_pos"]]
    neg = d[~d["carb_pos"]]
    fs = (pos["call"] == "S").sum()          # false-susceptible (gene+, pheno S)
    fr = (neg["call"] == "R").sum()          # false-resistant   (gene-, pheno R)
    p1, lo1, hi1 = wilson(fs, len(pos))
    p2, lo2, hi2 = wilson(fr, len(neg))
    print(f"carbapenemase-POSITIVE : n = {len(pos)}")
    print(f"   phenotypically SUSCEPTIBLE : {fs}  ({p1:.1%}, 95% CI {lo1:.1%}-{hi1:.1%})")
    print(f"   >>> genomic surveillance would OVER-call resistance in these")
    print(f"carbapenemase-NEGATIVE : n = {len(neg)}")
    print(f"   phenotypically RESISTANT   : {fr}  ({p2:.1%}, 95% CI {lo2:.1%}-{hi2:.1%})")
    print(f"   >>> genomic surveillance would MISS these")

    # ---------- by enzyme family ----------
    print("\n" + "=" * 76)
    print("DISCORDANCE BY CARBAPENEMASE FAMILY (the mechanistic prediction)")
    print("  expectation: OXA-48-like HIGH (weak meropenem hydrolysis),")
    print("               NDM/KPC LOW (strong hydrolysis)")
    print("=" * 76)
    rows = []
    for fam, g in pos.groupby("carb_family"):
        if len(g) < a.min_n:
            continue
        k = (g["call"] == "S").sum()
        p, lo, hi = wilson(k, len(g))
        rows.append(dict(family=fam, n=len(g), n_susceptible=k,
                         pct_susceptible=100 * p, ci_lo=100 * lo, ci_hi=100 * hi))
    fam_tab = pd.DataFrame(rows).sort_values("pct_susceptible", ascending=False)
    print(fam_tab.round(1).to_string(index=False))

    # ---------- the mechanism: porin state within enzyme ----------
    print("\n" + "=" * 76)
    print("MECHANISM: does porin state explain discordance within each family?")
    print("=" * 76)
    rows = []
    for fam in fam_tab["family"]:
        sub = pos[pos["carb_family"] == fam]
        for pr, g in sub.groupby("porin"):
            if len(g) < a.min_n:
                continue
            k = (g["call"] == "S").sum()
            p, lo, hi = wilson(k, len(g))
            rows.append(dict(family=fam, porin=pr, n=len(g), n_susc=k,
                             pct_susceptible=100 * p, ci_lo=100 * lo, ci_hi=100 * hi))
    mech = pd.DataFrame(rows)
    if len(mech):
        print(mech.round(1).to_string(index=False))
        # formal test within each family
        print("\n  Fisher tests (intact porin vs any porin disruption), within family:")
        for fam in mech["family"].unique():
            sub = pos[pos["carb_family"] == fam].copy()
            sub["disrupted"] = sub["porin"] != "none"
            t = pd.crosstab(sub["disrupted"], sub["call"] == "S")
            if t.shape == (2, 2) and t.values.min() >= 0:
                try:
                    orr, pv = stats.fisher_exact(t.values)
                    print(f"    {fam:<14} OR = {orr:6.2f}   p = {pv:.3g}"
                          f"   (n = {len(sub)})")
                except Exception:
                    pass

    # ---------- allele resolution ----------
    print("\n" + "=" * 76)
    print("ALLELE-LEVEL (single residues change kinetics ~27-fold in OXA-48-like)")
    print("=" * 76)
    rows = []
    for al, g in pos.groupby("carb_allele"):
        if len(g) < a.min_n:
            continue
        k = (g["call"] == "S").sum()
        p, lo, hi = wilson(k, len(g))
        rows.append(dict(allele=al[:34], n=len(g), pct_susceptible=100 * p,
                         ci_lo=100 * lo, ci_hi=100 * hi))
    al_tab = pd.DataFrame(rows).sort_values("pct_susceptible", ascending=False)
    print(al_tab.round(1).head(20).to_string(index=False))

    # ---------- clinical framing ----------
    print("\n" + "=" * 76)
    print("SURVEILLANCE IMPACT")
    print("=" * 76)
    tp = (pos["call"] == "R").sum(); fp = fs
    fn = (neg["call"] == "R").sum(); tn = (neg["call"] == "S").sum()
    ppv = tp / (tp + fp) if (tp + fp) else np.nan
    npv = tn / (tn + fn) if (tn + fn) else np.nan
    sens = tp / (tp + fn) if (tp + fn) else np.nan
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    print("Treating 'carbapenemase gene present' as a resistance test:")
    print(f"   PPV (gene+ really resistant)   : {ppv:.1%}")
    print(f"   NPV (gene- really susceptible) : {npv:.1%}")
    print(f"   sensitivity                    : {sens:.1%}")
    print(f"   specificity                    : {spec:.1%}")
    print("\n  A PPV well below 100% means purely genomic surveillance")
    print("  systematically OVER-calls carbapenem resistance.")

    out = d[["genome_id", "ST", "carb_family", "carb_allele", "carb_pos",
             "porin", "mic", "lc", "rc", "call"]]
    out.to_csv(OUTDIR / f"discordance_{a.drug}{a.tag}.tsv", sep="\t", index=False)
    fam_tab.to_csv(OUTDIR / f"discordance_by_family_{a.drug}{a.tag}.tsv", sep="\t", index=False)
    if len(mech):
        mech.to_csv(OUTDIR / f"discordance_mechanism_{a.drug}{a.tag}.tsv", sep="\t", index=False)
    print(f"\n[out] {OUTDIR}/discordance_{a.drug}.tsv")


if __name__ == "__main__":
    main()
