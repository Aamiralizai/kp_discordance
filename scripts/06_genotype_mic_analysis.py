#!/usr/bin/env python
"""06_genotype_mic_analysis.py — the first real scientific result.

Joins Kleborate genotypes to BV-BRC meropenem MICs and tests the mechanistic
hypothesis the toy model predicted:

    carbapenemase alone      -> moderate MIC increase
    porin loss alone         -> little/no increase
    BOTH together            -> SUPER-ADDITIVE increase

Works on partial Kleborate output (safe to run mid-run).

Usage:
    python scripts/06_genotype_mic_analysis.py
    python scripts/06_genotype_mic_analysis.py --drug imipenem
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

KLEB = Path("data/processed/kleborate_klebsiella_pneumo_complex_output.tsv")
KLEB_ALT = Path("data/processed/kleborate_results.tsv")
PHEN = Path("data/processed/phenotypes_clean.tsv")
OUTDIR = Path("results/tables"); FIGDIR = Path("results/figures")

C_STRAIN = "strain"
C_SPECIES = "enterobacterales__species__species"
C_QC = "general__contig_stats__QC_warnings"
C_ST = "klebsiella_pneumo_complex__mlst__ST"
C_CARB = "klebsiella_pneumo_complex__amr__Bla_Carb_acquired"
C_OMP = "klebsiella_pneumo_complex__amr__Omp_mutations"
C_ESBL = "klebsiella_pneumo_complex__amr__Bla_ESBL_acquired"


def has_carb(v):
    """Any acquired carbapenemase present?"""
    if not isinstance(v, str):
        return False
    s = v.strip()
    return s not in {"", "-", "nan", "NA", "None"}


def carb_family(v):
    """Coarse carbapenemase family: KPC / NDM / OXA-48-like / VIM / IMP / other."""
    if not has_carb(v):
        return "none"
    s = v.upper()
    fams = []
    if "KPC" in s: fams.append("KPC")
    if "NDM" in s: fams.append("NDM")
    if "OXA-48" in s or "OXA-181" in s or "OXA-232" in s or "OXA-244" in s:
        fams.append("OXA-48-like")
    if "VIM" in s: fams.append("VIM")
    if "IMP" in s: fams.append("IMP")
    if not fams:
        return "other"
    return "+".join(sorted(set(fams)))


def porin_status(v):
    """Classify OmpK35/36 disruption from Kleborate 3.x HGVS-style Omp_mutations.

    Observed encodings:
      OmpK36:p.134_135insGlyAsp   -> GD loop-3 insertion (narrowed pore)
      OmpK36:p.136_137insThrAsp   -> TD loop-3 insertion
      OmpK35:p.Glu42fs            -> frameshift, truncation
      OmpK36:p.Ala183fs           -> frameshift, OmpK36 loss
      OmpK36:c.25C>T              -> known premature stop, OmpK36 loss
      OmpK35:c.G175del            -> deletion, truncation

    Returns a GRADED class, ordered by expected permeability impact:
      none < ompK35_loss < ompK36_loop3 < ompK36_loss
    (OmpK36 disruption dominates when both are present.)
    """
    if not isinstance(v, str) or v.strip() in {"", "-", "nan", "NA", "None"}:
        return "none"
    k36_loss = k36_loop3 = k35_loss = False
    for tok in re.split(r"[;,]", v):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"(OmpK3[56])\s*:\s*(.+)", tok, flags=re.I)
        if not m:
            continue
        gene, chg = m.group(1).upper(), m.group(2)
        truncating = bool(re.search(r"(fs|del|dup|\*|Ter)", chg, flags=re.I)) \
            or bool(re.search(r"c\.25C>T", chg, flags=re.I))   # known premature stop
        insertion = bool(re.search(r"\d+_\d+ins", chg, flags=re.I))
        if gene == "OMPK36":
            if truncating:
                k36_loss = True
            elif insertion:
                k36_loop3 = True
        elif gene == "OMPK35":
            if truncating:
                k35_loss = True
    if k36_loss:
        return "ompK36_loss"
    if k36_loop3:
        return "ompK36_loop3"
    if k35_loss:
        return "ompK35_loss"
    return "none"


PORIN_ORDER = ["none", "ompK35_loss", "ompK36_loop3", "ompK36_loss"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="meropenem")
    ap.add_argument("--min-group", type=int, default=5)
    a = ap.parse_args()

    kleb_path = KLEB if KLEB.exists() else KLEB_ALT
    if not kleb_path.exists():
        sys.exit(f"ERROR: no Kleborate table. Looked for {KLEB} and {KLEB_ALT}.\n"
                 "Run scripts/05_merge_kleborate.py first.")
    if not PHEN.exists():
        sys.exit(f"ERROR: {PHEN} not found. Run scripts/03_build_dataset.py first.")

    OUTDIR.mkdir(parents=True, exist_ok=True); FIGDIR.mkdir(parents=True, exist_ok=True)

    kb = pd.read_csv(kleb_path, sep="\t", dtype=str, low_memory=False)
    if C_CARB not in kb.columns:
        sys.exit(f"ERROR: {kleb_path} lacks '{C_CARB}'.\n"
                 "You may have merged the hAMRonization table by mistake — "
                 "use the fixed 05_merge_kleborate.py.")
    print(f"Kleborate genotypes : {len(kb)} genomes  ({kleb_path.name})")

    # ---- QC filter ----
    n0 = len(kb)
    if C_SPECIES in kb.columns:
        kb = kb[kb[C_SPECIES].astype(str).str.contains("pneumoniae", case=False, na=False)]
    print(f"  after species filter (K. pneumoniae): {len(kb)}  [dropped {n0 - len(kb)}]")

    # ---- genotype features ----
    kb["carbapenemase"] = kb[C_CARB].apply(has_carb)
    kb["carb_family"] = kb[C_CARB].apply(carb_family)
    kb["porin"] = kb[C_OMP].apply(porin_status) if C_OMP in kb.columns else "none"
    kb["porin_loss"] = kb["porin"].isin(["ompK36_loop3", "ompK36_loss"])
    kb["any_omp"] = kb["porin"] != "none"
    kb["ST"] = kb[C_ST] if C_ST in kb.columns else "NA"
    kb["genome_id"] = kb[C_STRAIN].astype(str).str.strip()

    # ---- phenotypes ----
    ph = pd.read_csv(PHEN, sep="\t", dtype=str, low_memory=False)
    abcol = next(c for c in ph.columns if c.lower().strip() == "antibiotic")
    gidcol = next(c for c in ph.columns if c.lower().strip() == "genome_id")
    ph = ph[ph[abcol].str.lower().str.strip() == a.drug.lower()].copy()
    ph["genome_id"] = ph[gidcol].astype(str).str.strip()
    ph["log2_mic"] = pd.to_numeric(ph.get("log2_mic"), errors="coerce")
    ph["mic_censor"] = ph.get("mic_censor", "").fillna("")
    ph = ph.dropna(subset=["log2_mic"])
    # one record per genome (median if duplicated)
    ph = ph.groupby("genome_id", as_index=False).agg(
        log2_mic=("log2_mic", "median"), n_records=("log2_mic", "size"))
    print(f"{a.drug} MICs        : {len(ph)} genomes")

    df = kb.merge(ph, on="genome_id", how="inner")
    print(f"JOINED              : {len(df)} genomes with genotype + {a.drug} MIC\n")
    if len(df) < 20:
        print("!! Very few joined genomes — Kleborate is probably still running.")
        print("!! Re-run when typing completes.\n")

    # ---- the 2x2 hypothesis test ----
    def grp(r):
        c, p = r["carbapenemase"], r["porin_loss"]
        return ("both" if c and p else "carbapenemase_only" if c
                else "porin_only" if p else "neither")
    df["group"] = df.apply(grp, axis=1)

    order = ["neither", "porin_only", "carbapenemase_only", "both"]
    print("=" * 74)
    print(f"{a.drug.upper()} MIC BY GENOTYPE GROUP")
    print("=" * 74)
    print(f"{'group':<22}{'n':>6}{'median log2':>13}{'median mg/L':>13}{'IQR log2':>18}")
    stats_rows = []
    for g in order:
        s = df.loc[df["group"] == g, "log2_mic"]
        if len(s) == 0:
            print(f"{g:<22}{0:>6}{'-':>13}{'-':>13}{'-':>18}")
            continue
        med = s.median(); q1, q3 = s.quantile(.25), s.quantile(.75)
        print(f"{g:<22}{len(s):>6}{med:>13.2f}{2**med:>13.1f}{f'{q1:.2f}-{q3:.2f}':>18}")
        stats_rows.append(dict(group=g, n=len(s), median_log2=med, median_mgL=2**med,
                               q1=q1, q3=q3))
    res = pd.DataFrame(stats_rows)

    # ---- super-additivity ----
    def med(g):
        s = df.loc[df["group"] == g, "log2_mic"]
        return (s.median(), len(s))
    m_n, n_n = med("neither"); m_p, n_p = med("porin_only")
    m_c, n_c = med("carbapenemase_only"); m_b, n_b = med("both")

    print("\n" + "=" * 74)
    print("SUPER-ADDITIVITY TEST  (the mechanistic prediction)")
    print("=" * 74)
    if min(n_n, n_p, n_c, n_b) < a.min_group:
        print(f"Insufficient data — every group needs >= {a.min_group} genomes.")
        print(f"  n: neither={n_n}, porin={n_p}, carb={n_c}, both={n_b}")
    else:
        e_p, e_c, e_b = m_p - m_n, m_c - m_n, m_b - m_n
        additive = e_p + e_c
        print(f"  porin effect alone        : {e_p:+.2f} log2 ({2**e_p:.1f}x)")
        print(f"  carbapenemase alone       : {e_c:+.2f} log2 ({2**e_c:.1f}x)")
        print(f"  additive expectation      : {additive:+.2f} log2 ({2**additive:.1f}x)")
        print(f"  OBSERVED both             : {e_b:+.2f} log2 ({2**e_b:.1f}x)")
        excess = e_b - additive
        print(f"  EXCESS over additive      : {excess:+.2f} log2 ({2**excess:.1f}x)")
        u = stats.mannwhitneyu(df.loc[df.group == "both", "log2_mic"],
                               df.loc[df.group == "carbapenemase_only", "log2_mic"],
                               alternative="greater")
        print(f"\n  both > carbapenemase_only : Mann-Whitney U p = {u.pvalue:.3g}")
        if excess > 0.5 and u.pvalue < 0.05:
            print("\n  ==> SUPER-ADDITIVE synergy detected. This is the mechanistic")
            print("      prediction confirmed on real genomes.")
        elif excess > 0.5:
            print("\n  ==> Super-additive in direction, not yet significant.")
        else:
            print("\n  ==> No super-additivity in this sample. Worth investigating:")
            print("      MIC censoring (>=64 caps the ceiling) can mask synergy.")

    # censoring caveat
    print("\nNOTE: censored MICs ('>=64') were treated as their numeric value.")
    print("      If many 'both' genomes are censored at the top of the scale,")
    print("      the true synergy is UNDERESTIMATED here.")

    # ---- carbapenemase family breakdown ----
    print("\n" + "=" * 74)
    print("BY CARBAPENEMASE FAMILY (carriers only)")
    print("=" * 74)
    sub = df[df["carbapenemase"]]
    if len(sub):
        fam = (sub.groupby("carb_family")["log2_mic"]
               .agg(n="size", median="median").sort_values("n", ascending=False))
        fam["median_mgL"] = 2 ** fam["median"]
        print(fam.round(2).to_string())

    # ---- porin class breakdown ----
    print("\n" + "=" * 74)
    print("BY PORIN DISRUPTION CLASS")
    print("=" * 74)
    pk = df.groupby("porin")["log2_mic"].agg(n="size", median="median")
    pk = pk.reindex([g for g in PORIN_ORDER if g in pk.index])
    pk["median_mgL"] = 2 ** pk["median"]
    print(pk.round(2).to_string())
    print("  (ordered by expected permeability impact)")

    # --- two-way: carbapenemase family x porin class (the key table) ---
    print("\n" + "=" * 74)
    print("TWO-WAY: carbapenemase family x porin class (median log2 MIC)")
    print("=" * 74)
    fams = ["none", "OXA-48-like", "KPC", "NDM"]
    sub2 = df[df["carb_family"].isin(fams)]
    if len(sub2):
        piv = sub2.pivot_table(index="carb_family", columns="porin",
                               values="log2_mic", aggfunc="median")
        cnt = sub2.pivot_table(index="carb_family", columns="porin",
                               values="log2_mic", aggfunc="size")
        piv = piv.reindex(index=[f for f in fams if f in piv.index],
                          columns=[c for c in PORIN_ORDER if c in piv.columns])
        cnt = cnt.reindex(index=piv.index, columns=piv.columns)
        print("median log2 MIC:"); print(piv.round(2).to_string())
        print("\nn per cell:"); print(cnt.fillna(0).astype(int).to_string())
        print("\n  Read across each row: does porin disruption raise MIC MORE")
        print("  in carbapenemase carriers than in the 'none' row? That gradient")
        print("  IS the enzyme x permeability interaction the model predicts.")

    # ---- save ----
    df_out = df[["genome_id", "ST", "carb_family", "porin", "carbapenemase",
                 "porin_loss", "group", "log2_mic"]]
    df_out.to_csv(OUTDIR / f"genotype_mic_{a.drug}.tsv", sep="\t", index=False)
    res.to_csv(OUTDIR / f"group_summary_{a.drug}.tsv", sep="\t", index=False)
    print(f"\n[out] {OUTDIR}/genotype_mic_{a.drug}.tsv")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        present = [g for g in order if (df["group"] == g).sum() > 0]
        data = [df.loc[df.group == g, "log2_mic"].values for g in present]
        fig, ax = plt.subplots(figsize=(8, 5))
        lbls = [g.replace("_", "\n") for g in present]
        try:
            bp = ax.boxplot(data, tick_labels=lbls, showfliers=False, patch_artist=True)
        except TypeError:      # matplotlib < 3.9
            bp = ax.boxplot(data, labels=lbls, showfliers=False, patch_artist=True)
        cols = {"neither": "#2b8cbe", "porin_only": "#d95f0e",
                "carbapenemase_only": "#238b45", "both": "#cb181d"}
        for patch, g in zip(bp["boxes"], present):
            patch.set_facecolor(cols.get(g, "#999")); patch.set_alpha(.65)
        for i, (g, d) in enumerate(zip(present, data), start=1):
            x = np.random.normal(i, 0.06, len(d))
            ax.plot(x, d, ".", color="k", alpha=.18, ms=3)
            ax.text(i, ax.get_ylim()[1], f"n={len(d)}", ha="center", va="bottom", fontsize=9)
        ax.set_ylabel(f"log2 {a.drug} MIC (mg/L)")
        ax.set_title(f"{a.drug.capitalize()} MIC by carbapenemase / porin genotype",
                     fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        fig.savefig(FIGDIR / f"genotype_mic_{a.drug}.png", dpi=150, bbox_inches="tight")
        print(f"[out] {FIGDIR}/genotype_mic_{a.drug}.png")
    except Exception as e:
        print(f"(figure skipped: {e})")


if __name__ == "__main__":
    main()
