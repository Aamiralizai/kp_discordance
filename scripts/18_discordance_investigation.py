#!/usr/bin/env python3
"""
18_discordance_investigation.py — what the discordant isolates actually are.

The editor's second route to a non-incremental contribution was to reveal
biology through the discrepancies themselves rather than only to quantify them.
This script interrogates the two discordant cells directly.

    resistant despite no detected carbapenemase
    susceptible despite a detected carbapenemase

FOUR ANALYSES

1. LINEAGE ENRICHMENT among resistant carbapenemase-negative isolates.
   Per sequence type, a Fisher test of resistance against all other
   carbapenemase-negative isolates, Benjamini-Hochberg corrected. A clone that
   is resistant far above background without a carbapenemase is either carrying
   a mechanism the typer does not report, or carrying a carbapenemase the
   assembly did not recover. Both are findable, and the two are distinguished
   by analysis 4.

2. BREAKPOINT PROXIMITY among susceptible carbapenemase-positive isolates.
   An isolate whose MIC sits at the top of the susceptible range is not
   evidence of a silent carbapenemase in the same sense as one several
   dilutions below. Reporting the distribution rather than the count is the
   difference between a finding and an artefact of where the breakpoint falls.

3. ALLELE-LEVEL STRUCTURE within family. Family-level aggregation can hide a
   single weakly hydrolysing allele driving an apparent family effect.

4. ASSEMBLY QUALITY as a competing explanation, when QC metrics are available.
   This is the analysis that connects the two editorial criticisms. A
   fragmented assembly can fail to recover a carbapenemase gene, which
   produces a resistant carbapenemase-negative isolate that is an assembly
   artefact. The same fragmentation can truncate ompK36 and be read as porin
   disruption. If discordance is associated with assembly quality, part of the
   reported effect is technical, and the rerun will change it. If it is not,
   the finding survives, and the paper can say so with evidence rather than
   assertion.

USAGE

    python scripts/18_discordance_investigation.py
    python scripts/18_discordance_investigation.py \
        --qc data/interim/assembly_qc_metrics.tsv
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

TAB = Path("results/tables")
SUP = Path("results/supplementary")

OMPK36_DISRUPTED = {"ompK36_loop3", "ompK36_loss", "ompK36_reduced_expression"}


def wilson(k: int, n: int, z: float = 1.959963984540054):
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = (z / d) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (round(100 * p, 1), round(100 * max(0.0, centre - half), 1),
            round(100 * min(1.0, centre + half), 1))


def bh(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg, returning q in the input order."""
    n = len(pvals)
    if n == 0:
        return []
    order = np.argsort(pvals)
    q = np.empty(n, dtype=float)
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = n - rank + 1
        prev = min(prev, pvals[idx] * n / i)
        q[idx] = prev
    return [float(x) for x in q]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=SUP / "S1_genome_level_data.tsv")
    ap.add_argument("--drug", default="meropenem", choices=["meropenem", "imipenem"])
    ap.add_argument("--qc", type=Path,
                    help="assembly_qc_metrics.tsv from 00b_qc_assembly.py")
    ap.add_argument("--min-st", type=int, default=10,
                    help="minimum isolates per ST to test for enrichment")
    ap.add_argument("--tag", default="")
    ap.add_argument("--outdir", type=Path, default=TAB)
    a = ap.parse_args()

    if not a.data.exists():
        sys.exit(f"Missing {a.data}")
    d = pd.read_csv(a.data, sep="\t", dtype=str, low_memory=False)
    col = f"{a.drug}_interpretation"
    d = d[d[col].isin(["S", "R"])].copy()
    d["y"] = (d[col] == "R").astype(int)
    d["cp"] = d["carbapenemase_present"].astype(str).str.lower().isin(
        {"true", "1", "yes"})
    d["sequence_type"] = d["sequence_type"].fillna("unknown")
    d["porin_class"] = d["porin_class"].fillna("unknown")
    d["ompk36"] = d["porin_class"].isin(OMPK36_DISRUPTED)

    a.outdir.mkdir(parents=True, exist_ok=True)
    neg, pos = d[~d.cp], d[d.cp]
    fr = neg[neg.y == 1]      # resistant, carbapenemase-negative
    fs = pos[pos.y == 0]      # susceptible, carbapenemase-positive

    print(f"drug : {a.drug}")
    print(f"binary set {len(d):,} | carbapenemase-negative {len(neg):,} "
          f"(resistant {len(fr):,}) | carbapenemase-positive {len(pos):,} "
          f"(susceptible {len(fs):,})")

    # ---------------------------------------------- 1. lineage enrichment
    print("\n" + "=" * 78)
    print("1. LINEAGE ENRICHMENT  —  resistance without a detected carbapenemase")
    print("=" * 78)
    rows = []
    for st, g in neg.groupby("sequence_type"):
        if len(g) < a.min_st:
            continue
        rest = neg[neg.sequence_type != st]
        table = [[int(g.y.sum()), len(g) - int(g.y.sum())],
                 [int(rest.y.sum()), len(rest) - int(rest.y.sum())]]
        orr, p = stats.fisher_exact(table)
        pct, lo, hi = wilson(int(g.y.sum()), len(g))
        rows.append(dict(sequence_type=st, n=len(g), n_resistant=int(g.y.sum()),
                         resistance_pct=pct, ci_lower=lo, ci_upper=hi,
                         odds_ratio=(round(orr, 3) if np.isfinite(orr) else np.inf),
                         p_value=p,
                         ompK36_disrupted_pct=round(100 * g.ompk36.mean(), 1)))
    enr = pd.DataFrame(rows)
    if not enr.empty:
        enr["q_value"] = bh(enr.p_value.tolist())
        enr = enr.sort_values("resistance_pct", ascending=False)
        bg = 100 * neg.y.mean()
        print(f"background resistance among carbapenemase-negative isolates: "
              f"{bg:.1f}%  (n = {len(neg):,})")
        print(f"\n{'ST':14s}{'n':>5s}{'R':>5s}{'R%':>7s}  {'95% CI':>13s}"
              f"{'OR':>8s}{'q':>9s}   ompK36-disrupted%")
        for _, r in enr.iterrows():
            star = " *" if r.q_value < 0.05 else "  "
            print(f"{r.sequence_type:14s}{r.n:5d}{r.n_resistant:5d}"
                  f"{r.resistance_pct:7.1f}  [{r.ci_lower:5.1f}-{r.ci_upper:5.1f}]"
                  f"{r.odds_ratio:8.2f}{r.q_value:9.2e}{star}"
                  f"   {r.ompK36_disrupted_pct:5.1f}")
        p = a.outdir / f"discordance_lineage_enrichment{a.tag}.tsv"
        enr.to_csv(p, sep="\t", index=False)
        print(f"\n-> {p}")
        sig = enr[(enr.q_value < 0.05) & (enr.odds_ratio > 1)]
        if not sig.empty:
            print(f"\n{len(sig)} lineage(s) resistant above background at q < 0.05 "
                  f"without a detected carbapenemase.")
            print("These are the candidates for an unreported mechanism or an "
                  "unrecovered gene.\nAnalysis 4 distinguishes the two.")

    # ------------------------------------------- 2. breakpoint proximity
    print("\n" + "=" * 78)
    print("2. BREAKPOINT PROXIMITY  —  susceptibility despite a carbapenemase")
    print("=" * 78)
    up = pd.to_numeric(fs[f"{a.drug}_MIC_upper_mgL"], errors="coerce")
    prox = (up.value_counts().sort_index().rename_axis("MIC_upper_mgL")
            .reset_index(name="n_isolates"))
    if not prox.empty:
        prox["pct_of_susceptible_gene_positive"] = (
            100 * prox.n_isolates / prox.n_isolates.sum()).round(1)
        prox["cumulative_pct"] = prox.pct_of_susceptible_gene_positive.cumsum().round(1)
        print(f"n with a numeric upper bound: {int(up.notna().sum())} of {len(fs)}")
        print(f"\n{'MIC upper (mg/L)':>18s}{'n':>6s}{'%':>8s}{'cum %':>8s}")
        for _, r in prox.iterrows():
            print(f"{r.MIC_upper_mgL:18.2f}{int(r.n_isolates):6d}"
                  f"{r.pct_of_susceptible_gene_positive:8.1f}"
                  f"{r.cumulative_pct:8.1f}")
        top = up.max()
        at_top = int((up == top).sum())
        print(f"\n{at_top} of {int(up.notna().sum())} "
              f"({100 * at_top / max(1, int(up.notna().sum())):.1f}%) sit at "
              f"{top:g} mg/L, the highest MIC still called susceptible.")
        print("A single doubling dilution would move these out of the "
              "susceptible category,\nso they are better described as "
              "borderline than as robustly susceptible.")
        p = a.outdir / f"discordance_breakpoint_proximity{a.tag}.tsv"
        prox.to_csv(p, sep="\t", index=False)
        print(f"-> {p}")

    # ------------------------------------------------ 3. allele structure
    print("\n" + "=" * 78)
    print("3. ALLELE-LEVEL STRUCTURE  —  susceptibility rate by allele")
    print("=" * 78)
    pos2 = pos.copy()
    pos2["allele"] = pos2["carbapenemase_allele"].fillna("unknown")
    rows = []
    for (fam, al), g in pos2.groupby(["carbapenemase_family", "allele"]):
        ns = int((g.y == 0).sum())
        pct, lo, hi = wilson(ns, len(g))
        rows.append(dict(family=fam, allele=al, n=len(g), n_susceptible=ns,
                         susceptible_pct=pct, ci_lower=lo, ci_upper=hi,
                         low_support="yes" if len(g) < a.min_st else "no"))
    al = pd.DataFrame(rows).sort_values(["family", "n"],
                                        ascending=[True, False])
    shown = al[al.low_support == "no"]
    print(f"{'family':18s}{'allele':14s}{'n':>5s}{'S':>5s}{'S%':>7s}  95% CI")
    for _, r in shown.iterrows():
        print(f"{r.family:18s}{r.allele:14s}{r.n:5d}{r.n_susceptible:5d}"
              f"{r.susceptible_pct:7.1f}  [{r.ci_lower:.1f}-{r.ci_upper:.1f}]")
    p = a.outdir / f"discordance_by_allele{a.tag}.tsv"
    al.to_csv(p, sep="\t", index=False)
    print(f"\n-> {p}  ({len(al) - len(shown)} alleles with n < {a.min_st} "
          f"written but not shown)")

    # ------------------------------------------------- 4. assembly quality
    print("\n" + "=" * 78)
    print("4. ASSEMBLY QUALITY AS A COMPETING EXPLANATION")
    print("=" * 78)
    if not a.qc:
        print("  Not run: no --qc metrics supplied.")
        print("  Produce them with scripts/00b_qc_assembly.py, then re-run this")
        print("  script with --qc. Until then the technical explanation for")
        print("  discordance is untested, and the manuscript should not claim")
        print("  otherwise.")
    else:
        if not a.qc.exists():
            sys.exit(f"Missing {a.qc}")
        q = pd.read_csv(a.qc, sep="\t", dtype=str, low_memory=False)
        key = "BVBRC_genome_id"
        q["sample"] = q["sample"].astype(str)
        d["_k"] = d[key].astype(str)
        m = d.merge(q, left_on="_k", right_on="sample", how="inner")
        print(f"  matched {len(m):,} of {len(d):,} genomes to QC metrics")
        if m.empty:
            print("  No overlap — check that QC sample names match "
                  f"{key} values.")
            return 0
        num = ["n50_bp", "n_contigs", "genome_size_bp", "contamination_pct",
               "completeness_pct"]
        for c in num:
            if c in m.columns:
                m[c] = pd.to_numeric(m[c], errors="coerce")
        m["discordant_neg"] = (~m.cp) & (m.y == 1)
        m["discordant_pos"] = (m.cp) & (m.y == 0)

        rows = []
        for label, mask in (("resistant, carbapenemase-negative", m.discordant_neg),
                            ("susceptible, carbapenemase-positive", m.discordant_pos)):
            for c in num:
                if c not in m.columns or m[c].notna().sum() == 0:
                    continue
                x, yv = m.loc[mask, c].dropna(), m.loc[~mask, c].dropna()
                if len(x) < 5 or len(yv) < 5:
                    continue
                u, p = stats.mannwhitneyu(x, yv, alternative="two-sided")
                rows.append(dict(group=label, metric=c, n_discordant=len(x),
                                 median_discordant=round(float(x.median()), 1),
                                 n_other=len(yv),
                                 median_other=round(float(yv.median()), 1),
                                 p_value=p))
        if rows:
            qc = pd.DataFrame(rows)
            qc["q_value"] = bh(qc.p_value.tolist())
            print(f"\n{'group':36s}{'metric':18s}{'median disc.':>13s}"
                  f"{'median other':>14s}{'q':>10s}")
            for _, r in qc.iterrows():
                star = " *" if r.q_value < 0.05 else ""
                print(f"{r.group:36s}{r.metric:18s}"
                      f"{r.median_discordant:13.1f}{r.median_other:14.1f}"
                      f"{r.q_value:10.2e}{star}")
            p = a.outdir / f"discordance_vs_assembly_quality{a.tag}.tsv"
            qc.to_csv(p, sep="\t", index=False)
            print(f"\n-> {p}")
            print("\n  A significantly lower N50 or higher contig count among")
            print("  resistant carbapenemase-negative isolates is evidence that")
            print("  some of that cell is unrecovered genes rather than novel")
            print("  biology. Report it either way.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
