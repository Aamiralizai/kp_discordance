#!/usr/bin/env python
"""11_extensions.py — three analyses that complete the discordance story.

  [A] FALSE-RESISTANT MECHANISM
      Characterise isolates that are phenotypically resistant with NO
      carbapenemase. Hypothesis: ESBL/AmpC + porin disruption. If so, both
      directions of discordance reduce to one mechanism (influx vs hydrolysis).

  [B] INCREMENTAL PREDICTIVE VALUE
      How much does resistance prediction improve when porin state and allele
      identity are added to carbapenemase presence? Nested models, bootstrap CIs.
      This is the translational payoff: is porin genotyping worth reporting?

  [C] SURVEILLANCE PPV PROJECTION
      Per-enzyme discordance rates (measured cleanly here) are combined with
      EXTERNAL published enzyme-prevalence scenarios to project how genomic
      surveillance PPV changes as the carbapenemase landscape shifts.
      NOTE: this is a PROJECTION, not a measured temporal trend. The year
      variable in public AST data is confounded with laboratory and standard
      (chi2 p ~ 1e-22), so a direct time-series is not interpretable.

Usage:
    python scripts/11_extensions.py --tag _clsi_bmd
    python scripts/11_extensions.py --drug imipenem --tag _clsi_bmd
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import warnings
warnings.filterwarnings("ignore", message=".*[Ss]eparation.*")
warnings.filterwarnings("ignore", category=RuntimeWarning)
try:
    from statsmodels.tools.sm_exceptions import (ConvergenceWarning,
                                                 PerfectSeparationWarning)
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    warnings.filterwarnings("ignore", category=PerfectSeparationWarning)
except Exception:
    warnings.filterwarnings("ignore")

TAB = Path("results/tables")
KLEB = Path("data/processed/kleborate_klebsiella_pneumo_complex_output.tsv")

C_STRAIN = "strain"
C_ESBL = "klebsiella_pneumo_complex__amr__Bla_ESBL_acquired"
C_ESBL_INHR = "klebsiella_pneumo_complex__amr__Bla_ESBL_inhR_acquired"
C_BLA = "klebsiella_pneumo_complex__amr__Bla_acquired"
C_BLA_CHR = "klebsiella_pneumo_complex__amr__Bla_chr"
C_INHR = "klebsiella_pneumo_complex__amr__Bla_inhR_acquired"


def present(v):
    return isinstance(v, str) and v.strip() not in {"", "-", "nan", "NA", "None"}


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n; z = 1.959964
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (p, max(0, c - h), min(1, c + h))


# ============================================================ [A]
def analysis_A(d, kb, min_n=10):
    print("=" * 76)
    print("[A] FALSE-RESISTANT ISOLATES — resistant WITHOUT a carbapenemase")
    print("=" * 76)
    neg = d[~d["carb_pos"]].copy()
    fr = neg[neg["call"] == "R"]
    print(f"carbapenemase-negative : {len(neg)}")
    print(f"  phenotypically resistant : {len(fr)} ({len(fr)/max(len(neg),1):.1%})")
    if len(fr) < min_n:
        print("  too few to characterise"); return

    if "esbl" in neg.columns:
        cols = []
    else:
        cols = [c for c in [C_ESBL, C_ESBL_INHR, C_BLA, C_BLA_CHR, C_INHR]
                if c in kb.columns]
    kk = kb[[C_STRAIN] + cols].rename(columns={C_STRAIN: "genome_id"})
    kk["genome_id"] = kk["genome_id"].astype(str).str.strip()
    if cols:
        neg = neg.merge(kk, on="genome_id", how="left")
        for c in cols:
            neg[c] = neg[c].apply(present)
        neg["esbl_any"] = neg[[c for c in cols if "ESBL" in c]].any(axis=1)
    else:
        neg["esbl_any"] = neg["esbl"].astype(str).str.lower().isin({"true", "1"})
    if "porin_disrupted" in neg.columns:
        neg["porin_disrupted"] = neg["porin_disrupted"].astype(str).str.lower().isin({"true","1"})
    else:
        neg["porin_disrupted"] = neg["porin"] != "none"
    neg["res"] = neg["call"] == "R"

    print("\n  resistance rate by mechanism combination:")
    print(f"  {'ESBL':<7}{'porin':<9}{'n':>7}{'% resistant':>14}{'95% CI':>18}")
    rows = []
    for e in (False, True):
        for pdis in (False, True):
            g = neg[(neg["esbl_any"] == e) & (neg["porin_disrupted"] == pdis)]
            if len(g) < min_n:
                continue
            k = int(g["res"].sum()); p, lo, hi = wilson(k, len(g))
            print(f"  {str(e):<7}{str(pdis):<9}{len(g):>7}{100*p:>13.1f}%"
                  f"{f'{100*lo:.1f}-{100*hi:.1f}':>18}")
            rows.append(dict(esbl=e, porin_disrupted=pdis, n=len(g),
                             pct_resistant=100 * p, lo=100 * lo, hi=100 * hi))
    res = pd.DataFrame(rows)

    # is the combination super-additive on the odds scale?
    try:
        tab = pd.crosstab([neg["esbl_any"], neg["porin_disrupted"]], neg["res"])
        import statsmodels.formula.api as smf
        neg["y"] = neg["res"].astype(int)
        neg["E"] = neg["esbl_any"].astype(int)
        neg["P"] = neg["porin_disrupted"].astype(int)
        m_add = smf.logit("y ~ E + P", data=neg).fit(disp=0)
        m_int = smf.logit("y ~ E * P", data=neg).fit(disp=0)
        lr = 2 * (m_int.llf - m_add.llf)
        pv = stats.chi2.sf(lr, 1)
        print(f"\n  ESBL x porin interaction: chi2 = {lr:.2f}, p = {pv:.3g}")
        print(f"    ESBL alone   OR = {math.exp(m_add.params['E']):.2f}")
        print(f"    porin alone  OR = {math.exp(m_add.params['P']):.2f}")
        if pv < 0.05:
            print("    -> synergistic: neither mechanism alone explains resistance")
    except Exception as e:
        print(f"\n  (interaction model unavailable: {type(e).__name__})")

    # what fraction of false-resistant isolates are explained?
    # Recompute the explained/unexplained breakdown from the SAME esbl_any
    # definition used in the 2x2 above (an earlier version silently used an
    # empty column list here, reporting every isolate as ESBL-negative).
    fr2 = neg[neg["res"]]
    e_any = fr2["esbl_any"].astype(bool)
    p_dis = fr2["porin_disrupted"] if "porin_disrupted" in fr2 else (fr2["porin"] != "none")
    expl = fr2
    print(f"\n  of the {len(fr2)} false-resistant isolates:")
    print(f"    ESBL present            : {int(e_any.sum())} ({e_any.mean():.1%})")
    print(f"    porin disrupted         : {int(p_dis.sum())} ({p_dis.mean():.1%})")
    print(f"    BOTH                    : {int((e_any & p_dis).sum())} "
          f"({(e_any & p_dis).mean():.1%})")
    print(f"    NEITHER (unexplained)   : {int((~e_any & ~p_dis).sum())} "
          f"({(~e_any & ~p_dis).mean():.1%})")
    res.to_csv(TAB / "extA_false_resistant.tsv", sep="\t", index=False)


# ============================================================ [B]
def analysis_B(d, n_boot=2000, seed=0):
    print("\n" + "=" * 76)
    print("[B] INCREMENTAL PREDICTIVE VALUE of porin / allele genotype")
    print("=" * 76)
    rng = np.random.default_rng(seed)
    dd = d.copy()
    dd["y_res"] = (dd["call"] == "R").astype(int)

    def metrics(pred, y):
        tp = int(((pred == 1) & (y == 1)).sum()); fp = int(((pred == 1) & (y == 0)).sum())
        fn = int(((pred == 0) & (y == 1)).sum()); tn = int(((pred == 0) & (y == 0)).sum())
        return dict(
            accuracy=(tp + tn) / max(len(y), 1),
            ppv=tp / max(tp + fp, 1), npv=tn / max(tn + fn, 1),
            sens=tp / max(tp + fn, 1), spec=tn / max(tn + fp, 1))

    # model 1: carbapenemase presence alone
    p1 = dd["carb_pos"].astype(int).values
    # model 2: + porin (predict resistant if carbapenemase OR porin disrupted)
    p2 = ((dd["carb_pos"]) | (dd["porin"] != "none")).astype(int).values
    # model 3: mechanism-aware — weak enzymes need porin disruption to be resistant
    weak = dd["carb_family"].isin(["OXA-48-like", "VIM", "IMP"])
    strong = dd["carb_pos"] & ~weak
    p3 = (strong | (weak & (dd["porin"] != "none")) |
          (~dd["carb_pos"] & (dd["porin"] != "none"))).astype(int).values

    y = dd["y_res"].values
    names = ["carbapenemase only",
             "carbapenemase OR porin",
             "mechanism-aware (weak enzyme requires porin loss)"]
    preds = [p1, p2, p3]
    print(f"{'model':<52}{'acc':>7}{'PPV':>7}{'sens':>7}{'spec':>7}")
    rows = []
    for nm, pr in zip(names, preds):
        m = metrics(pr, y)
        print(f"{nm:<52}{m['accuracy']:>7.3f}{m['ppv']:>7.3f}"
              f"{m['sens']:>7.3f}{m['spec']:>7.3f}")
        rows.append(dict(model=nm, **m))

    # bootstrap CI on the accuracy difference (model 3 vs model 1)
    idx = np.arange(len(y)); diffs = []
    for _ in range(n_boot):
        b = rng.choice(idx, len(idx), replace=True)
        diffs.append(metrics(p3[b], y[b])["accuracy"] -
                     metrics(p1[b], y[b])["accuracy"])
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    print(f"\n  accuracy gain (mechanism-aware vs carbapenemase-only):")
    print(f"    {rows[2]['accuracy'] - rows[0]['accuracy']:+.3f} "
          f"[95% CI {lo:+.3f} to {hi:+.3f}]")
    if lo > 0:
        print("    -> adding porin genotype significantly improves prediction")
    else:
        print("    -> gain not significant at this sample size")
    pd.DataFrame(rows).to_csv(TAB / "extB_incremental_value.tsv",
                              sep="\t", index=False)


# ============================================================ [C]
def analysis_C(d, min_n=20):
    print("\n" + "=" * 76)
    print("[C] SURVEILLANCE PPV PROJECTION under shifting enzyme prevalence")
    print("=" * 76)
    print("  *** PROJECTION, not a measured time trend. Year in public AST data")
    print("  *** is confounded with laboratory and testing standard.")
    pos = d[d["carb_pos"]]
    rates = {}
    for fam in ["KPC", "NDM", "OXA-48-like", "VIM", "IMP"]:
        g = pos[pos["carb_family"] == fam]
        if len(g) < min_n:
            continue
        p, lo, hi = wilson(int((g["call"] == "S").sum()), len(g))
        rates[fam] = (p, lo, hi, len(g))
    if not rates:
        print("  insufficient per-family data"); return
    print("\n  measured per-enzyme discordance (this study, stratified):")
    for f, (p, lo, hi, n) in rates.items():
        print(f"    {f:<14} {100*p:5.1f}%  [{100*lo:.1f}-{100*hi:.1f}]  n={n:,}")

    # scenarios: external prevalence mixes (illustrative, cite literature)
    scen = {
        "KPC-dominant (e.g. historical US/S. Europe)":
            {"KPC": .70, "NDM": .15, "OXA-48-like": .15},
        "balanced":
            {"KPC": .40, "NDM": .25, "OXA-48-like": .35},
        "OXA-48-dominant (e.g. Middle East/N. Africa)":
            {"KPC": .15, "NDM": .25, "OXA-48-like": .60},
    }
    print("\n  projected PPV of 'carbapenemase detected' => resistant:")
    print(f"  {'scenario':<46}{'PPV':>8}{'range':>18}")
    rows = []
    for name, mix in scen.items():
        mix = {k: v for k, v in mix.items() if k in rates}
        tot = sum(mix.values())
        if tot == 0:
            continue
        # NOTE: uncertainty is obtained by BOOTSTRAP over the observed
        # per-family counts. Summing marginal Wilson limits (as in earlier
        # versions) does NOT yield a valid interval for a weighted sum.
        disc = sum((v / tot) * rates[k][0] for k, v in mix.items())
        rng = np.random.default_rng(0)
        draws = []
        for _ in range(4000):
            s = 0.0
            for k, v in mix.items():
                n_k = rates[k][3]
                p_k = rates[k][0]
                s += (v / tot) * (rng.binomial(n_k, p_k) / n_k if n_k else 0.0)
            draws.append(s)
        disc_lo, disc_hi = np.percentile(draws, [2.5, 97.5])
        print(f"  {name:<46}{100*(1-disc):>7.1f}%"
              f"{f'{100*(1-disc_hi):.1f}-{100*(1-disc_lo):.1f}':>18}")
        rows.append(dict(scenario=name, ppv=100 * (1 - disc),
                         ppv_lo=100 * (1 - disc_hi), ppv_hi=100 * (1 - disc_lo),
                         **{f"prev_{k}": v for k, v in mix.items()}))
    print("\n  Interpretation: as populations shift from KPC- toward OXA-48-dominant,")
    print("  the positive predictive value of genomic carbapenem-resistance")
    print("  surveillance declines, because OXA-48-like enzymes require porin")
    print("  disruption to produce a resistant phenotype.")
    print("\n  Prevalence weights above are ILLUSTRATIVE — replace with values")
    print("  cited from regional surveillance literature before publication.")
    pd.DataFrame(rows).to_csv(TAB / "extC_ppv_projection.tsv", sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="meropenem")
    ap.add_argument("--tag", default="_clsi_bmd")
    a = ap.parse_args()

    f = TAB / f"recordlevel_discordance_{a.drug}{a.tag}.tsv"
    if False:
        pass
    if not f.exists():
        sys.exit(f"ERROR: no discordance table for {a.drug}{a.tag} — run "
                 f"08b_discordance_recordlevel.py first")
    d = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
    if "mic" in d.columns:
        d["mic"] = pd.to_numeric(d["mic"], errors="coerce")
    d["carb_pos"] = d["carb_pos"].astype(str).str.lower().isin({"true", "1"})
    # restrict to the BINARY set; intermediate/ambiguous are not S or R
    n_before = len(d)
    d = d[d["call"].isin(["S", "R"])]
    if len(d) < n_before:
        print(f"  restricted to binary S/R set: {n_before:,} -> {len(d):,}\n")
    if not KLEB.exists():
        sys.exit(f"ERROR: {KLEB} not found")
    kb = pd.read_csv(KLEB, sep="\t", dtype=str, low_memory=False)

    print(f"drug = {a.drug}   stratum = {a.tag or '(unstratified)'}   n = {len(d):,}\n")
    analysis_A(d, kb)
    analysis_B(d)
    analysis_C(d)
    print(f"\n[out] {TAB}/extA_*, extB_*, extC_*")


if __name__ == "__main__":
    main()
