#!/usr/bin/env python3
"""
17_rule_derivation.py — predictive rules from carbapenemase, porin and lineage.

The discordance analysis establishes that genotype and phenotype disagree in
both directions, and that the disagreement depends on carbapenemase family and
porin state. That is an association result. It does not tell a laboratory what
to predict from a given genotype, which is the question that decides whether
any of this is usable.

This script answers that question directly. It defines nested rule sets of
increasing complexity, and for each one reports the predictive values a
laboratory would actually obtain.

RULE TIERS

    T0  carbapenemase present            the current default: gene -> resistant
    T1  carbapenemase family
    T2  carbapenemase family x porin class
    T3  carbapenemase family x porin class x lineage (CG258 membership)

Each tier is a decision list. A stratum predicts resistance when the resistance
rate among its training isolates exceeds --threshold. A stratum absent from the
training data falls back to the prediction of its parent tier, so every isolate
receives a prediction and no tier can quietly refuse the hard cases.

WHY THE VALIDATION IS STRUCTURED THIS WAY

A rule fitted and scored on the same isolates is optimistic, and the finer
strata here are small enough for that optimism to matter. Two cross-validated
estimates are therefore reported alongside the apparent one.

    stratified k-fold   folds balanced on outcome. Measures optimism from
                        fitting stratum-specific thresholds.

    grouped k-fold      folds split by sequence type, so no ST appears in both
                        training and test. Public genome collections are
                        heavily clonal; without this, a rule can be scored on
                        near-identical isolates of a clone it was fitted on,
                        which inflates every statistic. This is the estimate to
                        quote.

The gap between the two is itself a result: it quantifies how much of the
rules' apparent performance is clonal redundancy rather than transferable
signal.

OUTPUT

    results/tables/rules_strata{tag}.tsv        every stratum, with Wilson CIs
    results/tables/rules_performance{tag}.tsv   tier performance, 3 estimates
    results/tables/rules_decision_list{tag}.tsv the fitted rules themselves
    results/tables/rules_comparison{tag}.tsv    each tier against the T0 baseline

USAGE

    python scripts/17_rule_derivation.py --drug meropenem
    python scripts/17_rule_derivation.py --drug imipenem --tag _imipenem
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TAB = Path("results/tables")
SUP = Path("results/supplementary")

# Clonal group 258: the globally dominant KPC-associated K. pneumoniae lineage.
# Grouping its single-locus variants prevents the lineage tier from splitting
# one clonal group across several sparse strata.
CG258 = {"ST258", "ST11", "ST340", "ST437", "ST512", "ST418", "ST833",
         "ST3984", "ST4496"}

# OmpK36-level disruption. Isolated OmpK35 loss is kept as its own level rather
# than pooled with these: pooling it was the error that produced a diluted
# porin effect in an earlier version of this analysis.
OMPK36_DISRUPTED = {"ompK36_loop3", "ompK36_loss", "ompK36_reduced_expression"}


# --------------------------------------------------------------- statistics
def wilson(k: int, n: int, z: float = 1.959963984540054):
    """Wilson score interval. Degrades gracefully at n = 0 and at k = 0 or n."""
    if n == 0:
        return (float("nan"), float("nan"), float("nan"))
    p = k / n
    d = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = (z / d) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (p, max(0.0, centre - half), min(1.0, centre + half))


def fmt_ci(k: int, n: int) -> tuple:
    p, lo, hi = wilson(k, n)
    return (round(100 * p, 1), round(100 * lo, 1), round(100 * hi, 1))


def metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    """Predictive values for a binary rule. y = 1 means resistant."""
    tp = int(((y_pred == 1) & (y_true == 1)).sum())
    fp = int(((y_pred == 1) & (y_true == 0)).sum())
    tn = int(((y_pred == 0) & (y_true == 0)).sum())
    fn = int(((y_pred == 0) & (y_true == 1)).sum())
    out = dict(n=len(y_true), TP=tp, FP=fp, TN=tn, FN=fn)
    for name, k, d in (("PPV", tp, tp + fp), ("NPV", tn, tn + fn),
                       ("sensitivity", tp, tp + fn), ("specificity", tn, tn + fp),
                       ("accuracy", tp + tn, len(y_true))):
        p, lo, hi = fmt_ci(k, d)
        out[name], out[f"{name}_lo"], out[f"{name}_hi"] = p, lo, hi
        out[f"{name}_denom"] = d
    return out


# ------------------------------------------------------------------ feature
def porin_level(v: str) -> str:
    if v in OMPK36_DISRUPTED:
        return "ompK36_disrupted"
    if v == "ompK35_loss":
        return "ompK35_only"
    if v == "none":
        return "intact"
    return "unknown"


def family_level(v: str, present: bool) -> str:
    if not present or v in ("none", "", None) or pd.isna(v):
        return "none"
    v = str(v)
    if "+" in v:
        return "multiple"
    if v in ("KPC", "OXA-48-like", "NDM", "VIM"):
        return v
    return "other"


TIERS = {
    "T0_gene_only": ["carb_pos"],
    "T1_family": ["family"],
    "T2_family_porin": ["family", "porin"],
    "T3_family_porin_lineage": ["family", "porin", "lineage"],
}
# Each tier falls back to the one before it when a stratum is unseen.
PARENT = {"T0_gene_only": None,
          "T1_family": "T0_gene_only",
          "T2_family_porin": "T1_family",
          "T3_family_porin_lineage": "T2_family_porin"}


def fit_rules(df: pd.DataFrame, keys: list[str], threshold: float) -> dict:
    """Map stratum -> (predict_resistant, n, n_resistant) from training data."""
    rules = {}
    for key, g in df.groupby(keys, dropna=False, observed=True):
        if not isinstance(key, tuple):
            key = (key,)
        n = len(g)
        nr = int(g["y"].sum())
        rules[key] = (1 if (nr / n) > threshold else 0, n, nr)
    return rules


def predict(df: pd.DataFrame, fitted: dict, threshold: float,
            prevalence_fallback: int) -> np.ndarray:
    """Apply the tier hierarchy, falling back to coarser tiers then prevalence."""
    order = ["T3_family_porin_lineage", "T2_family_porin",
             "T1_family", "T0_gene_only"]
    out = np.empty(len(df), dtype=int)
    for i, (_, row) in enumerate(df.iterrows()):
        pred = None
        for tier in order:
            if tier not in fitted:
                continue
            key = tuple(row[k] for k in TIERS[tier])
            if key in fitted[tier]:
                pred = fitted[tier][key][0]
                break
        out[i] = prevalence_fallback if pred is None else pred
    return out


def chain_for(tier: str) -> list[str]:
    """The tier itself plus every coarser tier it can fall back to."""
    chain, t = [], tier
    while t is not None:
        chain.append(t)
        t = PARENT[t]
    return chain


def evaluate_tier(df: pd.DataFrame, tier: str, threshold: float,
                  n_splits: int, seed: int) -> tuple[dict, pd.DataFrame]:
    """Apparent and three cross-validated estimates for one tier.

    Three resampling schemes are reported because they isolate different
    sources of optimism, and the differences between them are informative.

      cv_stratified         Outcome-balanced folds, isolates split at random.
                            Measures only the optimism of fitting
                            stratum-specific decisions. Ignores clonality, so
                            it is an upper bound.

      cv_grouped_by_ST      Folds split by sequence type AND balanced on
                            outcome (StratifiedGroupKFold). No ST appears in
                            both training and test, so a rule cannot be scored
                            on near-identical isolates of a clone it was fitted
                            on. This is the estimate to quote.

      cv_grouped_unbalanced Folds split by sequence type with no outcome
                            balancing (plain GroupKFold). Reported because in
                            a heavily clonal collection this lets whole
                            resistant clones land in a single fold: observed
                            fold prevalence ranges roughly 0.25-0.62 here
                            against a pooled 0.38. That swing is a property of
                            the resampling scheme, not of the rules, so this
                            figure is a pessimistic bound rather than an
                            estimate of transferability.

    Also returns a per-stratum stability table: how often each stratum's
    decision flipped across the grouped folds. A rule that flips is a rule
    whose direction depends on which clones happened to be sampled.
    """
    chain = chain_for(tier)
    y = df["y"].to_numpy()
    base_rate = int(y.mean() > threshold)

    fitted_full = {t: fit_rules(df, TIERS[t], threshold) for t in chain}
    app = metrics(y, predict(df, fitted_full, threshold, base_rate))

    from sklearn.model_selection import (GroupKFold, StratifiedGroupKFold,
                                         StratifiedKFold)

    flips: dict[tuple, int] = {}
    n_folds_seen: dict[tuple, int] = {}

    def cv(splitter, groups=None, record=False):
        pred = np.empty(len(df), dtype=int)
        it = (splitter.split(df, y, groups) if groups is not None
              else splitter.split(df, y))
        for tr, te in it:
            dtr, dte = df.iloc[tr], df.iloc[te]
            f = {t: fit_rules(dtr, TIERS[t], threshold) for t in chain}
            pred[te] = predict(dte, f, threshold,
                               int(dtr["y"].mean() > threshold))
            if record:
                for key, (p, _, _) in f[tier].items():
                    full = fitted_full[tier].get(key)
                    if full is None:
                        continue
                    n_folds_seen[key] = n_folds_seen.get(key, 0) + 1
                    if p != full[0]:
                        flips[key] = flips.get(key, 0) + 1
        return metrics(y, pred)

    out = {
        "apparent": app,
        "cv_stratified": cv(StratifiedKFold(n_splits=n_splits, shuffle=True,
                                            random_state=seed)),
        "cv_grouped_by_ST": cv(StratifiedGroupKFold(n_splits=n_splits,
                                                    shuffle=True,
                                                    random_state=seed),
                               groups=df["sequence_type"], record=True),
        "cv_grouped_unbalanced": cv(GroupKFold(n_splits=n_splits),
                                    groups=df["sequence_type"]),
    }

    stab = pd.DataFrame([
        dict(tier=tier,
             stratum=" | ".join(f"{k}={v}" for k, v
                                in zip(TIERS[tier], key)),
             n=fitted_full[tier][key][1],
             n_resistant=fitted_full[tier][key][2],
             resistance_pct=round(100 * fitted_full[tier][key][2]
                                  / fitted_full[tier][key][1], 1),
             decision=("RESISTANT" if fitted_full[tier][key][0]
                       else "SUSCEPTIBLE"),
             folds_present=n_folds_seen.get(key, 0),
             folds_decision_flipped=flips.get(key, 0),
             stable=("no" if flips.get(key, 0) else "yes"))
        for key in fitted_full[tier]
    ]).sort_values(["folds_decision_flipped", "n"], ascending=[False, False])

    return out, stab


def evaluate_subset(df: pd.DataFrame, tier: str, threshold: float,
                    n_splits: int, seed: int, label: str) -> list[dict]:
    """Grouped-CV performance restricted to one carbapenemase stratum.

    Pooling carbapenemase-positive and -negative isolates hides the structure
    that matters: a rule can be highly transferable in one and not in the
    other, and a single pooled figure reports neither.
    """
    if len(df) < n_splits * 2 or df["y"].nunique() < 2:
        return []
    from sklearn.model_selection import StratifiedGroupKFold
    chain = chain_for(tier)
    y = df["y"].to_numpy()
    pred = np.empty(len(df), dtype=int)
    sgk = StratifiedGroupKFold(n_splits=n_splits, shuffle=True,
                               random_state=seed)
    try:
        splits = list(sgk.split(df, y, df["sequence_type"]))
    except ValueError:
        return []
    for tr, te in splits:
        dtr, dte = df.iloc[tr], df.iloc[te]
        f = {t: fit_rules(dtr, TIERS[t], threshold) for t in chain}
        pred[te] = predict(dte, f, threshold, int(dtr["y"].mean() > threshold))
    rows = [dict(subset=label, tier=tier, estimate="cv_grouped_by_ST",
                 **metrics(y, pred))]
    fitted = {t: fit_rules(df, TIERS[t], threshold) for t in chain}
    rows.append(dict(subset=label, tier=tier, estimate="apparent",
                     **metrics(y, predict(df, fitted, threshold,
                                          int(y.mean() > threshold)))))
    return rows


# --------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path,
                    default=SUP / "S1_genome_level_data.tsv")
    ap.add_argument("--drug", default="meropenem", choices=["meropenem", "imipenem"])
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="stratum resistance rate above which the rule predicts R")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260101)
    ap.add_argument("--min-support", type=int, default=10,
                    help="strata below this n are reported but flagged")
    ap.add_argument("--tag", default="")
    ap.add_argument("--outdir", type=Path, default=TAB)
    a = ap.parse_args()

    if not a.data.exists():
        sys.exit(f"Missing {a.data}.\nRun the record-level pipeline first; this "
                 f"script never recomputes the analysis set.")

    d = pd.read_csv(a.data, sep="\t", dtype=str, low_memory=False)
    col = f"{a.drug}_interpretation"
    if col not in d.columns:
        sys.exit(f"{a.data} has no column {col!r}.")

    # binary analysis set only — intermediate and unresolvable isolates have no
    # S/R outcome and must not be imputed into one.
    d = d[d[col].isin(["S", "R"])].copy()
    d["y"] = (d[col] == "R").astype(int)

    d["carb_pos"] = (d["carbapenemase_present"].astype(str).str.lower()
                     .isin({"true", "1", "yes"}).astype(int))
    d["family"] = [family_level(f, bool(p)) for f, p
                   in zip(d["carbapenemase_family"], d["carb_pos"])]
    d["porin"] = d["porin_class"].fillna("unknown").map(porin_level)
    d["sequence_type"] = d["sequence_type"].fillna("unknown")
    d["lineage"] = np.where(d["sequence_type"].isin(CG258), "CG258", "non-CG258")

    print(f"drug            : {a.drug}")
    print(f"binary set      : {len(d):,} genomes "
          f"({int(d['y'].sum()):,} resistant, {len(d) - int(d['y'].sum()):,} susceptible)")
    print(f"sequence types  : {d['sequence_type'].nunique():,}")
    print(f"rule threshold  : stratum resistance rate > {a.threshold:.0%}")
    print(f"folds           : {a.folds}")

    a.outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------- stratum-level table
    rows = []
    for tier, keys in TIERS.items():
        for key, g in d.groupby(keys, dropna=False, observed=True):
            if not isinstance(key, tuple):
                key = (key,)
            n, nr = len(g), int(g["y"].sum())
            rate, lo, hi = fmt_ci(nr, n)
            rows.append(dict(
                tier=tier,
                stratum=" | ".join(f"{k}={v}" for k, v in zip(keys, key)),
                **{k: v for k, v in zip(keys, key)},
                n=n, n_resistant=nr, n_susceptible=n - nr,
                resistance_pct=rate, ci_lower=lo, ci_upper=hi,
                predicts_resistant=int((nr / n) > a.threshold),
                low_support=("yes" if n < a.min_support else "no"),
                n_sequence_types=g["sequence_type"].nunique()))
    strata = pd.DataFrame(rows).sort_values(
        ["tier", "n"], ascending=[True, False])
    p = a.outdir / f"rules_strata{a.tag}.tsv"
    strata.to_csv(p, sep="\t", index=False)
    print(f"\nstrata table    : {p}  ({len(strata)} strata)")

    # -------------------------------------------------- tier performance
    perf_rows, dl_rows, stab_all, sub_rows = [], [], [], []
    pos = d[d.carb_pos == 1]
    neg = d[d.carb_pos == 0]
    for tier in TIERS:
        res, stab = evaluate_tier(d, tier, a.threshold, a.folds, a.seed)
        stab_all.append(stab)
        for estimate, m in res.items():
            perf_rows.append(dict(tier=tier, estimate=estimate,
                                  n_strata=int((strata.tier == tier).sum()), **m))
        sub_rows += evaluate_subset(pos, tier, a.threshold, a.folds, a.seed,
                                    "carbapenemase_positive")
        sub_rows += evaluate_subset(neg, tier, a.threshold, a.folds, a.seed,
                                    "carbapenemase_negative")
        fitted = fit_rules(d, TIERS[tier], a.threshold)
        for key, (pred, n, nr) in sorted(fitted.items(),
                                         key=lambda kv: -kv[1][1]):
            rate, lo, hi = fmt_ci(nr, n)
            dl_rows.append(dict(
                tier=tier,
                rule=" AND ".join(f"{k}={v}" for k, v
                                  in zip(TIERS[tier], key)),
                then_predict="RESISTANT" if pred else "SUSCEPTIBLE",
                n=n, n_resistant=nr, resistance_pct=rate,
                ci_lower=lo, ci_upper=hi,
                low_support=("yes" if n < a.min_support else "no")))

    perf = pd.DataFrame(perf_rows)
    p = a.outdir / f"rules_performance{a.tag}.tsv"
    perf.to_csv(p, sep="\t", index=False)
    print(f"performance     : {p}")

    dl = pd.DataFrame(dl_rows)
    p = a.outdir / f"rules_decision_list{a.tag}.tsv"
    dl.to_csv(p, sep="\t", index=False)
    print(f"decision list   : {p}")

    stability = pd.concat(stab_all, ignore_index=True)
    p = a.outdir / f"rules_stability{a.tag}.tsv"
    stability.to_csv(p, sep="\t", index=False)
    print(f"stability       : {p}")

    subsets = pd.DataFrame(sub_rows)
    p = a.outdir / f"rules_by_carbapenemase_status{a.tag}.tsv"
    subsets.to_csv(p, sep="\t", index=False)
    print(f"subset analysis : {p}")

    # ------------------------------------------- comparison against T0
    ESTIMATES = ("apparent", "cv_stratified", "cv_grouped_by_ST",
                 "cv_grouped_unbalanced")
    comp = []
    for estimate in ESTIMATES:
        base = perf[(perf.tier == "T0_gene_only")
                    & (perf.estimate == estimate)].iloc[0]
        for tier in TIERS:
            r = perf[(perf.tier == tier) & (perf.estimate == estimate)].iloc[0]
            comp.append(dict(
                estimate=estimate, tier=tier,
                PPV=r.PPV, delta_PPV=round(r.PPV - base.PPV, 1),
                NPV=r.NPV, delta_NPV=round(r.NPV - base.NPV, 1),
                sensitivity=r.sensitivity,
                delta_sensitivity=round(r.sensitivity - base.sensitivity, 1),
                specificity=r.specificity,
                delta_specificity=round(r.specificity - base.specificity, 1),
                accuracy=r.accuracy,
                delta_accuracy=round(r.accuracy - base.accuracy, 1),
                misclassified=int(r.FP + r.FN),
                misclassified_baseline=int(base.FP + base.FN),
                net_reclassified_correctly=int((base.FP + base.FN)
                                               - (r.FP + r.FN))))
    comp = pd.DataFrame(comp)
    p = a.outdir / f"rules_comparison{a.tag}.tsv"
    comp.to_csv(p, sep="\t", index=False)
    print(f"comparison      : {p}")

    # ------------------------------------------------------------ console
    print("\n" + "=" * 78)
    print("TIER PERFORMANCE  (%, 95% Wilson CI)")
    print("=" * 78)
    for tier in TIERS:
        print(f"\n{tier}")
        for estimate in ESTIMATES:
            r = perf[(perf.tier == tier) & (perf.estimate == estimate)].iloc[0]
            print(f"  {estimate:21s} "
                  f"PPV {r.PPV:5.1f} [{r.PPV_lo:4.1f}-{r.PPV_hi:4.1f}]  "
                  f"NPV {r.NPV:5.1f} [{r.NPV_lo:4.1f}-{r.NPV_hi:4.1f}]  "
                  f"sens {r.sensitivity:5.1f}  spec {r.specificity:5.1f}  "
                  f"acc {r.accuracy:5.1f}  "
                  f"err {int(r.FP + r.FN):4d}")

    print("\n" + "=" * 78)
    print("GAIN OVER T0 (gene-only), ST-grouped cross-validation")
    print("=" * 78)
    sub = comp[comp.estimate == "cv_grouped_by_ST"]
    for _, r in sub.iterrows():
        print(f"  {r.tier:26s} dPPV {r.delta_PPV:+6.1f}  dNPV {r.delta_NPV:+6.1f}  "
              f"dacc {r.delta_accuracy:+6.1f}  "
              f"net correctly reclassified {r.net_reclassified_correctly:+5d}")

    print("\n" + "=" * 78)
    print("PERFORMANCE SPLIT BY CARBAPENEMASE STATUS (ST-grouped CV)")
    print("=" * 78)
    if not subsets.empty:
        for label in ("carbapenemase_positive", "carbapenemase_negative"):
            s = subsets[(subsets.subset == label)
                        & (subsets.estimate == "cv_grouped_by_ST")]
            if s.empty:
                continue
            print(f"\n  {label}  (n = {int(s.iloc[0].n):,})")
            for _, r in s.iterrows():
                print(f"    {r.tier:26s} PPV {r.PPV:5.1f} "
                      f"[{r.PPV_lo:4.1f}-{r.PPV_hi:4.1f}]  "
                      f"NPV {r.NPV:5.1f} [{r.NPV_lo:4.1f}-{r.NPV_hi:4.1f}]  "
                      f"acc {r.accuracy:5.1f}")

    print("\n" + "=" * 78)
    print("RULE STABILITY ACROSS ST-GROUPED FOLDS (T2)")
    print("=" * 78)
    st2 = stability[stability.tier == "T2_family_porin"]
    unstable = st2[st2.folds_decision_flipped > 0]
    if unstable.empty:
        print("  every T2 stratum kept its decision in every fold")
    else:
        print(f"  {'stratum':48s} {'n':>5s} {'R%':>6s} {'flips':>6s}  decision")
        for _, r in unstable.iterrows():
            print(f"  {r.stratum:48s} {r.n:5d} {r.resistance_pct:6.1f} "
                  f"{r.folds_decision_flipped:3d}/{r.folds_present:<2d}  "
                  f"{r.decision}")
        print("\n  A stratum that flips is one whose rate sits near the "
              f"{a.threshold:.0%} decision\n  boundary, so which clones are "
              "sampled determines the rule's direction.\n  Such a stratum "
              "should be reported as its observed rate with a confidence\n"
              "  interval, not as a dichotomous rule.")

    print("\n" + "=" * 78)
    print("STRATA WITH THE STRONGEST RESISTANCE SIGNAL (T2, n >= "
          f"{a.min_support})")
    print("=" * 78)
    t2 = strata[(strata.tier == "T2_family_porin")
                & (strata.n >= a.min_support)].sort_values(
        "resistance_pct", ascending=False)
    print(f"  {'stratum':48s} {'n':>5s} {'R%':>6s}  95% CI")
    for _, r in t2.iterrows():
        print(f"  {r.stratum:48s} {r.n:5d} {r.resistance_pct:6.1f}  "
              f"[{r.ci_lower:.1f}-{r.ci_upper:.1f}]")

    n_low = int((strata.low_support == "yes").sum())
    if n_low:
        print(f"\n{n_low} strata have n < {a.min_support} and are flagged "
              f"low_support in the output. Rules resting on them are not "
              f"reportable as predictive values.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
