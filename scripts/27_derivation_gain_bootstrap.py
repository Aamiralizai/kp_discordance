#!/usr/bin/env python3
"""27_derivation_gain_bootstrap.py — does the porin rule's gain survive clonality?

WHY THIS EXISTS
The derivation set reports that adding porin status within enzyme families
raises positive predictive value from 90.0 to 95.4 percent and accuracy from
92.5 to 93.4. Those are cross-validated with folds split by sequence type, which
stops a rule being scored on its own training clone. It does not put an interval
on the gain itself, and the gain is the claim the paper makes.

Public collections are clonal, so isolates of one sequence type are not
independent observations. This resamples whole sequence types with replacement
and reports a percentile interval for the paired difference the porin rule makes
on the same isolates, using the out-of-fold predictions from the same
sequence-type-grouped cross-validation the paper reports. It is the derivation
counterpart of the BioProject bootstrap in 19_external_validation.py, so the
internal and external increments are stated on the same footing.

WHAT THIS INTERVAL DOES AND DOES NOT COVER
It resamples clusters around a fixed set of out-of-fold predictions produced by
one fold assignment. It therefore covers the sampling of sequence types, and it
does not refit the rules inside each replicate, so it does not propagate the
uncertainty of rule fitting or the arbitrariness of the fold split. --repeats
addresses the second of those directly: it repeats the whole grouped
cross-validation under a different seed each time and reports the spread of the
gain across fold assignments, which says whether the point estimate is an
artefact of one split.

USAGE (from the project root)
    python3 scripts/27_derivation_gain_bootstrap.py --drug meropenem

OUTPUT
    results/tables/derivation_gain_bootstrap{tag}.tsv
"""
import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

spec = importlib.util.spec_from_file_location(
    "rule_derivation", ROOT / "scripts" / "17_rule_derivation.py")
rd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rd)


def oof_predictions(df, tier, threshold, n_splits, seed):
    """Out-of-fold predictions under sequence-type-grouped cross-validation."""
    from sklearn.model_selection import StratifiedGroupKFold
    y = df["y"].to_numpy()
    chain = rd.chain_for(tier)
    pred = np.empty(len(df), dtype=int)
    sgk = StratifiedGroupKFold(n_splits=n_splits, shuffle=True,
                               random_state=seed)
    for tr, te in sgk.split(df, y, df["sequence_type"]):
        dtr, dte = df.iloc[tr], df.iloc[te]
        fitted = {t: rd.fit_rules(dtr, rd.TIERS[t], threshold) for t in chain}
        pred[te] = rd.predict(dte, fitted, threshold,
                              int(dtr["y"].mean() > threshold))
    return pred


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--drug", default="meropenem")
    ap.add_argument("--data",
                    default="results/tables/recordlevel_discordance_meropenem_clsi_bmd.tsv")
    ap.add_argument("--tag", default="")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260101)   # 17_rule_derivation.py default
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--repeats", type=int, default=20,
                    help="repeat the whole grouped CV under this many seeds, "
                         "to show the gain is not an artefact of one fold split")
    a = ap.parse_args()

    path = ROOT / a.data
    if not path.exists():
        sys.exit(f"27_derivation_gain_bootstrap: missing {path}")
    d = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    call = "call" if "call" in d.columns else f"{a.drug}_interpretation"
    d = d[d[call].isin(["S", "R"])].copy()
    d["y"] = (d[call] == "R").astype(int)
    d["carb_pos"] = (d["carb_pos"].astype(str).str.lower()
                     .isin({"true", "1", "yes"}).astype(int))
    d["family"] = [rd.family_level(f, bool(p)) for f, p
                   in zip(d["carb_family"], d["carb_pos"])]
    d["porin"] = d["porin"].fillna("unknown").map(rd.porin_level)
    d["sequence_type"] = d["ST"].fillna("unknown")
    d["lineage"] = np.where(d["sequence_type"].isin(rd.CG258),
                            "CG258", "non-CG258")
    d = d.reset_index(drop=True)
    print(f"  binary set {len(d):,} genomes, "
          f"{d['sequence_type'].nunique():,} sequence types")

    tiers = ("T0_gene_only", "T2_family_porin")
    preds = {t: oof_predictions(d, t, a.threshold, a.folds, a.seed)
             for t in tiers}
    y = d["y"].to_numpy()
    groups = d["sequence_type"].to_numpy()
    idx_by_group = {g: np.flatnonzero(groups == g) for g in np.unique(groups)}
    keys = list(idx_by_group)
    rng = np.random.default_rng(a.seed)

    def score(sel):
        out = {}
        for name, p in (("a", preds[tiers[0]][sel]), ("b", preds[tiers[1]][sel])):
            yy = y[sel]
            tp = int(((p == 1) & (yy == 1)).sum()); fp = int(((p == 1) & (yy == 0)).sum())
            tn = int(((p == 0) & (yy == 0)).sum())
            out[f"{name}_acc"] = 100.0 * (tp + tn) / len(yy) if len(yy) else np.nan
            out[f"{name}_ppv"] = 100.0 * tp / (tp + fp) if tp + fp else np.nan
        return out

    obs = score(np.arange(len(y)))
    diffs = {"accuracy": [], "PPV": []}
    for _ in range(a.boot):
        draw = rng.choice(len(keys), size=len(keys), replace=True)
        sel = np.concatenate([idx_by_group[keys[i]] for i in draw])
        m = score(sel)
        if np.isnan(m["a_ppv"]) or np.isnan(m["b_ppv"]):
            continue
        diffs["accuracy"].append(m["b_acc"] - m["a_acc"])
        diffs["PPV"].append(m["b_ppv"] - m["a_ppv"])

    rows = []
    for metric, k in (("accuracy", "acc"), ("PPV", "ppv")):
        v = np.asarray(diffs[metric])
        rows.append(dict(
            metric=metric, rule_a=tiers[0], rule_b=tiers[1],
            value_a=round(obs[f"a_{k}"], 1), value_b=round(obs[f"b_{k}"], 1),
            difference=round(obs[f"b_{k}"] - obs[f"a_{k}"], 1),
            ci_lo=round(float(np.percentile(v, 2.5)), 1),
            ci_hi=round(float(np.percentile(v, 97.5)), 1),
            n_sequence_types=len(keys), n_isolates=int(len(y)),
            n_resamples=int(len(v)),
            method="sequence-type cluster bootstrap on grouped-CV predictions"))
    # fold sensitivity: repeat the entire grouped cross-validation under
    # different seeds, so the reported gain is not conditional on one split
    if a.repeats > 1:
        rep = []
        for i in range(a.repeats):
            pr = {t: oof_predictions(d, t, a.threshold, a.folds, a.seed + i)
                  for t in tiers}
            yy = y
            vals = {}
            for name, p in (("a", pr[tiers[0]]), ("b", pr[tiers[1]])):
                tp = int(((p == 1) & (yy == 1)).sum())
                fp = int(((p == 1) & (yy == 0)).sum())
                tn = int(((p == 0) & (yy == 0)).sum())
                vals[f"{name}_acc"] = 100.0 * (tp + tn) / len(yy)
                vals[f"{name}_ppv"] = 100.0 * tp / (tp + fp) if tp + fp else np.nan
            rep.append(dict(seed=a.seed + i,
                            accuracy_gain=vals["b_acc"] - vals["a_acc"],
                            PPV_gain=vals["b_ppv"] - vals["a_ppv"]))
        rp = pd.DataFrame(rep)
        rp.to_csv(ROOT / "results" / "tables" /
                  f"derivation_gain_fold_sensitivity{a.tag}.tsv",
                  sep="\t", index=False)
        print(f"\n  Fold sensitivity across {a.repeats} grouped-CV splits:")
        for col, lab in (("PPV_gain", "PPV gain"),
                         ("accuracy_gain", "accuracy gain")):
            v = rp[col]
            print(f"    {lab:<14} median {v.median():+.1f}  "
                  f"range {v.min():+.1f} to {v.max():+.1f}  "
                  f"splits favouring the porin rule {int((v > 0).sum())}/{len(v)}")

    out = ROOT / "results" / "tables" / f"derivation_gain_bootstrap{a.tag}.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print("\n  Gain from the family-specific porin rule, derivation set:")
    for r in rows:
        print(f"    {r['metric']:<9} {r['value_a']:5.1f} -> {r['value_b']:5.1f}  "
              f"difference {r['difference']:+5.1f} "
              f"(95% CI {r['ci_lo']:+.1f} to {r['ci_hi']:+.1f})")
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
