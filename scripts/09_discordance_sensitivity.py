#!/usr/bin/env python
"""09_discordance_sensitivity.py — stress-test the discordance finding.

Five checks, in order of how much damage each could do:

  1. PROVENANCE  — is discordance clustered in a few studies/labs/countries/years?
                   (If yes, we are measuring a laboratory, not biology.)
  2. BREAKPOINT  — CLSI (S<=1,R>=4) vs EUCAST (S<=2,R>8). Does the rate move?
  3. CLONALITY   — deduplicate ST x country x year. Outbreak clusters are not
                   independent observations; CIs must widen.
  4. LINEAGE     — logistic regression of susceptibility on porin state adjusted
                   for ST. Is the porin effect just ST258?
  5. MULTIPLICITY— Benjamini-Hochberg FDR across all family/allele tests.

Usage:
    python scripts/09_discordance_sensitivity.py
    python scripts/09_discordance_sensitivity.py --drug imipenem
"""
import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from kp_defs import BREAKPOINTS as KP_BREAKPOINTS
# Reuse the record-level classifier rather than maintaining a second
# implementation that can drift from it.
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location(
    "_rl", str(Path(__file__).resolve().parent / "08b_discordance_recordlevel.py"))
_rl = _ilu.module_from_spec(_spec)
try:
    _spec.loader.exec_module(_rl)
    classify_interval_shared = _rl.classify_interval
except Exception:
    classify_interval_shared = None
from kp_defs import UNVERIFIED_BREAKPOINTS

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

sys.path.insert(0, str(Path(__file__).parent))
# (the legacy 08_discordance_analysis module is archived as invalid and is
# deliberately NOT imported)

RAW_PHEN = Path("data/raw/phenotypes/kp_amr_anchor.tsv")
DISC = Path("results/tables")
OUTDIR = Path("results/tables")

def bp_sets_for(drug):
    """Drug-SPECIFIC breakpoints. An earlier version applied the meropenem
    EUCAST values (S<=2, R>8) to every drug, which is wrong: EUCAST imipenem
    uses a different resistant breakpoint."""
    out = {}
    d = drug.strip().lower()
    for std, block in KP_BREAKPOINTS.items():
        if d in block:
            out[std] = dict(block[d])
    if not out:
        print(f"  (no VERIFIED breakpoints for {drug}; see config/breakpoints.csv)")
    return out




def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n; z = 1.959964
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (p, max(0, c - h), min(1, c + h))


def classify(mic, lc, rc, bp_s, bp_r):
    if mic is None or (isinstance(mic, float) and math.isnan(mic)):
        return "ambiguous"
    if lc:
        return "S" if mic <= bp_s else "ambiguous"
    if rc:
        return "R" if mic >= bp_r else "ambiguous"
    if mic <= bp_s: return "S"
    if mic >= bp_r: return "R"
    return "I"


def _exp(x, cap=1e6):
    """math.exp with a cap. Under separation the point estimate and its CI can
    overflow; report a bound rather than crashing."""
    try:
        v = math.exp(x)
    except OverflowError:
        return math.inf
    return v if v < cap else math.inf


def _fmt_or(b, se):
    """Format OR and 95% CI, showing '>1e6' / '<1e-6' where unbounded."""
    def f(v):
        if v is math.inf or v > 1e6:
            return ">1e6"
        if v < 1e-6:
            return "<1e-6"
        return f"{v:.4f}"
    return f(_exp(b)), f(_exp(b - 1.96 * se)), f(_exp(b + 1.96 * se))


def firth_logit(X, y, max_iter=200, tol=1e-6):
    """Firth penalised logistic regression (Jeffreys-prior penalty).

    Standard maximum likelihood fails under complete or quasi-complete
    separation, which occurs here because some porin strata contain no
    susceptible isolates. Firth's method removes first-order bias and yields
    finite estimates in exactly these situations.
    """
    X = np.asarray(X, float); y = np.asarray(y, float)
    n, k = X.shape
    beta = np.zeros(k)
    for _ in range(max_iter):
        eta = X @ beta
        pr = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        W = pr * (1 - pr)
        XtWX = X.T @ (X * W[:, None])
        try:
            XtWX_inv = np.linalg.inv(XtWX + 1e-10 * np.eye(k))
        except np.linalg.LinAlgError:
            return None
        # hat values give the Jeffreys penalty term
        H = (X * W[:, None]) @ XtWX_inv @ X.T
        h = np.clip(np.diag(H), 0, 1)
        U = X.T @ (y - pr + h * (0.5 - pr))
        step = XtWX_inv @ U
        # damp large steps for stability
        nrm = np.linalg.norm(step)
        if nrm > 5:
            step *= 5 / nrm
        beta_new = beta + step
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new
    else:
        return None
    eta = X @ beta
    pr = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
    W = pr * (1 - pr)
    try:
        cov = np.linalg.inv(X.T @ (X * W[:, None]) + 1e-10 * np.eye(k))
    except np.linalg.LinAlgError:
        return None
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    return beta, se


def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    n = len(p); order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(ranked, 0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="meropenem")
    ap.add_argument("--min-n", type=int, default=10)
    ap.add_argument("--tag", default="_clsi_bmd",
                    help="table suffix written by 08b_discordance_recordlevel.py")
    a = ap.parse_args()

    disc_f = DISC / f"recordlevel_discordance_{a.drug}{a.tag}.tsv"
    if not disc_f.exists():
        sys.exit(f"\nMissing {disc_f}.\nRun 08b_discordance_recordlevel.py "
                 f"first. Legacy discordance_*.tsv outputs are INVALID (see "
                 f"archive/deprecated/) and are no longer used as a fallback.\n")
    print(f"input: {disc_f}")
    d = pd.read_csv(disc_f, sep="\t", dtype=str, low_memory=False)
    for c in ("mic", "lo", "hi"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    # Every boolean-like column must be converted: the table is read with
    # dtype=str, so a stored False arrives as the STRING "False", which is
    # truthy. Omitting lo_inc/hi_inc previously caused right-censored records
    # at a EUCAST resistant breakpoint to be discarded as ambiguous.
    for c in ["lc", "rc", "carb_pos", "lo_inc", "hi_inc", "porin_disrupted",
              "esbl"]:
        if c in d.columns:
            d[c] = d[c].astype(str).str.strip().str.lower().isin({"true", "1"})
    # KEEP the full set: the breakpoint comparison must reclassify ALL
    # isolates, not only those already binary under the reference standard.
    d_all = d.copy()
    n0 = len(d)
    d = d[d["call"].isin(["S", "R"])].copy()
    print(f"binary S/R set: {n0:,} -> {len(d):,} "
          f"(excluded {n0-len(d):,} intermediate/unresolvable under the "
          f"reference standard)")
    if not {"lo_inc", "hi_inc"}.issubset(d_all.columns):
        sys.exit(
            f"\nSTALE INPUT: {disc_f.name} has no lo_inc/hi_inc columns.\n"
            f"It was produced before interval-bound inclusivity was recorded.\n"
            f"Reclassifying it would treat '>8' as ambiguous rather than\n"
            f"resistant and silently discard right-censored records.\n\n"
            f"Regenerate it:\n"
            f"  python scripts/08b_discordance_recordlevel.py --drug {a.drug} \\\n"
            f"      --standard CLSI --method 'Broth dilution' --tag {a.tag}\n")
    d["susceptible"] = d["call"] == "S"
    pos = d[d["carb_pos"]].copy()
    global REG_ROWS, LOO_ROWS
    REG_ROWS, LOO_ROWS = [], []
    print("=" * 76)
    print(f"SENSITIVITY ANALYSIS — {a.drug.upper()}")
    print(f"carbapenemase-positive isolates: {len(pos)}")
    print("=" * 76)

    # ============ 1. PROVENANCE ============
    print("\n[1] PROVENANCE / BATCH EFFECTS")
    print("-" * 76)
    if RAW_PHEN.exists():
        raw = pd.read_csv(RAW_PHEN, sep="\t", dtype=str, low_memory=False).fillna("")
        abc = next((c for c in raw.columns if c.lower().strip() == "antibiotic"), None)
        gid = next((c for c in raw.columns if c.lower().strip() == "genome_id"), None)
        if abc and gid:
            raw = raw[raw[abc].str.lower().str.strip() == a.drug.lower()]
            meta_cols = [c for c in raw.columns if any(
                k in c.lower() for k in ["laboratory_typing", "testing_standard",
                                         "isolation_country", "collection_year",
                                         "bioproject", "pmid", "year"])]
            print(f"metadata columns found: {meta_cols if meta_cols else 'NONE'}")
            if meta_cols:
                raw = raw[[gid] + meta_cols].drop_duplicates(subset=[gid])
                raw = raw.rename(columns={gid: "genome_id"})
                m = pos.merge(raw, on="genome_id", how="left")
                globals()["_PROV"] = raw          # reuse for clonality dedup
                for c in meta_cols:
                    vc = m[c].fillna("").replace("", "<missing>")
                    top = vc.value_counts().head(8)
                    if len(top) < 2:
                        continue
                    print(f"\n  --- {c} ---")
                    rows, kept_levels = [], []
                    for lvl in top.index:
                        g = m[vc == lvl]            # match on the FULL label
                        if len(g) < a.min_n:
                            continue
                        kept_levels.append(lvl)     # keep full label for subsetting
                        k = int(g["susceptible"].sum())
                        p, lo, hi = wilson(k, len(g))
                        # truncate only for DISPLAY, never for matching
                        rows.append(dict(level=(str(lvl)[:26] + "..") if len(str(lvl)) > 28
                                         else str(lvl),
                                         n=len(g), n_susc=k,
                                         pct=100 * p, lo=100 * lo, hi=100 * hi))
                    if rows:
                        t = pd.DataFrame(rows)
                        print(t.round(1).to_string(index=False))
                        # chi2 across levels
                        sub = m[vc.isin(kept_levels)]
                        try:
                            ct = pd.crosstab(vc[sub.index], m.loc[sub.index, "susceptible"])
                            if ct.shape[0] > 1 and ct.shape[1] == 2:
                                chi2, pv, _, _ = stats.chi2_contingency(ct)
                                flag = "  <-- HETEROGENEOUS" if pv < 0.05 else ""
                                print(f"    chi2 p = {pv:.3g}{flag}")
                        except Exception:
                            pass
                print("\n  Interpretation: strong heterogeneity across labs/standards")
                print("  means part of the discordance is METHODOLOGICAL, not biological.")
        else:
            print("  could not locate antibiotic/genome_id columns in raw file")
    else:
        print(f"  {RAW_PHEN} not found — skipping provenance check")

    # ============ 2. BREAKPOINT SENSITIVITY ============
    bp_sets = bp_sets_for(a.drug)
    print("\n[2] BREAKPOINT STANDARD (CLSI vs EUCAST)")
    print(f"    breakpoints applied for {a.drug}: " +
          "; ".join(f"{k} S{v.get('S_op','<=')}{v['S']} R{v.get('R_op','>=')}{v['R']}"
                    for k, v in bp_sets.items()))
    print("-" * 76)
    if {"lo", "hi"}.issubset(d_all.columns):
        print("  (reclassifying the FULL set under each standard, not the")
        print("   already-binary subset — otherwise the comparison is biased)")
        def classify_iv(lo, hi, s, r, s_op="<=", r_op=">=",
                        lo_inc=True, hi_inc=True):
            """Delegates to the record-level classifier so the two cannot drift."""
            if pd.isna(lo) or pd.isna(hi):
                return "ambiguous"
            if classify_interval_shared is None:
                raise RuntimeError("record-level classifier unavailable")
            return classify_interval_shared((float(lo), float(hi),
                                             bool(lo_inc), bool(hi_inc)),
                                            s, r, s_op, r_op)

        rows = []
        for name, bp in bp_sets.items():
            _li = d_all["lo_inc"] if "lo_inc" in d_all else [True]*len(d_all)
            _hi = d_all["hi_inc"] if "hi_inc" in d_all else [True]*len(d_all)
            calls = [classify_iv(lo, hi, bp["S"], bp["R"],
                                 bp.get("S_op", "<="), bp.get("R_op", ">="), li, hi_)
                     for lo, hi, li, hi_ in zip(d_all["lo"], d_all["hi"], _li, _hi)]
            dd = d_all.copy(); dd["call2"] = calls
            bset = dd[dd["call2"].isin(["S", "R"])]
            p_ = bset[bset["carb_pos"]]; n_ = bset[~bset["carb_pos"]]
            fs = int((p_["call2"] == "S").sum()); fr = int((n_["call2"] == "R").sum())
            r1 = wilson(fs, len(p_)); r2 = wilson(fr, len(n_))
            rows.append(dict(standard=name, S=bp["S"], R=bp["R"],
                             n_binary=len(bset),
                             n_intermediate=int((dd["call2"] == "I").sum()),
                             n_ambiguous=int((dd["call2"] == "ambiguous").sum()),
                             gene_pos_n=len(p_),
                             pct_false_susceptible=round(100 * r1[0], 1),
                             fs_lo=round(100 * r1[1], 1), fs_hi=round(100 * r1[2], 1),
                             gene_neg_n=len(n_),
                             pct_false_resistant=round(100 * r2[0], 1),
                             fr_lo=round(100 * r2[1], 1), fr_hi=round(100 * r2[2], 1)))
        bt = pd.DataFrame(rows)
        print(bt.to_string(index=False))
        if len(bt) < 2:
            have = ", ".join(bt["standard"]) if len(bt) else "none"
            print(f"\n  !! Only ONE standard has verified breakpoints ({have}).")
            print("  !! The CLSI-versus-EUCAST comparison cannot be computed.")
            print("  !! Verify the missing rows in config/breakpoints.csv "
                  "(set verified=YES")
            print("  !! after checking the value AND the comparator against the")
            print("  !! document you will cite), then rerun.")
        else:
            delta = abs(bt["pct_false_susceptible"].iloc[0]
                        - bt["pct_false_susceptible"].iloc[1])
            print(f"\n  absolute shift in false-susceptible rate: "
                  f"{delta:.1f} percentage points")
            if delta > 3:
                print("  !! Breakpoint choice materially changes the headline.")
                print("  !! Report BOTH, and state the exact document versions used.")
            else:
                print("  Result is robust to breakpoint standard.")
        # per-family under each standard: the ordering matters more than the level
        print("\n  per-family false-susceptible rate under each standard:")
        for name, bp in bp_sets.items():
            _li = d_all["lo_inc"] if "lo_inc" in d_all else [True]*len(d_all)
            _hi = d_all["hi_inc"] if "hi_inc" in d_all else [True]*len(d_all)
            calls = [classify_iv(lo, hi, bp["S"], bp["R"],
                                 bp.get("S_op", "<="), bp.get("R_op", ">="), li, hi_)
                     for lo, hi, li, hi_ in zip(d_all["lo"], d_all["hi"], _li, _hi)]
            dd = d_all.copy(); dd["call2"] = calls
            sub = dd[dd["carb_pos"] & dd["call2"].isin(["S", "R"])]
            parts = []
            for fam, g in sub.groupby("carb_family"):
                if len(g) < a.min_n:
                    continue
                parts.append(f"{fam} {100*(g['call2']=='S').mean():.1f}% (n={len(g)})")
            print(f"    {name:<8}" + "; ".join(parts))
        bt.to_csv(OUTDIR / f"sensitivity_breakpoints_{a.drug}{a.tag}.tsv",
                  sep="\t", index=False)
    else:
        print("  (legacy table without interval bounds — rerun "
              "08b_discordance_recordlevel.py with --breakpoints EUCAST instead)")
        bt = pd.DataFrame()

    # ============ 3. CLONALITY ============
    print("\n[3] CLONALITY — deduplicating ST x country x year")
    print("-" * 76)
    prov = globals().get("_PROV")
    if prov is not None:
        pos = pos.merge(prov, on="genome_id", how="left", suffixes=("", "_prov"))
    key_cols = ["ST"]
    for extra in ["isolation_country", "collection_year", "country", "year"]:
        cc = [c for c in pos.columns if extra == c.lower()]
        if cc and cc[0] not in key_cols:
            key_cols.append(cc[0])
    pos["_clone_key"] = (pos[key_cols].fillna("NA").astype(str)
                         .apply(lambda r: "|".join(r), axis=1))
    dedup = pos.drop_duplicates(subset=["_clone_key", "carb_family", "porin"])
    k0 = int(pos["susceptible"].sum()); p0, lo0, hi0 = wilson(k0, len(pos))
    k1 = int(dedup["susceptible"].sum()); p1, lo1, hi1 = wilson(k1, len(dedup))
    print(f"  all isolates      : n = {len(pos):5d}  susceptible = {100*p0:.1f}% "
          f"[{100*lo0:.1f}-{100*hi0:.1f}]")
    print(f"  deduplicated      : n = {len(dedup):5d}  susceptible = {100*p1:.1f}% "
          f"[{100*lo1:.1f}-{100*hi1:.1f}]")
    print(f"  dedup key: {key_cols}")
    print(f"  CI width: {100*(hi0-lo0):.1f} -> {100*(hi1-lo1):.1f} percentage points")
    if abs(p0 - p1) > 0.05:
        print("  !! Estimate shifts >5pp after dedup — clonal expansion was inflating it.")
    else:
        print("  Estimate stable after deduplication.")

    # per-family after dedup
    rows = []
    for fam, g in dedup.groupby("carb_family"):
        if len(g) < a.min_n: continue
        k = int(g["susceptible"].sum()); p, lo, hi = wilson(k, len(g))
        rows.append(dict(family=fam, n=len(g), pct=100 * p, lo=100 * lo, hi=100 * hi))
    if rows:
        print("\n  per family (deduplicated):")
        print(pd.DataFrame(rows).sort_values("pct", ascending=False)
              .round(1).to_string(index=False))

    # ============ 4. LINEAGE ADJUSTMENT ============
    print("\n[4] LINEAGE ADJUSTMENT — is the porin effect just ST?")
    print("    exposure: OmpK36-level disruption (same as the record-level")
    print("    Fisher test); isolated OmpK35 loss is NOT counted as disrupted")
    print("-" * 76)
    try:
        import statsmodels.api as sm
        import statsmodels.formula.api as smf
        for fam in ["OXA-48-like", "KPC"]:
            sub = pos[pos["carb_family"] == fam].copy()
            n_family = len(sub)
            # SAME CONTRAST AS THE RECORD-LEVEL FISHER TEST: truly intact
            # porins versus OmpK36-level disruption. Isolates whose only change
            # is OmpK35 truncation are excluded from BOTH arms. Pooling them
            # into the reference group estimates a different quantity and gives
            # materially different odds ratios.
            sub = sub[(sub["porin"] == "none") | sub["porin_disrupted"]].copy()
            n_excluded = n_family - len(sub)
            if len(sub) < 30 or sub["susceptible"].nunique() < 2:
                print(f"  {fam}: insufficient data after restricting to the "
                      f"intact-vs-OmpK36 contrast "
                      f"(n = {len(sub)} of {n_family})")
                continue
            print(f"  {fam}  (intact vs OmpK36-disrupted; {n_excluded} "
                  f"OmpK35-only isolates excluded; n = {len(sub)} of {n_family})")
            sub["disrupted"] = sub["porin_disrupted"].astype(int)
            sub["y"] = sub["susceptible"].astype(int)
            top = sub["ST"].value_counts().head(10).index
            sub["ST_grp"] = np.where(sub["ST"].isin(top), sub["ST"], "other")
            def _record(model, label, orr, lo, hi, pv, method):
                REG_ROWS.append(dict(drug=a.drug, family=fam, model=label,
                                     odds_ratio=orr, ci_lo=lo, ci_hi=hi,
                                     p_value=pv, n=len(sub),
                                     events=int(sub["y"].sum()),
                                     n_excluded_ompK35_only=n_excluded,
                                     n_family=n_family, method=method))

            def firth_report(formula, label, n_ev):
                """Firth fallback; reports OR, 95% CI and p or states failure."""
                try:
                    import patsy
                    ymat, Xmat = patsy.dmatrices(formula, sub, return_type="dataframe")
                    names = list(Xmat.columns)
                    r = firth_logit(Xmat.values, ymat.values.ravel())
                    if r is None:
                        print(f"    {label:<16} NOT ESTIMABLE (Firth also failed)")
                        return None
                    beta, se = r
                    i = names.index("disrupted")
                    b, s = beta[i], se[i]
                    orr, lo, hi = _fmt_or(b, s)
                    z = b / s if s > 0 else np.nan
                    pv = 2 * stats.norm.sf(abs(z)) if s > 0 else np.nan
                    print(f"    {label:<16} Firth OR = {orr} "
                          f"(95% CI {lo}-{hi}), p = {pv:.3g}, "
                          f"n = {len(sub)}, events = {n_ev}")
                    _record(None, label, math.exp(b) if abs(b) < 700 else float("inf"),
                            math.exp(b - 1.96 * s) if abs(b) < 700 else float("nan"),
                            math.exp(b + 1.96 * s) if abs(b) < 700 else float("nan"),
                            pv, "Firth penalised")
                    return "firth"
                except (OverflowError, FloatingPointError, ValueError,
                        np.linalg.LinAlgError) as e:
                    print(f"    {label:<16} NOT ESTIMABLE "
                          f"(Firth failed: {type(e).__name__})")
                    return None
                except Exception as e:
                    print(f"    {label:<16} NOT ESTIMABLE "
                          f"(Firth failed: {type(e).__name__})")
                    return None

            def fit_report(formula, label):
                """Fit, CHECK CONVERGENCE, and report fully. A model that did
                not converge is reported as not estimable rather than having
                its coefficients printed."""
                try:
                    m = smf.logit(formula, data=sub).fit(disp=0)
                except Exception:
                    print(f"    {label:<16} ML failed (likely complete "
                          f"separation) — refitting with Firth")
                    return firth_report(formula, label, int(sub["y"].sum()))
                conv = bool(m.mle_retvals.get("converged", False))
                n_ev = int(sub["y"].sum())
                if not conv:
                    print(f"    {label:<16} ML did not converge "
                          f"(n={len(sub)}, events={n_ev}) — refitting with Firth")
                    return firth_report(formula, label, n_ev)
                b = m.params["disrupted"]; se = m.bse["disrupted"]
                orr, lo, hi = _fmt_or(b, se)
                print(f"    {label:<16} OR = {orr} (95% CI {lo}-{hi}), "
                      f"p = {m.pvalues['disrupted']:.3g}, "
                      f"n = {len(sub)}, events = {n_ev}, converged = yes")
                _record(m, label, _exp(b), _exp(b - 1.96 * se), _exp(b + 1.96 * se),
                        m.pvalues["disrupted"], "ML, converged")
                return m

            m0 = fit_report("y ~ disrupted", "unadjusted")
            m1 = fit_report("y ~ disrupted + C(ST_grp)", "ST-adjusted")
            if (m0 is not None and m1 is not None
                    and m0 != "firth" and m1 != "firth"):
                c0 = m0.params["disrupted"]; c1 = m1.params["disrupted"]
                if m1.pvalues["disrupted"] < 0.05 and np.sign(c1) == np.sign(c0):
                    print("                     -> association persists after "
                          "lineage adjustment")
                else:
                    print("                     -> association attenuates after "
                          "lineage adjustment")
            elif m1 is not None:
                print("                     -> adjusted estimate obtained by "
                      "Firth penalised regression (report as such)")
            else:
                print("                     -> adjusted estimate NOT reportable")
    except ImportError:
        print("  statsmodels not installed — pip install statsmodels")

    # ============ 4b. CLUSTERED AND LEAVE-ONE-OUT ANALYSES ============
    print("\n[4b] NON-INDEPENDENCE — clustered SEs and leave-one-out")
    print("-" * 76)
    try:
        import statsmodels.api as sm
        import statsmodels.formula.api as smf
        for fam in ("OXA-48-like", "KPC"):
            sub = pos[pos["carb_family"] == fam].copy()
            # SAME contrast as the record-level Fisher test: truly intact
            # porins versus OmpK36-level disruption. Isolates whose only
            # change is OmpK35 truncation are excluded from BOTH arms.
            n_family = len(sub)
            sub = sub[(sub["porin"] == "none") | sub["porin_disrupted"]].copy()
            if len(sub) < 30 or sub["susceptible"].nunique() < 2:
                print(f"  {fam}: insufficient data after restricting to the "
                      f"intact-vs-OmpK36 contrast (n={len(sub)} of {n_family})")
                continue
            print(f"  {fam}  (intact vs OmpK36-disrupted; {n_family - len(sub)} "
                  f"OmpK35-only isolates excluded; n = {len(sub)})")
            sub["y"] = sub["susceptible"].astype(int)
            sub["disrupted"] = sub["porin_disrupted"].astype(int)
            # cluster on the coarsest available proxy for study/collection
            cluster_col = None
            for cand in ("platform", "year", "standard", "ST"):
                if cand in sub.columns and sub[cand].nunique() > 1:
                    cluster_col = cand; break
            try:
                m_naive = smf.logit("y ~ disrupted", data=sub).fit(disp=0)
                naive_se = m_naive.bse["disrupted"]
                if cluster_col:
                    # groups must be an integer-coded array aligned to the
                    # model's own index, not a raw string Series
                    grp = pd.Categorical(sub[cluster_col].astype(str)).codes
                    m_cl = smf.logit("y ~ disrupted", data=sub).fit(
                        disp=0, cov_type="cluster",
                        cov_kwds={"groups": np.asarray(grp)})
                    cl_se = m_cl.bse["disrupted"]
                    print(f"  {fam:<14} coef = {m_naive.params['disrupted']:+.3f} | "
                          f"naive SE {naive_se:.3f} -> cluster SE {cl_se:.3f} "
                          f"(by {cluster_col}, {sub[cluster_col].nunique()} groups), "
                          f"p = {m_cl.pvalues['disrupted']:.3g}")
                    if cl_se > 1.5 * naive_se:
                        print("                 !! SEs inflate substantially — "
                              "non-independence matters here")
                else:
                    print(f"  {fam}: no clustering variable available")
            except Exception as e:
                print(f"  {fam}: clustered model failed ({type(e).__name__}) — "
                      f"often separation")
            # leave-one-ST-out
            if "ST" in sub.columns and sub["ST"].nunique() > 2:
                ors = []
                for st in sub["ST"].value_counts().head(8).index:
                    s2 = sub[sub["ST"] != st]
                    t2 = pd.crosstab(s2["disrupted"], s2["susceptible"])
                    if t2.shape == (2, 2):
                        orr, _ = stats.fisher_exact(t2.values)
                        ors.append((st, orr, len(s2)))
                if ors:
                    lo = min(o[1] for o in ors); hi = max(o[1] for o in ors)
                    print(f"                 leave-one-ST-out OR range: "
                          f"{lo:.3f}-{hi:.3f} over {len(ors)} exclusions")
                    for st, orr, n_ in ors:
                        LOO_ROWS.append(dict(drug=a.drug, family=fam,
                                             excluded_ST=st, odds_ratio=orr,
                                             n_remaining=n_))
    except ImportError:
        print("  statsmodels not installed")

    # ============ 5. MULTIPLICITY ============
    print("\n[5] MULTIPLE TESTING (Benjamini-Hochberg FDR)")
    print("    contrast: intact vs OmpK36-disrupted (OmpK35-only excluded),")
    print("    identical to the record-level Fisher tests")
    print("-" * 76)
    tests = []
    for fam, g in pos.groupby("carb_family"):
        # SAME restricted contrast as sections 3-4 and the record-level Fisher
        # test: intact versus OmpK36-disrupted, OmpK35-only excluded.
        g = g[(g["porin"] == "none") | g["porin_disrupted"]].copy()
        if len(g) < a.min_n: continue
        g["disrupted"] = g["porin_disrupted"]
        ct = pd.crosstab(g["disrupted"], g["susceptible"])
        if ct.shape == (2, 2):
            orr, pv = stats.fisher_exact(ct.values)
            tests.append(dict(test=f"{fam}: porin effect", n=len(g), OR=orr, p=pv))
    if tests:
        t = pd.DataFrame(tests)
        t["q_BH"] = bh_fdr(t["p"].values)
        t["significant_q05"] = t["q_BH"] < 0.05
        print(t.round(4).to_string(index=False))
        print(f"\n  {int(t['significant_q05'].sum())} of {len(t)} survive FDR correction")

    # save
    bt.to_csv(OUTDIR / f"sensitivity_breakpoints_{a.drug}.tsv", sep="\t", index=False)
    if tests:
        t.to_csv(OUTDIR / f"sensitivity_fdr_{a.drug}.tsv", sep="\t", index=False)
    if REG_ROWS:
        pd.DataFrame(REG_ROWS).to_csv(
            OUTDIR / f"sensitivity_regression_{a.drug}{a.tag}.tsv",
            sep="\t", index=False)
        print(f"    regression estimates -> sensitivity_regression_"
              f"{a.drug}{a.tag}.tsv ({len(REG_ROWS)} rows)")
    if LOO_ROWS:
        pd.DataFrame(LOO_ROWS).to_csv(
            OUTDIR / f"sensitivity_leaveoneout_{a.drug}{a.tag}.tsv",
            sep="\t", index=False)
        print(f"    leave-one-ST-out    -> sensitivity_leaveoneout_"
              f"{a.drug}{a.tag}.tsv ({len(LOO_ROWS)} rows)")
    print(f"\n[out] {OUTDIR}/sensitivity_*_{a.drug}{a.tag}.tsv")


if __name__ == "__main__":
    main()
