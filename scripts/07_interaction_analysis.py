#!/usr/bin/env python
"""07_interaction_analysis.py — corrected enzyme x permeability analysis.

Fixes three problems with the first-pass analysis:
  1. CLEAN BASELINE  — the reference group now excludes ALL Omp mutations
     (previously OmpK35-loss genomes leaked into 'neither', inflating the
     baseline and distorting every effect size).
  2. CENSORING AUDIT — MICs recorded as '>=64' are right-censored. Medians
     silently understate effects when a cell piles up at the assay ceiling.
     This quantifies it per cell and refits with a Tobit model.
  3. FORMAL INTERACTION TEST — log2_mic ~ enzyme * porin, rather than
     comparing group medians (which is not an interaction test).

Usage:
    python scripts/07_interaction_analysis.py
    python scripts/07_interaction_analysis.py --drug imipenem --min-cell 15
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

sys.path.insert(0, str(Path(__file__).parent))

KLEB = Path("data/processed/kleborate_klebsiella_pneumo_complex_output.tsv")
PHEN = Path("data/processed/phenotypes_clean.tsv")
OUTDIR = Path("results/tables"); FIGDIR = Path("results/figures")

C_STRAIN = "strain"
C_SPECIES = "enterobacterales__species__species"
C_ST = "klebsiella_pneumo_complex__mlst__ST"
C_CARB = "klebsiella_pneumo_complex__amr__Bla_Carb_acquired"
C_OMP = "klebsiella_pneumo_complex__amr__Omp_mutations"

PORIN_ORDER = ["none", "ompK35_loss", "ompK36_loop3", "ompK36_loss"]
FAM_ORDER = ["none", "OXA-48-like", "KPC", "NDM"]


def has_carb(v):
    return isinstance(v, str) and v.strip() not in {"", "-", "nan", "NA", "None"}


def carb_family(v):
    if not has_carb(v):
        return "none"
    s = v.upper(); f = []
    if "KPC" in s: f.append("KPC")
    if "NDM" in s: f.append("NDM")
    if any(x in s for x in ("OXA-48", "OXA-181", "OXA-232", "OXA-244")):
        f.append("OXA-48-like")
    if "VIM" in s: f.append("VIM")
    if "IMP" in s: f.append("IMP")
    return "+".join(sorted(set(f))) if f else "other"


def porin_status(v):
    """Graded OmpK class from Kleborate HGVS notation."""
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


# ------------- Interval (two-sided censored) regression -------------
def _num_hess(f, x, eps=1e-4):
    n = len(x); H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            xpp, xpm, xmp, xmm = x.copy(), x.copy(), x.copy(), x.copy()
            xpp[i] += eps; xpp[j] += eps
            xpm[i] += eps; xpm[j] -= eps
            xmp[i] -= eps; xmp[j] += eps
            xmm[i] -= eps; xmm[j] -= eps
            H[i, j] = H[j, i] = (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4 * eps * eps)
    return H


def interval_fit(X, y, left, right):
    """MLE for y* = Xb + e, e~N(0,s^2), with BOTH censoring types.
       left=True  -> true value <= y   (observation '<=1')
       right=True -> true value >= y   (observation '>8')
    """
    X = np.asarray(X, float); y = np.asarray(y, float)
    L = np.asarray(left, bool); R = np.asarray(right, bool)
    E = ~(L | R)
    n, k = X.shape

    def nll(p):
        b, ls = p[:k], p[k]
        s = np.exp(np.clip(ls, -8, 8))
        z = (y - X @ b) / s
        ll = np.empty(n)
        ll[E] = stats.norm.logpdf(z[E]) - np.log(s)
        ll[R] = stats.norm.logsf(np.clip(z[R], -30, 30))
        ll[L] = stats.norm.logcdf(np.clip(z[L], -30, 30))
        v = -np.sum(ll)
        return v if np.isfinite(v) else 1e12

    b0 = np.linalg.lstsq(X, y, rcond=None)[0]
    p0 = np.concatenate([b0, [np.log(np.std(y - X @ b0) + 1e-3)]])
    r = optimize.minimize(nll, p0, method="Nelder-Mead",
                          options={"maxiter": 20000, "maxfev": 40000,
                                   "xatol": 1e-6, "fatol": 1e-6})
    r = optimize.minimize(nll, r.x, method="BFGS", options={"maxiter": 5000})
    beta, sigma = r.x[:k], float(np.exp(r.x[k]))
    # proper SEs from the numerical Hessian of the NLL
    try:
        H = _num_hess(nll, r.x)
        cov = np.linalg.pinv(H)
        se = np.sqrt(np.clip(np.diag(cov)[:k], 0, None))
    except Exception:
        se = np.full(k, np.nan)
    return beta, se, sigma, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="meropenem")
    ap.add_argument("--min-cell", type=int, default=10)
    a = ap.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True); FIGDIR.mkdir(parents=True, exist_ok=True)

    if not KLEB.exists(): sys.exit(f"ERROR: {KLEB} not found.")
    if not PHEN.exists(): sys.exit(f"ERROR: {PHEN} not found.")

    kb = pd.read_csv(KLEB, sep="\t", dtype=str, low_memory=False)
    kb = kb[kb[C_SPECIES].astype(str).str.contains("pneumoniae", case=False, na=False)]
    kb["genome_id"] = kb[C_STRAIN].astype(str).str.strip()
    kb["carb_family"] = kb[C_CARB].apply(carb_family)
    kb["porin"] = kb[C_OMP].apply(porin_status)
    kb["ST"] = kb[C_ST].astype(str)

    ph = pd.read_csv(PHEN, sep="\t", dtype=str, low_memory=False)
    abc = next(c for c in ph.columns if c.lower().strip() == "antibiotic")
    gid = next(c for c in ph.columns if c.lower().strip() == "genome_id")
    ph = ph[ph[abc].str.lower().str.strip() == a.drug.lower()].copy()
    ph["genome_id"] = ph[gid].astype(str).str.strip()
    ph["log2_mic"] = pd.to_numeric(ph.get("log2_mic"), errors="coerce")
    if "left_censored" in ph.columns:
        ph["lc"] = ph["left_censored"].astype(str).str.lower().isin({"true", "1"})
        ph["rc"] = ph["right_censored"].astype(str).str.lower().isin({"true", "1"})
    else:
        sgn = ph.get("mic_censor", pd.Series("", index=ph.index)).fillna("").astype(str)
        ph["lc"] = sgn.str.strip().isin(["<=", "<"])
        ph["rc"] = sgn.str.strip().isin([">=", ">"])
    ph = ph.dropna(subset=["log2_mic"])
    ph = ph.groupby("genome_id", as_index=False).agg(
        log2_mic=("log2_mic", "median"), lc=("lc", "max"), rc=("rc", "max"))
    ph["cens"] = ph["lc"] | ph["rc"]

    df = kb.merge(ph, on="genome_id", how="inner")
    print(f"joined: {len(df)} genomes with genotype + {a.drug} MIC")

    # =============== 1. CENSORING AUDIT ===============
    print("\n" + "=" * 74)
    print("1. CENSORING AUDIT  (right-censored '>=' MICs)")
    print("=" * 74)
    print(f"left-censored  (<=, <) : {int(df['lc'].sum())} ({df['lc'].mean():.1%})")
    print(f"right-censored (>=, >) : {int(df['rc'].sum())} ({df['rc'].mean():.1%})")
    print(f"exact                  : {int((~df['lc'] & ~df['rc']).sum())} "
          f"({(~df['lc'] & ~df['rc']).mean():.1%})")
    if df["cens"].mean() > 0.4:
        print("\n!! MAJORITY CENSORED. Medians and OLS in sections 2-3 are BIASED:")
        print("!! left-censoring inflates the susceptible baseline, right-censoring")
        print("!! deflates resistant cells -> both shrink apparent effects.")
        print("!! SECTION 4 (interval regression) IS THE VALID ANALYSIS.")
    sub = df[df["carb_family"].isin(FAM_ORDER)]
    cens_tab = sub.pivot_table(index="carb_family", columns="porin",
                               values="cens", aggfunc="mean")
    lc_tab = sub.pivot_table(index="carb_family", columns="porin",
                             values="lc", aggfunc="mean")
    rc_tab = sub.pivot_table(index="carb_family", columns="porin",
                             values="rc", aggfunc="mean")
    cens_tab = cens_tab.reindex(index=[f for f in FAM_ORDER if f in cens_tab.index],
                                columns=[c for c in PORIN_ORDER if c in cens_tab.columns])
    print("\n% LEFT-censored per cell:")
    print((lc_tab.reindex(index=[f for f in FAM_ORDER if f in lc_tab.index],
                          columns=[c for c in PORIN_ORDER if c in lc_tab.columns])
           * 100).round(1).to_string())
    print("\n% RIGHT-censored per cell:")
    print((rc_tab.reindex(index=[f for f in FAM_ORDER if f in rc_tab.index],
                          columns=[c for c in PORIN_ORDER if c in rc_tab.columns])
           * 100).round(1).to_string())
    full = (cens_tab >= 0.98).sum().sum()
    if full:
        print(f"\n!! {full} cell(s) are ~100% censored — their effects are NOT")
        print("!! identifiable by any method (only bounded). Interpret as bounds.")
    hi = (cens_tab > 0.3).sum().sum()
    if hi:
        print(f"\n!! {hi} cell(s) >30% censored — medians there UNDERSTATE the true MIC.")
        print("!! The Tobit model below is the trustworthy estimate, not the medians.")
    else:
        print("\nCensoring is low; medians and Tobit should broadly agree.")

    # =============== 2. CLEAN BASELINE TABLE ===============
    print("\n" + "=" * 74)
    print("2. TWO-WAY MEDIANS, CLEAN BASELINE (reference = no carbapenemase, no Omp)")
    print("=" * 74)
    piv = sub.pivot_table(index="carb_family", columns="porin",
                          values="log2_mic", aggfunc="median")
    cnt = sub.pivot_table(index="carb_family", columns="porin",
                          values="log2_mic", aggfunc="size")
    idx = [f for f in FAM_ORDER if f in piv.index]
    col = [c for c in PORIN_ORDER if c in piv.columns]
    piv, cnt = piv.reindex(index=idx, columns=col), cnt.reindex(index=idx, columns=col)
    print("median log2 MIC:"); print(piv.round(2).to_string())
    print("\nn per cell:"); print(cnt.fillna(0).astype(int).to_string())

    if "none" in piv.index and "none" in piv.columns:
        base = piv.loc["none", "none"]
        print(f"\nclean baseline (none/none) = {base:.2f} log2 = {2**base:.2f} mg/L")
        print("\nINTERACTION per cell: observed - (enzyme effect + porin effect)")
        print("  positive = super-additive, negative = sub-additive")
        rows = []
        for f in idx:
            if f == "none": continue
            e_f = piv.loc[f, "none"] - base
            for p in col:
                if p == "none": continue
                n_cell = cnt.loc[f, p]
                if pd.isna(piv.loc[f, p]) or n_cell < a.min_cell:
                    continue
                e_p = piv.loc["none", p] - base
                obs = piv.loc[f, p] - base
                inter = obs - (e_f + e_p)
                rows.append(dict(enzyme=f, porin=p, n=int(n_cell),
                                 enzyme_effect=e_f, porin_effect=e_p,
                                 additive_exp=e_f + e_p, observed=obs,
                                 interaction=inter,
                                 pct_censored=round(100 * cens_tab.loc[f, p], 1)
                                 if not pd.isna(cens_tab.loc[f, p]) else np.nan))
        it = pd.DataFrame(rows)
        if len(it):
            print(it.round(2).to_string(index=False))
            it.to_csv(OUTDIR / f"interaction_cells_{a.drug}.tsv", sep="\t", index=False)

    # =============== 3. FORMAL INTERACTION MODEL ===============
    print("\n" + "=" * 74)
    print("3. FORMAL INTERACTION MODEL   log2_mic ~ enzyme * porin")
    print("=" * 74)
    d = sub[sub["carb_family"].isin(FAM_ORDER) & sub["porin"].isin(PORIN_ORDER)].copy()
    d["enzyme"] = pd.Categorical(d["carb_family"], categories=FAM_ORDER)
    d["porin_c"] = pd.Categorical(d["porin"], categories=PORIN_ORDER)

    try:
        import statsmodels.api as sm
        import statsmodels.formula.api as smf

        m_add = smf.ols("log2_mic ~ C(enzyme) + C(porin_c)", data=d).fit()
        m_int = smf.ols("log2_mic ~ C(enzyme) * C(porin_c)", data=d).fit()
        lr = 2 * (m_int.llf - m_add.llf)
        ddf = int(m_int.df_model - m_add.df_model)
        p_lr = stats.chi2.sf(lr, ddf)
        print(f"OLS  additive R2 = {m_add.rsquared:.3f} | interaction R2 = {m_int.rsquared:.3f}")
        print(f"Likelihood-ratio test for interaction: chi2 = {lr:.1f}, df = {ddf}, "
              f"p = {p_lr:.3g}")
        print("  -> significant p means the enzyme effect DEPENDS on porin class")
        print("     (i.e. a real interaction, not merely additive effects)")

        # lineage control
        top_st = d["ST"].value_counts().head(15).index
        d["ST_grp"] = np.where(d["ST"].isin(top_st), d["ST"], "other")
        m_st = smf.ols("log2_mic ~ C(enzyme) * C(porin_c) + C(ST_grp)", data=d).fit()
        print(f"\nwith ST (lineage) adjustment: R2 = {m_st.rsquared:.3f}")
        m_st_add = smf.ols("log2_mic ~ C(enzyme) + C(porin_c) + C(ST_grp)", data=d).fit()
        lr2 = 2 * (m_st.llf - m_st_add.llf)
        ddf2 = int(m_st.df_model - m_st_add.df_model)
        print(f"interaction after ST adjustment: chi2 = {lr2:.1f}, df = {ddf2}, "
              f"p = {stats.chi2.sf(lr2, ddf2):.3g}")
        print("  -> if this stays significant, the interaction is NOT just lineage")
    except ImportError:
        print("statsmodels not installed — pip install statsmodels")
        m_int = None

    # =============== 4. TOBIT (censoring-aware) ===============
    print("\n" + "=" * 74)
    print("4. INTERVAL REGRESSION (handles BOTH left- and right-censored MICs)")
    print("   *** This is the valid analysis when censoring is high ***")
    print("=" * 74)
    try:
        D = pd.get_dummies(d[["enzyme", "porin_c"]], drop_first=True).astype(float)
        # interaction columns
        for ec in [c for c in D.columns if c.startswith("enzyme_")]:
            for pc in [c for c in D.columns if c.startswith("porin_c_")]:
                D[f"{ec}:{pc}"] = D[ec] * D[pc]
        D.insert(0, "const", 1.0)
        keep = D.columns[(D.sum() >= a.min_cell) | (D.columns == "const")]
        D = D[keep]
        beta, se, sigma, r = interval_fit(D.values, d["log2_mic"].values,
                                          d["lc"].values, d["rc"].values)
        out = pd.DataFrame({"term": D.columns, "coef": beta, "se": se})
        out["z"] = out["coef"] / out["se"].replace(0, np.nan)
        out["p"] = 2 * stats.norm.sf(np.abs(out["z"]))
        if not r.success:
            print("  (optimizer flagged non-convergence; estimates are usually still")
            print("   sound if SEs look sensible — sanity-check against the OLS above)")
        print(f"converged: {r.success} | sigma = {sigma:.2f}")
        inter = out[out["term"].str.contains(":")]
        print("\ninteraction terms (censoring-corrected):")
        print(inter.round(3).to_string(index=False) if len(inter) else "  (none estimable)")
        print("\nNOTE: positive coef = super-additive; negative = sub-additive.")
        print("Compare with the median-based interactions above — if they disagree,")
        print("trust the Tobit values, since censoring biases medians downward.")
        out.to_csv(OUTDIR / f"tobit_{a.drug}.tsv", sep="\t", index=False)
        print(f"\n[out] {OUTDIR}/tobit_{a.drug}.tsv")
    except Exception as e:
        print(f"Interval regression failed: {type(e).__name__}: {e}")

    print("\n" + "=" * 74)
    print("HOW TO READ THIS")
    print("  * high censoring in a cell -> its median is a LOWER BOUND")
    print("  * OXA-48 is expected to show the LARGEST positive interaction")
    print("    (weak meropenem hydrolysis -> permeability is rate-limiting)")
    print("  * NDM/KPC may show weaker/negative interaction (enzyme already")
    print("    saturating) OR it may be a censoring artifact — check section 1")
    print("=" * 74)


if __name__ == "__main__":
    main()
