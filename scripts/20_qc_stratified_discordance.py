#!/usr/bin/env python3
"""20_qc_stratified_discordance.py — discordance by assembly QC verdict.

Exposes how each direction of genotype-phenotype discordance depends on the
PASS/WARN assembly-quality split, which metric drives the WARN flags, and the
lineage composition of the WARN carbapenemase-negative resistant subgroup.

Reads the locked primary record-level classification and the per-genome QC
table; writes:
  results/supplementary/S16_qc_stratified_discordance.tsv
  results/supplementary/S16b_warn_carbneg_profile.tsv
"""
from __future__ import annotations
import csv, math
from collections import Counter
from pathlib import Path

REC = Path("results/tables/recordlevel_discordance_meropenem_clsi_bmd.tsv")
QC = Path("data/interim/assembly_qc_metrics.tsv")
SUP = Path("results/supplementary"); SUP.mkdir(parents=True, exist_ok=True)

def wilson(k, n, z=1.96):
    if n == 0: return (0.0, 0.0, 0.0)
    p = k / n; d = 1 + z * z / n
    c = p + z * z / (2 * n); m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (100 * p, (c - m) / d * 100, (c + m) / d * 100)

def fisher(a, b, c, d):
    try:
        from scipy.stats import fisher_exact
        orr, p = fisher_exact([[a, b], [c, d]]); return orr, p
    except Exception:
        return float("nan"), float("nan")

qc = {r["sample"]: r for r in csv.DictReader(open(QC), delimiter="\t")}
def verdict(g): return qc.get(g, {}).get("verdict", "NO_METRIC")
rows = list(csv.DictReader(open(REC), delimiter="\t"))

def istrue(v): return str(v) in ("True", "1", "true")

out = [["direction", "qc_verdict", "n", "events", "pct", "ci_lo", "ci_hi"]]
# carbapenemase-positive susceptible
cp = [r for r in rows if istrue(r["carb_pos"]) and r["call"] in ("S", "R")]
# carbapenemase-negative resistant
cn = [r for r in rows if not istrue(r["carb_pos"]) and r["call"] in ("S", "R")]

def block(subset, direction, event_call):
    recs = []
    for g in ("PASS", "WARN", "PASS+WARN"):
        if g == "PASS+WARN":
            ss = subset
        else:
            ss = [r for r in subset if verdict(r["genome_id"]) == g]
        ev = [r for r in ss if r["call"] == event_call]
        p, lo, hi = wilson(len(ev), len(ss))
        out.append([direction, g, len(ss), len(ev), f"{p:.1f}", f"{lo:.1f}", f"{hi:.1f}"])
        recs.append((g, len(ev), len(ss)))
    return recs

cp_recs = block(cp, "carbapenemase-positive susceptible", "S")
cn_recs = block(cn, "carbapenemase-negative resistant", "R")

# Fisher PASS vs WARN for each direction
def pass_warn(recs):
    d = {g: (ev, n) for g, ev, n in recs}
    (we, wn), (pe, pn) = d["WARN"], d["PASS"]
    orr, p = fisher(we, wn - we, pe, pn - pe)
    return orr, p
cp_or, cp_p = pass_warn(cp_recs)
cn_or, cn_p = pass_warn(cn_recs)

with (SUP / "S16_qc_stratified_discordance.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerows(out)
    fh.write(f"# carbapenemase-positive susceptible PASS-vs-WARN Fisher OR={cp_or:.2f} p={cp_p:.2e}\n")
    fh.write(f"# carbapenemase-negative resistant PASS-vs-WARN Fisher OR={cn_or:.2f} p={cn_p:.2e}\n")

# WARN carbapenemase-negative resistant profile: warning metric, lineage, completeness
warn_cnr = [r for r in cn if r["call"] == "R" and verdict(r["genome_id"]) == "WARN"]
mcols = ["genome_size_bp", "gc_content_pct", "n50_bp", "n_contigs",
         "completeness_pct", "contamination_pct", "n_cds"]
mc = Counter()
for r in warn_cnr:
    q = qc[r["genome_id"]]
    for m in mcols:
        if q.get(m + "_verdict") == "WARN": mc[m] += 1
st_counts = Counter(r.get("ST", "?") for r in warn_cnr)
comp = sorted(float(qc[r["genome_id"]]["completeness_pct"]) for r in warn_cnr
              if qc[r["genome_id"]].get("completeness_pct"))
med_comp = comp[len(comp) // 2] if comp else float("nan")
with (SUP / "S16b_warn_carbneg_profile.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric_or_lineage", "warn_count_or_n"])
    w.writerow(["# WARN carbapenemase-negative resistant isolates", len(warn_cnr)])
    w.writerow(["# median CheckM completeness (%)", f"{med_comp:.1f}"])
    for m, n in mc.most_common():
        w.writerow([f"warn_metric:{m}", n])
    for s, n in st_counts.most_common(10):
        w.writerow([f"lineage_ST:{s}", n])

print("wrote S16_qc_stratified_discordance.tsv and S16b_warn_carbneg_profile.tsv")
print(f"carb-pos susceptible: PASS/WARN OR={cp_or:.2f} p={cp_p:.1e}")
print(f"carb-neg resistant:   PASS/WARN OR={cn_or:.2f} p={cn_p:.1e}")

# ---- S16c: carbapenemase-negative resistance by clonal group x QC verdict ----
# CG258 and CG101 are DISTINCT clonal groups. ST101 belongs to CG101, not CG258.
CG258 = {"258", "11", "512", "340", "437", "3107"}
def clonal_group(st):
    base = (st or "").replace("ST", "").split("-")[0]
    if base in CG258:
        return "CG258"
    if base == "101":
        return "CG101"
    return "other"
cn_all = [r for r in rows if not istrue(r["carb_pos"]) and r["call"] in ("S", "R")]
s16c = [["clonal_group", "qc_verdict", "n", "resistant", "pct", "ci_lo", "ci_hi"]]
for grp in ("CG258", "CG101", "other"):
    for v in ("PASS", "WARN"):
        ss = [r for r in cn_all if clonal_group(r.get("ST", "")) == grp
              and verdict(r["genome_id"]) == v]
        rr = [r for r in ss if r["call"] == "R"]
        p, lo, hi = wilson(len(rr), len(ss))
        s16c.append([grp, v, len(ss), len(rr), f"{p:.1f}", f"{lo:.1f}", f"{hi:.1f}"])
with (SUP / "S16c_carbneg_by_lineage_qc.tsv").open("w", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(s16c)

# ---- missing-QC selection-bias check (C2) ----
# Compare carbapenemase carriage of typed genomes that lack a QC metric row
# against all typed genomes, using the pre-filter Kleborate table if present.
from pathlib import Path as _P
pf = _P("data/processed/kleborate_klebsiella_pneumo_complex_output.prefilter.tsv")
if pf.exists():
    PAMR = "klebsiella_pneumo_complex__amr__Bla_Carb_acquired"
    qc_ids = set(qc)
    kb = list(csv.DictReader(pf.open(), delimiter="\t"))
    def is_carb(r):
        v = (r.get(PAMR, "") or "").strip()
        return bool(v) and v not in ("-",)
    miss = [r for r in kb if r["strain"] not in qc_ids]
    n_miss = len(miss); n_all = len(kb)
    miss_carb = sum(1 for r in miss if is_carb(r))
    all_carb = sum(1 for r in kb if is_carb(r))
    with (SUP / "missing_qc_composition.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["group", "n", "carbapenemase_positive", "carb_pos_pct"])
        w.writerow(["typed_without_QC_metric", n_miss, miss_carb,
                    f"{100*miss_carb/n_miss:.1f}"])
        w.writerow(["all_typed", n_all, all_carb, f"{100*all_carb/n_all:.1f}"])
    print(f"missing-QC: {n_miss} genomes, {100*miss_carb/n_miss:.1f}% carb+ vs "
          f"{100*all_carb/n_all:.1f}% in all typed")
print("wrote S16c_carbneg_by_lineage_qc.tsv and missing_qc_composition.tsv")

# ---- S16d: is the carbapenemase-positive PASS/WARN difference compositional? --
# After the truncated-download repair the carbapenemase-positive direction, not
# the carbapenemase-negative one, differs by QC verdict (13.2% PASS vs 6.6%
# WARN). Susceptibility in that direction is driven by enzyme family and porin
# status, so the question is whether PASS and WARN assemblies simply differ in
# their family, porin and lineage composition. This fits the crude model and
# then the same model with those covariates added.
try:
    import pandas as pd, statsmodels.api as sm
    cp = [r for r in rows if istrue(r["carb_pos"]) and r["call"] in ("S", "R")
          and verdict(r["genome_id"]) in ("PASS", "WARN")]
    d = pd.DataFrame(dict(
        y=[1 if r["call"] == "S" else 0 for r in cp],
        warn=[1 if verdict(r["genome_id"]) == "WARN" else 0 for r in cp],
        family=[r["carb_family"] or "none" for r in cp],
        porin=[r["porin"] or "none" for r in cp],
        cg=[clonal_group(r.get("ST", "")) for r in cp]))
    out = [["model", "term", "odds_ratio", "ci_lo", "ci_hi", "p_value", "n"]]
    for name, cols in (("crude", []),
                       ("adjusted for family, porin and clonal group",
                        ["family", "porin", "cg"])):
        X = pd.get_dummies(d[["warn"] + cols], columns=cols, drop_first=True)
        X = sm.add_constant(X.astype(float))
        res = sm.Logit(d["y"].astype(float), X).fit(disp=0, method="bfgs", maxiter=500)
        b, se = res.params["warn"], res.bse["warn"]
        out.append([name, "WARN vs PASS", f"{math.exp(b):.2f}",
                    f"{math.exp(b - 1.96 * se):.2f}", f"{math.exp(b + 1.96 * se):.2f}",
                    f"{res.pvalues['warn']:.3g}", len(d)])
    with (SUP / "S16d_carbpos_qc_adjusted.tsv").open("w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(out)
    for r in out[1:]:
        print(f"carb-pos WARN vs PASS, {r[0]}: OR {r[2]} ({r[3]}-{r[4]}) p={r[5]} n={r[6]}")
    print("wrote S16d_carbpos_qc_adjusted.tsv")
except Exception as e:                                   # pragma: no cover
    print(f"S16d skipped: {e}")
