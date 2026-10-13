#!/usr/bin/env python3
"""21_table4_determinants.py — regenerate Table 4 from the locked primary set.

Table 4 reports five beta-lactam determinants among carbapenemase-negative
isolates in the primary meropenem stratum, each with the number present, the
resistant proportion when present, an unadjusted odds ratio against its absence,
and a Fisher p. It also reports the ESBL-by-OmpK36 interaction.

OmpK36 disruption and acquired ESBL are taken from the locked record-level
classification (so they match the main analysis). Chromosomal beta-lactamase,
SHV mutations and inhibitor-resistant beta-lactamase are read from the Kleborate
output and joined on genome id. Writes:
  results/tables/table4_determinants.tsv
  results/tables/table4_esbl_ompk36_interaction.tsv
"""
from __future__ import annotations
import csv, math
from pathlib import Path

REC = Path("results/tables/recordlevel_discordance_meropenem_clsi_bmd.tsv")
KB = Path("data/processed/kleborate_klebsiella_pneumo_complex_output.tsv")
TAB = Path("results/tables")

P = "klebsiella_pneumo_complex__amr__"
COL_CHR = P + "Bla_chr"
COL_SHV = P + "SHV_mutations"
COL_INHR = [P + "Bla_inhR_acquired", P + "Bla_ESBL_inhR_acquired"]

def present(v: str) -> bool:
    return bool(v) and v.strip() not in ("", "-", "0", "nan", "NA", "None")

def fisher(a, b, c, d):
    from scipy.stats import fisher_exact
    return fisher_exact([[a, b], [c, d]])

def istrue(v): return str(v) in ("True", "1", "true")

rec = list(csv.DictReader(open(REC), delimiter="\t"))
# carbapenemase-negative, binary primary set
cn = [r for r in rec if not istrue(r["carb_pos"]) and r["call"] in ("S", "R")]
ids = {r["genome_id"]: r for r in cn}

# join Kleborate determinant columns
kb = {}
with KB.open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        kb[r["strain"]] = r

def flag_from_kb(gid, cols):
    r = kb.get(gid)
    if r is None:
        return None
    return any(present(r.get(c, "")) for c in cols)

# build per-isolate determinant flags
det = {"OmpK36 disruption": {}, "Acquired ESBL": {},
       "Chromosomal beta-lactamase": {}, "SHV mutations": {},
       "Inhibitor-resistant beta-lactamase": {}}
resist = {}
missing_kb = 0
for gid, r in ids.items():
    resist[gid] = (r["call"] == "R")
    det["OmpK36 disruption"][gid] = istrue(r["porin_disrupted"])
    det["Acquired ESBL"][gid] = istrue(r["esbl"])
    chr_ = flag_from_kb(gid, [COL_CHR])
    shv = flag_from_kb(gid, [COL_SHV])
    inhr = flag_from_kb(gid, COL_INHR)
    if chr_ is None:
        missing_kb += 1
    det["Chromosomal beta-lactamase"][gid] = bool(chr_)
    det["SHV mutations"][gid] = bool(shv)
    det["Inhibitor-resistant beta-lactamase"][gid] = bool(inhr)

rows = [["determinant", "present_n", "resistant_pct_present", "OR", "p"]]
for name, flags in det.items():
    gids = list(ids)
    a = sum(1 for g in gids if flags[g] and resist[g])          # present & R
    b = sum(1 for g in gids if flags[g] and not resist[g])      # present & S
    c = sum(1 for g in gids if not flags[g] and resist[g])      # absent & R
    d = sum(1 for g in gids if not flags[g] and not resist[g])  # absent & S
    n_present = a + b
    pct = 100 * a / n_present if n_present else 0.0
    orr, p = fisher(a, b, c, d)
    rows.append([name, n_present, f"{pct:.1f}", f"{orr:.2f}", f"{p:.2g}"])

with (TAB / "table4_determinants.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerows(rows)
    fh.write(f"# carbapenemase-negative binary isolates: {len(ids)}\n")
    fh.write(f"# genomes without a Kleborate record (chr/SHV/inhR treated absent): {missing_kb}\n")

# ESBL x OmpK36 interaction on resistance, via logistic regression LR test
try:
    import numpy as np, statsmodels.api as sm
    gids = list(ids)
    y = np.array([1 if resist[g] else 0 for g in gids])
    e = np.array([1 if det["Acquired ESBL"][g] else 0 for g in gids])
    o = np.array([1 if det["OmpK36 disruption"][g] else 0 for g in gids])
    X_full = sm.add_constant(np.column_stack([e, o, e * o]))
    X_red = sm.add_constant(np.column_stack([e, o]))
    m_full = sm.Logit(y, X_full).fit(disp=0)
    m_red = sm.Logit(y, X_red).fit(disp=0)
    lr = 2 * (m_full.llf - m_red.llf)
    from scipy.stats import chi2
    p_int = chi2.sf(lr, 1)
    method = "logistic interaction LR test (1 df)"
except Exception as ex:
    # fall back to Breslow-Day-free: compare 2x2 within strata is not a single stat;
    # report the simple independence note
    lr = float("nan"); p_int = float("nan"); method = f"unavailable ({ex})"

with (TAB / "table4_esbl_ompk36_interaction.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["test", "statistic", "df", "p"])
    w.writerow([method, f"{lr:.2f}", 1, f"{p_int:.3f}"])

print("wrote table4_determinants.tsv and table4_esbl_ompk36_interaction.tsv")
for r in rows: print("  ", r)
print(f"interaction: {method} chi2={lr:.2f} p={p_int:.3f}")
