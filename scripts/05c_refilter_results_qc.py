#!/usr/bin/env python3
"""
05c_refilter_results_qc.py — apply assembly QC to already-generated results.

The supported route is to filter the genotype table with 05b_apply_qc_filter.py
and then run run_all.sh, which regenerates everything from the phenotype
records. Use this script only when the raw records are not to hand and the
record-level tables already exist, for example when reconciling a revision
against a previously deposited result set.

WHAT IT DOES

For every recordlevel_discordance_*.tsv it drops genomes whose assembly QC
verdict is FAIL, then rebuilds the four tables that 08b derives from it, using
the same logic and the same minimum-n rule:

    recordlevel_by_family_*     susceptibility by carbapenemase family
    recordlevel_mechanism_*     porin class within family
    recordlevel_fisher_*        intact vs OmpK36-disrupted, BH-adjusted
    recordlevel_audit_*         stage counts, which Figure 1 reads

WHAT IT CANNOT DO, AND WHY IT MATTERS

The audit's first six rows count phenotype RECORDS, before records are
collapsed to one row per genome. Those counts cannot be recomputed from a
per-genome table, so they are carried forward unchanged and flagged in the
output. They are slightly too high, because records belonging to a
QC-failed genome are still counted. The genome-level rows below them are
recomputed and are correct.

Run the full pipeline before submitting, so the deposited results are
internally consistent. This script exists to show what changes, not to replace
that.

USAGE
    python scripts/05c_refilter_results_qc.py
    python scripts/05c_refilter_results_qc.py --keep-warn false
"""
from __future__ import annotations

import argparse
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from kp_defs import PORIN_LABEL, PORIN_ORDER  # noqa: E402

TAB = Path("results/tables")
QC = Path("data/interim/assembly_qc_metrics.tsv")

# Rows above this marker count records, not genomes, and are carried forward.
RECORD_STAGES = 7


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = (z / d) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (100 * p, 100 * max(0.0, c - h), 100 * min(1.0, c + h))


def rebuild(d: pd.DataFrame, min_n: int, out_stem: str, n_removed: int,
            audit_src: Path | None) -> dict:
    """Rebuild the derived tables from a filtered per-genome table."""
    d = d.copy()
    d["carb_pos"] = d["carb_pos"].astype(str).str.lower().isin({"true", "1"})
    d["porin_disrupted"] = (d["porin_disrupted"].astype(str).str.lower()
                            .isin({"true", "1"}))
    binary = d[d["call"].isin(["S", "R"])]
    pos, neg = binary[binary.carb_pos], binary[~binary.carb_pos]

    # --- by family
    rows = []
    for fam, g in pos.groupby("carb_family"):
        if len(g) < min_n:
            continue
        k = int((g["call"] == "S").sum())
        p, lo, hi = wilson(k, len(g))
        rows.append(dict(family=fam, n=len(g), n_susceptible=k,
                         pct=round(p, 1), ci_lo=round(lo, 1), ci_hi=round(hi, 1)))
    fam_tab = (pd.DataFrame(rows).sort_values("pct", ascending=False)
               if rows else pd.DataFrame())

    # --- porin within family
    rows = []
    for fam in (fam_tab["family"] if len(fam_tab) else []):
        sub = pos[pos.carb_family == fam]
        for pr in PORIN_ORDER:
            g = sub[sub.porin == pr]
            if len(g) < min_n:
                continue
            k = int((g["call"] == "S").sum())
            p, lo, hi = wilson(k, len(g))
            rows.append(dict(family=fam, porin=pr, label=PORIN_LABEL[pr],
                             n=len(g), n_susc=k, pct=round(p, 1),
                             ci_lo=round(lo, 1), ci_hi=round(hi, 1)))
    mech = pd.DataFrame(rows)

    # --- Fisher, intact vs OmpK36-disrupted
    ftab = []
    for fam in (mech["family"].unique() if len(mech) else []):
        sub = pos[pos.carb_family == fam]
        a11 = int(((sub.porin == "none") & (sub["call"] == "R")).sum())
        a12 = int(((sub.porin == "none") & (sub["call"] == "S")).sum())
        a21 = int((sub.porin_disrupted & (sub["call"] == "R")).sum())
        a22 = int((sub.porin_disrupted & (sub["call"] == "S")).sum())
        if min(a11 + a12, a21 + a22) < min_n:
            continue
        orr, pv = stats.fisher_exact([[a11, a12], [a21, a22]])
        ftab.append(dict(family=fam, R_intact=a11, S_intact=a12,
                         R_disrupted=a21, S_disrupted=a22,
                         odds_ratio=orr, p_value=pv))
    ft = pd.DataFrame(ftab)
    if len(ft):
        pv = ft["p_value"].values
        order = np.argsort(pv)
        n_t = len(ft)
        q = np.minimum.accumulate(
            (np.sort(pv) * n_t / (np.arange(n_t) + 1))[::-1])[::-1]
        qq = np.empty(n_t)
        qq[order] = np.clip(q, 0, 1)
        ft["q_BH"] = qq

    # --- audit
    if audit_src and audit_src.exists():
        au = pd.read_csv(audit_src, sep="\t", dtype=str)
        carried = au.head(RECORD_STAGES).copy()
        kb = {r["stage"]: r["n"] for _, r in au.iterrows()}
        audit = [(s, carried.iloc[i]["n"]) for i, s in
                 enumerate(carried["stage"])]
        for s in ("Kleborate genomes", "confirmed K. pneumoniae",
                  "genotyping failures"):
            if s in kb:
                audit.append((s, kb[s]))
        audit.append(("failed assembly quality control (QualiBact v1.0)",
                      n_removed))
    else:
        audit = [("failed assembly quality control (QualiBact v1.0)", n_removed)]
    audit += [
        ("genomes with linked genotype", len(d)),
        ("susceptible", int((d["call"] == "S").sum())),
        ("intermediate", int((d["call"] == "I").sum())),
        ("resistant", int((d["call"] == "R").sum())),
        ("unresolvable", int((d["call"] == "ambiguous").sum())),
        ("binary analysis set", int(d["call"].isin(["S", "R"]).sum())),
    ]

    # --- write
    d.to_csv(TAB / f"recordlevel_discordance_{out_stem}.tsv", sep="\t", index=False)
    pd.DataFrame(audit, columns=["stage", "n"]).to_csv(
        TAB / f"recordlevel_audit_{out_stem}.tsv", sep="\t", index=False)
    if len(fam_tab):
        fam_tab.to_csv(TAB / f"recordlevel_by_family_{out_stem}.tsv",
                       sep="\t", index=False)
    if len(mech):
        mech.to_csv(TAB / f"recordlevel_mechanism_{out_stem}.tsv",
                    sep="\t", index=False)
    if len(ft):
        ft.to_csv(TAB / f"recordlevel_fisher_{out_stem}.tsv",
                  sep="\t", index=False)

    fs = int((pos["call"] == "S").sum())
    fr = int((neg["call"] == "R").sum())
    return dict(linked=len(d), binary=len(binary),
                pos=len(pos), pos_susc=fs, neg=len(neg), neg_res=fr)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc", type=Path, default=QC)
    ap.add_argument("--tables", type=Path, default=TAB)
    ap.add_argument("--min-n", type=int, default=10)
    ap.add_argument("--keep-warn", default="true",
                    help="keep WARN genomes (true, the default) or not (false)")
    ap.add_argument("--backup-dir", type=Path,
                    default=Path("results/tables_prefilter"))
    a = ap.parse_args()

    if not a.qc.exists():
        sys.exit(f"Missing {a.qc}. Run 00e_fetch_bvbrc_qc.py or 00b_qc_assembly.py.")
    qc = pd.read_csv(a.qc, sep="\t", dtype=str, low_memory=False)
    keep_set = {"PASS"} | ({"WARN"} if a.keep_warn.lower() == "true" else set())
    verdict = dict(zip(qc["sample"].astype(str), qc["verdict"].astype(str)))
    print(f"QC results : {len(verdict):,} genomes")
    print(f"keeping    : {', '.join(sorted(keep_set))}")

    files = sorted(a.tables.glob("recordlevel_discordance_*.tsv"))
    if not files:
        sys.exit(f"No recordlevel_discordance_*.tsv under {a.tables}")

    # one-time backup so the unfiltered set is never lost
    if not a.backup_dir.exists():
        shutil.copytree(a.tables, a.backup_dir)
        print(f"backed up unfiltered tables -> {a.backup_dir}")

    print()
    for f in files:
        stem = f.name[len("recordlevel_discordance_"):-len(".tsv")]
        src = a.backup_dir / f.name
        d = pd.read_csv(src if src.exists() else f, sep="\t",
                        dtype=str, low_memory=False)
        n0 = len(d)
        d["_v"] = d["genome_id"].astype(str).map(verdict).fillna("no_qc_result")
        kept = d[d["_v"].isin(keep_set)].drop(columns=["_v"])
        removed = n0 - len(kept)
        audit_src = a.backup_dir / f"recordlevel_audit_{stem}.tsv"
        r = rebuild(kept, a.min_n, stem, removed,
                    audit_src if audit_src.exists() else None)
        ps = wilson(r["pos_susc"], r["pos"])
        nr = wilson(r["neg_res"], r["neg"])
        print(f"{stem:28s} {n0:5,} -> {r['linked']:5,}  (-{removed:3d})   "
              f"gene+ susceptible {ps[0]:5.1f}% [{ps[1]:4.1f}-{ps[2]:4.1f}]   "
              f"gene- resistant {nr[0]:5.1f}% [{nr[1]:4.1f}-{nr[2]:4.1f}]")

    print(f"\nThe first {RECORD_STAGES} audit rows count RECORDS and are carried")
    print("forward unchanged. They are slightly too high, because records")
    print("belonging to a QC-failed genome are still counted. Regenerate them")
    print("by running the full pipeline before submission.")
    print("\nNext: scripts 09 to 16, then 17 and 18.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
