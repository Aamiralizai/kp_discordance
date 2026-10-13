#!/usr/bin/env python3
"""22_compare_truncation_fix.py — every headline estimate computed from genotypes
called on the incomplete assembly files and from genotypes called on the verified
assemblies (same phenotypes, same QC-retained genomes; only the assembly from
which 3,443 genotypes were called differs).

  python scripts/22_compare_truncation_fix.py results_pre_redownload results
Writes results/supplementary/S17_truncation_fix_comparison.tsv
"""
import json, sys, csv
from pathlib import Path

before_dir, after_dir = Path(sys.argv[1]), Path(sys.argv[2])
B = {a["analysis"]: a for a in json.load(open(before_dir / "manuscript_results.json"))["analyses"]}
A = {a["analysis"]: a for a in json.load(open(after_dir / "manuscript_results.json"))["analyses"]}

EN = "\u2013"

def ci(x):
    return f"{x[0]:.1f} ({x[1]:.1f}{EN}{x[2]:.1f})" if x else ""

rows = [["analysis", "quantity", "incomplete_assemblies", "verified_assemblies"]]
for name in A:
    a, b = A[name], B.get(name, {})
    def add(q, fb, fa):
        rows.append([name, q, fb, fa])
    add("binary set n", b.get("binary_set", ""), a.get("binary_set", ""))
    add("carbapenemase-positive n", b.get("carbapenemase_positive_n", ""), a.get("carbapenemase_positive_n", ""))
    add("carbapenemase-positive susceptible % (95% CI)",
        ci(b.get("carbapenemase_positive_susceptible_pct_ci")), ci(a.get("carbapenemase_positive_susceptible_pct_ci")))
    add("carbapenemase-negative n", b.get("carbapenemase_negative_n", ""), a.get("carbapenemase_negative_n", ""))
    add("carbapenemase-negative resistant % (95% CI)",
        ci(b.get("carbapenemase_negative_resistant_pct_ci")), ci(a.get("carbapenemase_negative_resistant_pct_ci")))
    for k in ("PPV", "NPV", "sensitivity", "specificity"):
        fb, fa = b.get(k, ""), a.get(k, "")
        add(k, f"{fb:.3f}" if isinstance(fb, float) else fb,
            f"{fa:.3f}" if isinstance(fa, float) else fa)
    for fam in sorted(set(a.get("by_family", {})) | set(b.get("by_family", {}))):
        fb, fa = b.get("by_family", {}).get(fam, {}), a.get("by_family", {}).get(fam, {})
        add(f"{fam}: susceptible % (95% CI), n",
            f"{ci(fb.get('pct_ci'))} n={fb.get('n','')}" if fb else "",
            f"{ci(fa.get('pct_ci'))} n={fa.get('n','')}" if fa else "")
    fb = {r["family"]: r for r in b.get("fisher_intact_vs_ompK36", [])}
    for r in a.get("fisher_intact_vs_ompK36", []):
        o = fb.get(r["family"], {})
        def fmt(d):
            if not d:
                return ""
            orr = ("not estimable" if d["odds_ratio"] != d["odds_ratio"]
                   else f"OR {d['odds_ratio']:.3f}, p = {d['p_value']:.1e}")
            return (f"{orr} (intact S/R {d['S_intact']}/{d['R_intact']}, "
                    f"disrupted S/R {d['S_disrupted']}/{d['R_disrupted']})")
        add(f"{r['family']}: porin contrast, intact vs OmpK36-disrupted", fmt(o), fmt(r))

out = after_dir / "supplementary" / "S17_truncation_fix_comparison.tsv"
with open(out, "w", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(rows)
w = max(len(r[1]) for r in rows)
cur = None
for r in rows[1:]:
    if r[0] != cur:
        cur = r[0]; print(f"\n== {cur}")
    print(f"  {r[1]:<{w}}  {str(r[2]):<45} -> {r[3]}")
print(f"\nwrote {out}")
