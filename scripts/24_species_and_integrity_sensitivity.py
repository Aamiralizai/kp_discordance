#!/usr/bin/env python3
"""24_species_and_integrity_sensitivity.py — two sensitivity analyses that
answer questions a reader will reasonably ask of this dataset.

1. SPECIES. Genomes are retained when Kleborate assigns a species whose name
   contains "pneumoniae", which keeps K. quasipneumoniae alongside
   K. pneumoniae. The title and the claims are about K. pneumoniae, so the
   primary analysis is repeated on K. pneumoniae sensu stricto.

2. ASSEMBLY INTEGRITY. Every genome in the primary set is cross-referenced
   against the effective integrity verdict (00g), so that the analysis set can
   be stated to rest on assemblies verified against the deposited record rather
   than assumed to.

Writes results/supplementary/S18_species_integrity_sensitivity.tsv
"""
import csv
import math
from pathlib import Path

REC = Path("results/tables/recordlevel_discordance_meropenem_clsi_bmd.tsv")
KLEB = Path("data/processed/kleborate_klebsiella_pneumo_complex_output.prefilter.tsv")
EFF = Path("results/qc_concordance/fasta_integrity_effective.tsv")
OUT = Path("results/supplementary/S18_species_integrity_sensitivity.tsv")
SPECIES_COL = "enterobacterales__species__species"


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return 100 * p, max((c - m) / d * 100, 0.0), (c + m) / d * 100


def rows(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith("#")],
                                   delimiter="\t"))


rec = [r for r in rows(REC) if r["call"] in ("S", "R")]
species = {r["strain"].strip(): r[SPECIES_COL] for r in rows(KLEB)}
verdict = {r["genome_id"]: r["verdict"] for r in rows(EFF)} if EFF.exists() else {}


def block(label, subset):
    pos = [r for r in subset if r["carb_pos"].lower() in ("true", "1")]
    neg = [r for r in subset if r["carb_pos"].lower() not in ("true", "1")]
    ps = sum(1 for r in pos if r["call"] == "S")
    nr = sum(1 for r in neg if r["call"] == "R")
    tp, tn = len(pos) - ps, len(neg) - nr
    p1, lo1, hi1 = wilson(ps, len(pos))
    p2, lo2, hi2 = wilson(nr, len(neg))
    return [label, len(subset), len(pos), ps, f"{p1:.1f}", f"{lo1:.1f}", f"{hi1:.1f}",
            len(neg), nr, f"{p2:.1f}", f"{lo2:.1f}", f"{hi2:.1f}",
            f"{tp / (tp + ps):.3f}", f"{(tp + tn) / len(subset):.3f}"]


out = [["analysis_set", "n_binary", "carbapenemase_positive_n",
        "carbapenemase_positive_susceptible_n", "pct", "ci_lo", "ci_hi",
        "carbapenemase_negative_n", "carbapenemase_negative_resistant_n",
        "pct_r", "ci_lo_r", "ci_hi_r", "PPV", "accuracy"]]
out.append(block("primary (K. pneumoniae species complex, as reported)", rec))
strict = [r for r in rec if species.get(r["genome_id"]) == "Klebsiella pneumoniae"]
out.append(block("K. pneumoniae sensu stricto", strict))

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(out)

# WHICH FAMILIES THE NON-pneumoniae GENOMES SIT IN.
# The headline estimates barely move when K. quasipneumoniae is dropped, which
# the table above shows. That is weaker than what is actually true, and the
# stronger statement is the one the title rests on, so it is computed here
# rather than asserted in the text. The family and porin findings rest on the
# OXA-48-like and VIM producers; if none of those is K. quasipneumoniae then
# that analysis is not merely robust to the species definition, it is the same
# isolates either way.
from collections import Counter

FOCUS = ("OXA-48-like", "VIM")
focus_rows = [r for r in rec
              if any(f in (r.get("carb_family") or "") for f in FOCUS)]
focus_quasi = [r for r in focus_rows
               if species.get(r["genome_id"]) != "Klebsiella pneumoniae"]
non_strict = [r for r in rec
              if species.get(r["genome_id"]) != "Klebsiella pneumoniae"]
comp = Counter(r.get("carb_family") or "none" for r in non_strict)
comp_path = OUT.parent / "S18b_non_pneumoniae_composition.tsv"
with open(comp_path, "w", newline="") as fh:
    wr = csv.writer(fh, delimiter="\t")
    wr.writerow(["carbapenemase_family", "n_non_pneumoniae"])
    for k, v in sorted(comp.items(), key=lambda kv: (-kv[1], kv[0])):
        wr.writerow([k, v])
    wr.writerow(["TOTAL", len(non_strict)])
print("\ncarbapenemase families of the non-pneumoniae genomes "
      f"({len(non_strict)} of {len(rec)}):")
for k, v in sorted(comp.items(), key=lambda kv: (-kv[1], kv[0])):
    print(f"   {v:>5}  {k}")
print(f"   isolates in a family containing OXA-48-like or VIM: "
      f"{len(focus_rows)}, of which non-pneumoniae {len(focus_quasi)}")
if focus_quasi:
    print("   NOTE: the family and porin analysis is NOT species-invariant; "
          "the manuscript may not say that it is.")
else:
    print("   every OXA-48-like and VIM producer is K. pneumoniae sensu "
          "stricto, so the family and porin analysis is identical under "
          "either species definition")
print(f"[out] {comp_path}")

sp = Counter(species.get(r["genome_id"], "unknown") for r in rec)
vd = Counter(verdict.get(r["genome_id"], "ABSENT") for r in rec)
print("species composition of the primary binary set:")
for k, v in sp.most_common():
    print(f"   {v:>5}  {k}")
print("\nassembly-integrity verdict of the primary binary set:")
for k, v in vd.most_common():
    print(f"   {v:>5}  {k}")
print(f"\n[out] {OUT}")
for r in out[1:]:
    print(f"   {r[0]:<52} n={r[1]:<5} carb+ {r[3]}/{r[2]} = {r[4]}%  "
          f"carb- {r[8]}/{r[7]} = {r[9]}%  acc {r[13]}")
