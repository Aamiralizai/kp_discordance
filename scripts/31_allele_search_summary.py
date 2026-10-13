#!/usr/bin/env python3
"""Rebuild the uncatalogued-allele search summary from the three run outputs.

The summary is derived here rather than written by hand so that it cannot
drift from hits.jsonl, short_orfs.jsonl and fragment_loci.jsonl.

Two screening bands were run. The main band covers 180 to 500 residues and
answers whether an uncatalogued carbapenemase is present. The short band
covers 60 to 180 residues and answers whether a catalogued carbapenemase is
present but truncated. The first version of this analysis ran only the main
band, which excluded every truncated gene and so could not have detected the
one mechanism the discordance hypothesis predicts.
"""
import collections
import csv
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from allele_classify import CATALOGUED, NOVEL, classify  # noqa: E402

# the scripts live in scripts/ and their data in results/allele_search/, so
# both are resolved from the repository root rather than from the file's own
# directory, which lets the script be moved without silently writing elsewhere
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "results", "allele_search")
SUPP = os.path.join(ROOT, "results", "supplementary")  # where every other supplementary table lives
OUT = os.path.join(SUPP, "S19a_allele_search_summary.tsv")
COHORTS = ["A_carbNeg_resistant", "B_carbPos_susceptible", "D_control_sample"]


def jsonl(name):
    path = os.path.join(WORK, name)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def carb_hits(rec):
    """Only carbapenemase-subclass references bear on the question.

    The main-band search scored every beta-lactamase-like protein against the
    whole reference set, so a record also contains exact matches to the native
    SHV and AmpC enzymes. Counting those as carbapenemases would report a
    carbapenemase in almost every genome, including the controls.
    """
    out = []
    for h in rec["hits"]:
        if h.get("subclass") != "CARBAPENEM":
            continue
        h = dict(h)
        # reclassify rather than trust the verdict stored by the run, so a
        # correction to the classifier reaches the summary without a rerun
        h["verdict"] = classify(h["identity"], h["ref_coverage"], h["subclass"])
        out.append(h)
    return out


def main():
    os.makedirs(SUPP, exist_ok=True)
    main_band = {r["genome_id"]: r for r in jsonl("hits.jsonl")}
    short_band = {r["genome_id"]: r for r in jsonl("short_orfs.jsonl")}
    loci = jsonl("fragment_loci.jsonl")

    # one row per genome that carries at least one carbapenemase gene fragment,
    # classified by where the fragment sits on its contig
    verdict = collections.defaultdict(set)
    for r in loci:
        if not r.get("ref"):
            continue
        if r["identity"] < 0.90:
            kind = "divergent_chromosomal_homologue"
        elif r["assembly_break"]:
            kind = "assembly_break"
        else:
            kind = "internal_truncation"
        verdict[r["genome"]].add(kind)

    rows = []
    for cohort in COHORTS:
        gids = [g for g, r in main_band.items() if r["cohort"] == cohort]
        exact = sum(1 for g in gids
                    if any(h["verdict"] in CATALOGUED for h in carb_hits(main_band[g])))
        novel = sum(1 for g in gids
                    if any(h["verdict"] in NOVEL for h in carb_hits(main_band[g])))
        # the background ceiling is the best non-carbapenemase-call hit each
        # genome achieves, not every hit pooled, since one genome contributing
        # hundreds of weak alignments would otherwise dominate the median
        bg = []
        for g in gids:
            rest = [h["identity"] for h in carb_hits(main_band[g])
                    if h["verdict"] not in CATALOGUED + NOVEL]
            if rest:
                bg.append(max(rest))
        short_orfs = sum(short_band[g]["short_orfs_screened"]
                         for g in gids if g in short_band)
        breaks = sum(1 for g in gids if "assembly_break" in verdict.get(g, ()))
        trunc = sum(1 for g in gids if "internal_truncation" in verdict.get(g, ()))
        rows.append(dict(
            cohort=cohort,
            genomes=len(gids),
            proteins_screened_main_band=sum(main_band[g]["screened"] for g in gids),
            short_orfs_screened=short_orfs,
            genomes_with_catalogued_carbapenemase_main_band=exact,
            genomes_with_novel_candidate_carbapenemase=novel,
            # the locus columns below draw on both screening bands
            genomes_with_carbapenemase_split_by_assembly_break=breaks,
            genomes_with_internally_truncated_carbapenemase=trunc,
            background_ceiling_median_pct=round(100 * statistics.median(bg), 1) if bg else "",
            background_ceiling_max_pct=round(100 * max(bg), 1) if bg else "",
        ))

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT)
    for r in rows:
        print(" ", r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
