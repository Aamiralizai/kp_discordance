#!/usr/bin/env python3
"""Regenerate the allele-search supplementary tables from the run outputs.

The first version of this table was produced by an ad hoc command rather than
a script, which left a supplementary table with no generator behind it. It is
rebuilt here so that a correction to the classifier reaches the table.

One row per genome in S19c and one row per gene fragment in S19b. The best
carbapenemase-subclass hit is reported with the
verdict the shared classifier assigns, and the short-band result is carried
alongside it so both screening windows appear in one place.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from allele_classify import classify  # noqa: E402

# the scripts live in scripts/ and their data in results/allele_search/, so
# both are resolved from the repository root rather than from the file's own
# directory, which lets the script be moved without silently writing elsewhere
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "results", "allele_search")
SUPP = os.path.join(ROOT, "results", "supplementary")  # where every other supplementary table lives
OUT = os.path.join(SUPP, "S19c_allele_search_per_genome.tsv")
OUT_LOCI = os.path.join(SUPP, "S19b_carbapenemase_gene_fragments.tsv")
COHORT_ORDER = {"A_carbNeg_resistant": 0, "B_carbPos_susceptible": 1,
                "D_control_sample": 2}


def jsonl(name):
    path = os.path.join(WORK, name)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def locus_kind(r):
    if r["identity"] < 0.90:
        return "divergent_chromosomal_homologue"
    return "assembly_break" if r["assembly_break"] else "internal_truncation"


def pick_reference(rows_for_genome):
    """Choose one reference allele per genome.

    Several KPC variants differ only outside the region a fragment covers, so
    one fragment can match a dozen of them at the same identity. Where two or
    more open reading frames of a genome match the same reference, the variant
    whose own length the fragments account for exactly is the one they
    reconstruct, so that reference is named. Otherwise the best single
    alignment is used.
    """
    by_ref = {}
    for r in rows_for_genome:
        by_ref.setdefault(r["ref"], []).append(r)
    multi = {k: v for k, v in by_ref.items()
             if len({(x["contig_len"], x["orf_start"]) for x in v}) > 1}
    pool = multi or by_ref

    def score(item):
        ref, rs = item
        aligned = sum(x["aligned_aa"] for x in rs)
        reflen = rs[0]["ref_len"]
        mean_id = sum(x["identity"] for x in rs) / len(rs)
        return (abs(aligned - reflen), -mean_id, ref)

    return min(pool.items(), key=score)[0]


def one_row_per_locus(rows):
    """One row per open reading frame, under one reference allele per genome."""
    per_genome = {}
    for r in rows:
        per_genome.setdefault(r["genome"], []).append(r)
    out = {}
    for gid, grows in per_genome.items():
        ref = pick_reference(grows)
        for r in grows:
            if r["ref"] != ref:
                continue
            key = (gid, r["contig_len"], r["orf_start"])
            prev = out.get(key)
            if prev is None or r["aligned_aa"] > prev["aligned_aa"]:
                out[key] = r
        # an open reading frame the chosen reference does not cover keeps its own
        for r in grows:
            key = (gid, r["contig_len"], r["orf_start"])
            if key not in out:
                prev = out.get(key)
                if prev is None or (r["identity"], r["aligned_aa"]) > \
                        (prev["identity"], prev["aligned_aa"]):
                    out[key] = r
    return list(out.values())


def write_loci(calls, cohort_of):
    rows = one_row_per_locus([r for r in jsonl("fragment_loci.jsonl") if r.get("ref")])
    if not rows:
        return
    rows.sort(key=lambda r: (COHORT_ORDER.get(cohort_of.get(r["genome"], ""), 9),
                             r["genome"],
                             int(r["ref_residues"].split("-")[0])))
    with open(OUT_LOCI, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cohort", "genome_id", "published_carbapenemase_call",
                    "reference_allele", "percent_identity", "aligned_residues",
                    "reference_residue_range", "reference_length",
                    "contig_length_bp", "orf_start", "orf_end", "strand",
                    "partial_at_orf_start", "partial_at_orf_end",
                    "bases_to_contig_start", "bases_to_contig_end",
                    "interpretation"])
        for r in rows:
            call = calls.get(r["genome"], ["NA", "NA"])
            w.writerow([cohort_of.get(r["genome"], "NA"), r["genome"],
                        call[1] if len(call) > 1 else call[0], r["ref"],
                        f"{100 * r['identity']:.1f}", r["aligned_aa"],
                        r["ref_residues"], r["ref_len"], r["contig_len"],
                        r["orf_start"], r["orf_end"],
                        "+" if r["strand"] > 0 else "-",
                        "yes" if r["partial_begin"] else "no",
                        "yes" if r["partial_end"] else "no",
                        r["bases_to_contig_start"], r["bases_to_contig_end"],
                        locus_kind(r)])
    print("wrote", OUT_LOCI, len(rows), "loci")


def main():
    os.makedirs(SUPP, exist_ok=True)
    calls = json.load(open(os.path.join(WORK, "carb_calls.json")))
    short = {r["genome_id"]: r for r in jsonl("short_orfs.jsonl")}
    loci = {}
    for r in jsonl("fragment_loci.jsonl"):
        if r.get("ref"):
            loci.setdefault(r["genome"], set()).add(locus_kind(r))

    rows = []
    for rec in jsonl("hits.jsonl"):
        gid = rec["genome_id"]
        carb = [h for h in rec["hits"] if h.get("subclass") == "CARBAPENEM"]
        for h in carb:
            h["verdict"] = classify(h["identity"], h["ref_coverage"], h["subclass"])
        best = max(carb, key=lambda h: (h["identity"], h["ref_coverage"])) if carb else None
        call = calls.get(gid, ["NA", "NA"])
        sh = short.get(gid, {})
        sh_best = max((h["identity"] for h in sh.get("carbapenemase_like", [])),
                      default=None)
        rows.append(dict(
            cohort=rec["cohort"],
            genome_id=gid,
            kleborate_family=call[0] if call else "NA",
            kleborate_allele=call[1] if len(call) > 1 else "NA",
            contigs=rec["contigs"],
            proteins_predicted=rec["proteins"],
            proteins_screened_main_band=rec["screened"],
            short_orfs_screened=sh.get("short_orfs_screened", ""),
            best_carbapenemase_reference=best["ref_symbol"] if best else "",
            best_identity_pct=round(100 * best["identity"], 1) if best else "",
            best_reference_coverage_pct=round(100 * best["ref_coverage"], 1) if best else "",
            protein_aa=best["aa_length"] if best else "",
            reference_aa=best["ref_aa_length"] if best else "",
            verdict_main_band=best["verdict"] if best else "no_carbapenemase_like_protein",
            best_short_orf_identity_pct=round(100 * sh_best, 1) if sh_best else "",
            fragment_locus_verdict=";".join(sorted(loci.get(gid, ()))) or "",
        ))

    rows.sort(key=lambda r: (COHORT_ORDER.get(r["cohort"], 9), r["genome_id"]))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT, len(rows), "rows")
    write_loci(calls, {r["genome_id"]: r["cohort"] for r in jsonl("hits.jsonl")})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
