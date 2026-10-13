#!/usr/bin/env python3
"""00g_effective_integrity_report.py — integrity verdict for the assembly that
was actually genotyped, built from measurements already on record.

WHY THIS EXISTS
00f_check_download_integrity.py re-reads every FASTA, which is the right thing
to do once but is wasteful afterwards. Both measurements already exist:

  results/qc_concordance/fasta_integrity_report.tsv  the original download audit
  data/raw/genomes_redownload/manifest.tsv           the verified re-retrieval

This script combines them, taking the re-retrieved copy wherever one exists
because that is the file the genotyper read, and applies one rule to both:

  COMPLETE    length == deposited length AND contigs == deposited contigs
              (where a deposited contig count exists) AND the file ends cleanly
  INCOMPLETE  any of those fails
  UNVERIFIED  no deposited record to compare against

The original audit classified on the missing final newline alone. That let
through 40 files which end cleanly but are short against the deposited record,
so the verdict column is recomputed here rather than carried over.

    python3 scripts/00g_effective_integrity_report.py
"""
import csv
from pathlib import Path

ORIG = Path("results/qc_concordance/fasta_integrity_report.tsv")
REDL = Path("data/raw/genomes_redownload/manifest.tsv")
OUT = Path("results/qc_concordance/fasta_integrity_effective.tsv")


def num(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def verdict(size, contigs, newline, exp_size, exp_ctg):
    if exp_size is None:
        return "UNVERIFIED"
    if size != exp_size:
        return "INCOMPLETE"
    if exp_ctg and contigs is not None and contigs != exp_ctg:
        return "INCOMPLETE"
    if newline is not None and not newline:
        return "INCOMPLETE"
    return "COMPLETE"


rows = {}
for r in csv.DictReader(open(ORIG), delimiter="\t"):
    rows[r["genome_id"]] = dict(
        genome_id=r["genome_id"], source="data/raw/genomes",
        local_size_bp=num(r["local_size_bp"]), local_contigs=num(r["local_contigs"]),
        ends_newline=num(r["ends_newline"]),
        bvbrc_size_bp=num(r["bvbrc_size_bp"]), bvbrc_contigs=num(r["bvbrc_contigs"]))

n_replaced = 0
if REDL.exists():
    for r in csv.DictReader(open(REDL), delimiter="\t"):
        if r.get("status") not in (None, "", "ok") and r.get("verified") != "yes":
            continue
        gid = r["genome_id"]
        got = num(r.get("got_length") or r.get("file_length"))
        gotc = num(r.get("got_contigs") or r.get("file_contigs"))
        if got is None:
            continue
        n_replaced += gid in rows
        rows[gid] = dict(
            genome_id=gid, source="data/raw/genomes_redownload",
            local_size_bp=got, local_contigs=gotc,
            ends_newline=num(r.get("ends_newline")) if r.get("ends_newline") else 1,
            bvbrc_size_bp=num(r["bvbrc_length"]), bvbrc_contigs=num(r["bvbrc_contigs"]))

out = []
for gid in sorted(rows):
    r = rows[gid]
    v = verdict(r["local_size_bp"], r["local_contigs"], r["ends_newline"],
                r["bvbrc_size_bp"], r["bvbrc_contigs"])
    ratio = (f"{r['local_size_bp'] / r['bvbrc_size_bp']:.4f}"
             if r["bvbrc_size_bp"] else "")
    out.append([gid, r["source"], r["local_size_bp"], r["local_contigs"],
                r["ends_newline"], r["bvbrc_size_bp"] or "",
                r["bvbrc_contigs"] or "", ratio, v])

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genome_id", "source_dir", "local_size_bp", "local_contigs",
                "ends_newline", "bvbrc_size_bp", "bvbrc_contigs", "size_ratio",
                "verdict"])
    w.writerows(out)

from collections import Counter
c = Counter(r[-1] for r in out)
print(f"  genomes            : {len(out):,}")
print(f"  re-retrieved copies: {n_replaced:,}")
for k in ("COMPLETE", "INCOMPLETE", "UNVERIFIED"):
    print(f"  {k:<11}: {c.get(k, 0):,}")
print(f"\n  [out] {OUT}")
