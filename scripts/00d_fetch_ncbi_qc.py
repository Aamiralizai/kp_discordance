#!/usr/bin/env python3
"""
00d_fetch_ncbi_qc.py — assembly QC metrics from NCBI, without the assemblies.

An alternative to computing QC locally with 00b_qc_assembly.py, for the case
where the assemblies were obtained from BV-BRC, which mirrors GenBank. NCBI
already publishes, per assembly, the full QualiBact metric set.

    total_sequence_length                       genome size
    gc_count / atgc_count                       GC content, full precision
    contig_n50, number_of_contigs               contiguity
    checkm_info.completeness / .contamination   CheckM v1.2.x
    annotation_info.stats.gene_counts           coding sequences

This matters for two reasons. It removes the CheckM2 run, which is the slow
step, and the completeness and contamination values come from CheckM, which is
the tool QualiBact v1.0 used when deriving its thresholds, so the comparison is
like for like.

DECLARE THIS IN THE METHODS. These statistics describe the assembly as
deposited at NCBI, not the local copy that was genotyped. For a BV-BRC mirror
of a GenBank record the two are the same sequence, but that is an assumption
and it should be stated rather than assumed silently. 00b_qc_assembly.py
computes the first five metrics directly from the local FASTA and is the way to
verify it on a sample.

Output matches 00b_qc_assembly.py so that 05b_apply_qc_filter.py and
18_discordance_investigation.py consume either interchangeably.

USAGE
    python scripts/00d_fetch_ncbi_qc.py
    python scripts/00d_fetch_ncbi_qc.py --accessions results/supplementary/S1_derivation_accessions.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/dataset_report"
UA = ("kp_discordance/1.0 (research pipeline; "
      "+https://github.com/Aamiralizai/kp_discordance)")

# QualiBact v1.0, Klebsiella pneumoniae. Mirrors config/qc_thresholds.csv; that
# file remains the source of truth and is what 00b_qc_assembly.py reads.
THRESHOLDS = {
    "genome_size_bp":    dict(fail_lo=4_900_000, warn_lo=5_200_000,
                              warn_hi=6_000_000, fail_hi=6_400_000),
    "gc_content_pct":    dict(fail_lo=55.89, warn_lo=56.60,
                              warn_hi=57.52, fail_hi=58.02),
    "n50_bp":            dict(fail_lo=58_000, warn_lo=97_000,
                              warn_hi=None, fail_hi=None),
    "n_contigs":         dict(fail_lo=None, warn_lo=None,
                              warn_hi=270, fail_hi=350),
    "completeness_pct":  dict(fail_lo=95.0, warn_lo=100.0,
                              warn_hi=None, fail_hi=None),
    "contamination_pct": dict(fail_lo=None, warn_lo=None,
                              warn_hi=2.0, fail_hi=5.0),
    "n_cds":             dict(fail_lo=4800, warn_lo=4800,
                              warn_hi=5700, fail_hi=6500),
}
ORDER = list(THRESHOLDS)


def classify(v, t) -> str:
    if v is None:
        return "not_evaluated"
    if t["fail_lo"] is not None and v < t["fail_lo"]:
        return "FAIL"
    if t["fail_hi"] is not None and v > t["fail_hi"]:
        return "FAIL"
    if t["warn_lo"] is not None and v < t["warn_lo"]:
        return "WARN"
    if t["warn_hi"] is not None and v > t["warn_hi"]:
        return "WARN"
    return "PASS"


def post(ids: list[str], by: str = "accessions",
         retries: int = 4, timeout: int = 120) -> list[dict]:
    """Query by assembly accession (by='accessions') or BioSample
    (by='biosample_ids'). BV-BRC records an assembly accession for only about
    half of these genomes, so BioSample is the fallback route, not an extra."""
    body = json.dumps({by: ids,
                       "returned_content": "COMPLETE",
                       "page_size": len(ids) * 2 + 10}).encode()
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                API, data=body,
                headers={"Content-Type": "application/json",
                         "Accept": "application/json", "User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return json.loads(fh.read().decode()).get("reports", [])
        except (urllib.error.URLError, urllib.error.HTTPError,
                TimeoutError, json.JSONDecodeError) as exc:
            if attempt == retries:
                print(f"\n  batch failed after {retries} attempts: {exc}",
                      file=sys.stderr)
                return []
            time.sleep(2 * attempt)
    return []


def extract(rep: dict) -> dict:
    st = rep.get("assembly_stats", {}) or {}
    cm = rep.get("checkm_info", {}) or {}
    ann = ((rep.get("annotation_info", {}) or {}).get("stats", {}) or {})
    genes = (ann.get("gene_counts", {}) or {})

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    size = num(st.get("total_sequence_length"))
    gc_n, at_n = num(st.get("gc_count")), num(st.get("atgc_count"))
    # Recompute GC at full precision. gc_percent is rounded to one decimal,
    # and the QualiBact bounds are specified to two.
    gc = (100.0 * gc_n / at_n) if (gc_n and at_n) else num(st.get("gc_percent"))

    return {
        "genome_size_bp": size,
        "gc_content_pct": round(gc, 4) if gc is not None else None,
        "n50_bp": num(st.get("contig_n50")),
        "n_contigs": num(st.get("number_of_contigs")),
        "completeness_pct": num(cm.get("completeness")),
        "contamination_pct": num(cm.get("contamination")),
        "n_cds": num(genes.get("protein_coding")),
        "checkm_version": cm.get("checkm_version", ""),
        "assembly_level": (rep.get("assembly_info", {}) or {}).get("assembly_level", ""),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", type=Path,
                    default=Path("results/supplementary/S1_derivation_accessions.tsv"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/interim/assembly_qc_metrics.tsv"))
    ap.add_argument("--batch", type=int, default=200)
    ap.add_argument("--sleep", type=float, default=0.25)
    a = ap.parse_args()

    if not a.accessions.exists():
        sys.exit(f"Missing {a.accessions}. Run 19a_fetch_derivation_accessions.py.")
    acc = pd.read_csv(a.accessions, sep="\t", dtype=str)
    # genome_id is what the genotype table's `strain` column holds, so the
    # output must be keyed on it for 05b_apply_qc_filter.py to join.
    gid_acc, gid_bio = {}, {}
    for _, r in acc.iterrows():
        g = str(r["genome_id"])
        if pd.notna(r.get("assembly_accession")):
            gid_acc[g] = str(r["assembly_accession"])
        if pd.notna(r.get("biosample_accession")):
            gid_bio[g] = str(r["biosample_accession"])
    gids = [str(g) for g in acc["genome_id"]]
    print(f"genomes            : {len(gids):,}")
    print(f"  with assembly acc: {len(gid_acc):,}")
    print(f"  with biosample   : {len(gid_bio):,}")

    def better(a_rep, b_rep):
        """Prefer the record that actually carries CheckM, then RefSeq."""
        if a_rep is None:
            return b_rep
        ac, bc = a_rep.get("checkm_info"), b_rep.get("checkm_info")
        if bool(ac) != bool(bc):
            return a_rep if ac else b_rep
        aa, ba = a_rep.get("accession", ""), b_rep.get("accession", "")
        if aa.startswith("GCF_") != ba.startswith("GCF_"):
            return a_rep if aa.startswith("GCF_") else b_rep
        return a_rep

    by_gid: dict[str, dict] = {}

    # pass 1 — assembly accession
    acc_to_gid = {v: k for k, v in gid_acc.items()}
    accs = list(acc_to_gid)
    for i in range(0, len(accs), a.batch):
        for rep in post(accs[i:i + a.batch], "accessions"):
            for key in (rep.get("accession"), rep.get("current_accession")):
                g = acc_to_gid.get(key)
                if g:
                    by_gid[g] = better(by_gid.get(g), rep)
                    break
        print(f"  accession pass {min(i + a.batch, len(accs)):>6,}/{len(accs):,}"
              f"  resolved {len(by_gid):,}", end="\r", flush=True)
        time.sleep(a.sleep)
    print()

    # pass 2 — BioSample, for everything pass 1 did not resolve
    todo = [g for g in gids if g not in by_gid and g in gid_bio]
    bio_to_gid = {gid_bio[g]: g for g in todo}
    bios = list(bio_to_gid)
    print(f"  unresolved after accession pass: {len(todo):,}")
    for i in range(0, len(bios), a.batch):
        for rep in post(bios[i:i + a.batch], "biosample_ids"):
            bs = ((rep.get("assembly_info", {}) or {})
                  .get("biosample", {}) or {}).get("accession")
            g = bio_to_gid.get(bs)
            if g:
                by_gid[g] = better(by_gid.get(g), rep)
        print(f"  biosample pass {min(i + a.batch, len(bios)):>6,}/{len(bios):,}"
              f"  resolved {len(by_gid):,}", end="\r", flush=True)
        time.sleep(a.sleep)
    print()

    rows, missing = [], []
    for gid in gids:
        rep = by_gid.get(gid)
        if rep is None:
            missing.append(gid)
            rows.append(dict(sample=gid,
                             assembly_accession=gid_acc.get(gid, ""),
                             verdict="ERROR", partial="yes",
                             error="no_ncbi_report"))
            continue
        m = extract(rep)
        row = dict(sample=gid, assembly_accession=rep.get("accession", ""),
                   error="", **m)
        worst, skipped = "PASS", 0
        for metric in ORDER:
            v = classify(m.get(metric), THRESHOLDS[metric])
            row[f"{metric}_verdict"] = v
            if v == "FAIL":
                worst = "FAIL"
            elif v == "WARN" and worst == "PASS":
                worst = "WARN"
            elif v == "not_evaluated":
                skipped += 1
        row["verdict"] = worst
        row["partial"] = "yes" if skipped else "no"
        row["n_metrics_not_evaluated"] = skipped
        rows.append(row)

    df = pd.DataFrame(rows)
    cols = (["sample", "assembly_accession", "verdict", "partial",
             "n_metrics_not_evaluated", "error"] + ORDER
            + [f"{m}_verdict" for m in ORDER]
            + ["checkm_version", "assembly_level"])
    df = df.reindex(columns=[c for c in cols if c in df.columns])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, sep="\t", index=False)

    n = len(df)
    print(f"\nretrieved  : {n - len(missing):,} of {n:,}")
    if missing:
        print(f"missing    : {len(missing):,} (recorded as ERROR)")
    print("\nverdicts")
    for v, k in df["verdict"].value_counts().items():
        print(f"  {v:14s} {k:6,}  ({100 * k / n:5.1f}%)")
    print("\nfailing metric")
    any_fail = False
    for m in ORDER:
        c = f"{m}_verdict"
        if c in df:
            k = int((df[c] == "FAIL").sum())
            if k:
                any_fail = True
                print(f"  {m:20s} {k:6,}")
    if not any_fail:
        print("  (none)")
    print("\nnot evaluated")
    for m in ORDER:
        c = f"{m}_verdict"
        if c in df:
            k = int((df[c] == "not_evaluated").sum())
            if k:
                print(f"  {m:20s} {k:6,}")
    print(f"\n-> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
