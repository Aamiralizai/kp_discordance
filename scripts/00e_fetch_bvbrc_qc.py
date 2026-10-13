#!/usr/bin/env python3
"""
00e_fetch_bvbrc_qc.py — assembly QC metrics from BV-BRC.

The preferred no-download route when the assemblies came from BV-BRC, which is
the case for this study. BV-BRC publishes, per genome, the complete QualiBact
metric set computed on the same assembly the pipeline genotyped.

    genome_length         genome size
    gc_content            GC content
    contig_n50, contigs   contiguity
    checkm_completeness   CheckM completeness
    checkm_contamination  CheckM contamination
    patric_cds            coding sequences

WHY THIS RATHER THAN NCBI (00d_fetch_ncbi_qc.py)

Two reasons, and the second is the important one.

Coverage. Roughly a tenth of these genomes have no retrievable NCBI assembly
record at all, and about half of the rest are GenBank-only submissions with no
RefSeq counterpart and therefore no CheckM result. Querying NCBI leaves both
gaps unfilled.

Correctness. NCBI's statistics describe the record deposited at NCBI. BV-BRC's
describe the assembly BV-BRC served, which is the one that was downloaded and
typed. Where a genome has been suppressed, replaced or revised at NCBI since
BV-BRC mirrored it, those are different sequences, and the genomes where they
differ are exactly the ones most likely to be of poor quality. Measuring the
wrong copy would hide the problem.

00b_qc_assembly.py remains the direct route. It computes the five sequence
metrics from the local FASTA and should be run on a sample to confirm these
agree.

Output matches 00b_qc_assembly.py, so 05b_apply_qc_filter.py and
18_discordance_investigation.py consume any of the three interchangeably.

USAGE
    python scripts/00e_fetch_bvbrc_qc.py
    python scripts/00e_fetch_bvbrc_qc.py --ids results/supplementary/S1_genome_level_data.tsv
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

API = "https://www.bv-brc.org/api/genome/"
UA = ("kp_discordance/1.0 (research pipeline; "
      "+https://github.com/Aamiralizai/kp_discordance)")

FIELDS = ("genome_id,genome_status,genome_quality,contigs,contig_n50,"
          "genome_length,gc_content,checkm_completeness,checkm_contamination,"
          "patric_cds,assembly_accession,sequencing_centers,plasmids")

# QualiBact v1.0, Klebsiella pneumoniae. config/qc_thresholds.csv is the source
# of truth; these mirror it for the no-download path.
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

MAP = {"genome_length": "genome_size_bp", "gc_content": "gc_content_pct",
       "contig_n50": "n50_bp", "contigs": "n_contigs",
       "checkm_completeness": "completeness_pct",
       "checkm_contamination": "contamination_pct", "patric_cds": "n_cds"}


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


def fetch(ids: list[str], retries: int = 4, timeout: int = 90) -> list[dict]:
    """GET, not POST. The RQL POST content type is rejected (HTTP 403) through
    some network paths, and the default Python user agent is rejected outright."""
    url = (API + "?in(genome_id,(" + ",".join(ids) + "))"
           + "&select(" + FIELDS + ")"
           + f"&limit({len(ids) + 10})&http_accept=application/json")
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json", "User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return json.loads(fh.read().decode())
        except (urllib.error.URLError, urllib.error.HTTPError,
                TimeoutError, json.JSONDecodeError) as exc:
            if attempt == retries:
                print(f"\n  batch failed after {retries} attempts: {exc}",
                      file=sys.stderr)
                return []
            time.sleep(2 * attempt)
    return []


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", type=Path,
                    default=Path("results/supplementary/S1_genome_level_data.tsv"))
    ap.add_argument("--id-column", default="BVBRC_genome_id")
    ap.add_argument("--out", type=Path,
                    default=Path("data/interim/assembly_qc_metrics.tsv"))
    ap.add_argument("--batch", type=int, default=100)
    ap.add_argument("--sleep", type=float, default=0.2)
    a = ap.parse_args()

    if not a.ids.exists():
        sys.exit(f"Missing {a.ids}")
    src = pd.read_csv(a.ids, sep="\t", dtype=str, low_memory=False)
    col = a.id_column if a.id_column in src.columns else "genome_id"
    if col not in src.columns:
        sys.exit(f"{a.ids} has neither {a.id_column!r} nor 'genome_id'.")
    gids = src[col].dropna().astype(str).unique().tolist()
    print(f"genomes : {len(gids):,}")

    got: dict[str, dict] = {}
    for i in range(0, len(gids), a.batch):
        for rec in fetch(gids[i:i + a.batch]):
            g = str(rec.get("genome_id", ""))
            if g:
                got[g] = rec
        print(f"  {min(i + a.batch, len(gids)):>6,}/{len(gids):,}  "
              f"retrieved {len(got):,}", end="\r", flush=True)
        time.sleep(a.sleep)
    print()

    rows, missing = [], []
    for g in gids:
        rec = got.get(g)
        if rec is None:
            missing.append(g)
            rows.append(dict(sample=g, verdict="ERROR", partial="yes",
                             error="no_bvbrc_record"))
            continue

        def num(k):
            v = rec.get(k)
            try:
                return float(v)
            except (TypeError, ValueError):
                return None

        m = {dst: num(srcf) for srcf, dst in MAP.items()}
        row = dict(sample=g, error="",
                   genome_status=rec.get("genome_status", ""),
                   bvbrc_genome_quality=rec.get("genome_quality", ""),
                   assembly_accession=rec.get("assembly_accession", "") or "",
                   sequencing_centers=rec.get("sequencing_centers", "") or "",
                   plasmids=rec.get("plasmids", ""), **m)
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
    cols = (["sample", "verdict", "partial", "n_metrics_not_evaluated", "error",
             "genome_status", "bvbrc_genome_quality", "assembly_accession",
             "sequencing_centers", "plasmids"] + ORDER
            + [f"{m}_verdict" for m in ORDER])
    df = df.reindex(columns=[c for c in cols if c in df.columns])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, sep="\t", index=False)

    n = len(df)
    print(f"\nretrieved : {n - len(missing):,} of {n:,}")
    if missing:
        print(f"missing   : {len(missing):,} (recorded as ERROR)")
    print("\nverdicts")
    for v, k in df["verdict"].value_counts().items():
        print(f"  {v:8s} {k:6,}  ({100 * k / n:5.1f}%)")
    print("\nfailing metric")
    for m in ORDER:
        c = f"{m}_verdict"
        if c in df:
            k = int((df[c] == "FAIL").sum())
            if k:
                print(f"  {m:20s} {k:6,}")
    ne = {m: int((df[f"{m}_verdict"] == "not_evaluated").sum())
          for m in ORDER if f"{m}_verdict" in df}
    ne = {k: v for k, v in ne.items() if v}
    if ne:
        print("\nnot evaluated")
        for m, k in ne.items():
            print(f"  {m:20s} {k:6,}")
    print(f"\n-> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
