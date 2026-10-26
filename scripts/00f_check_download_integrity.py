#!/usr/bin/env python3
"""00f_check_download_integrity.py — is every downloaded assembly complete?

WHY THIS EXISTS
A bulk FTPS download of 7,360 assemblies left 3,453 of them (46.9%) truncated:
the transfer stopped in mid-sequence and the client still exited successfully.
A truncated FASTA is still a valid FASTA, so it passed 00_qc_fasta.sh, and the
QualiBact metrics applied in 00b/00e describe the assembly as deposited in
BV-BRC rather than the copy on disk, so it passed quality control too. The
genotyper then typed the incomplete file and missed the genes on the contigs
that never arrived.

Neither existing check can catch this. Only a comparison against the deposited
record can, which is what this script does.

WHAT IT CHECKS, PER GENOME
    total sequence length == BV-BRC genome_length
    number of records     == BV-BRC contigs
    the file ends with a newline (a cut transfer almost never does)

A genome is reported as truncated when the length is short or the final newline
is missing. Genomes with no BV-BRC record are reported as such rather than
silently passed.

USAGE (from the project root)
    python3 scripts/00f_check_download_integrity.py \
        --genomes data/raw/genomes \
        --out results/qc_concordance/fasta_integrity_report.tsv

    --ids FILE   check only these genome ids (default: every .fna found)
    --jobs N     parallel FASTA readers (default: cpu_count)

OUTPUT
    fasta_integrity_report.tsv   one row per genome
    truncated_ids.txt            written next to the report, the ids to repair

REPAIRING WHAT THIS FINDS
    python3 scripts/01d_redownload_verified.py \
        results/qc_concordance/truncated_ids.txt data/raw/genomes_redownload
    bash scripts/04c_fix_truncated_genomes.sh
"""
import argparse
import csv
import json
import os
import subprocess
import shutil
import sys
import time
import urllib.request
from multiprocessing import Pool

API = "https://www.bv-brc.org/api"
UA = "kp-discordance/1.0 (download integrity check)"


def local_metrics(path):
    """(genome_id, total_bp, n_records, ends_with_newline)."""
    gid = os.path.basename(path)[:-4]
    try:
        with open(path, "rb") as fh:
            data = fh.read()
    except OSError:
        return gid, 0, 0, 0
    total = 0
    recs = data.split(b">")[1:]
    for rec in recs:
        nl = rec.find(b"\n")
        if nl < 0:
            continue
        total += len(rec[nl + 1:].replace(b"\n", b"").replace(b"\r", b""))
    return gid, total, len(recs), int(data.endswith(b"\n"))


def http_get(url, timeout=120):
    if shutil.which("curl"):
        r = subprocess.run(["curl", "-sS", "-m", str(timeout), "-A", UA,
                            "-H", "Accept: application/json", url],
                           capture_output=True)
        return r.stdout if r.returncode == 0 else b""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA,
                                                   "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read()
    except Exception:
        return b""


def deposited(ids):
    """BV-BRC genome_length and contigs for each id."""
    meta = {}
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        q = (f"in(genome_id,({','.join(chunk)}))"
             f"&select(genome_id,genome_length,contigs)&limit(1000)")
        for attempt in range(4):
            try:
                for row in json.loads(http_get(f"{API}/genome/?{q}")):
                    meta[row["genome_id"]] = row
                break
            except Exception:
                time.sleep(5 * (attempt + 1))
        print(f"  deposited metadata: {len(meta):,}/{len(ids):,}", flush=True)
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", action="append", default=None,
                    help="genome directory; repeat it to audit the EFFECTIVE set, "
                         "where a later directory overrides an earlier one, e.g. "
                         "--genomes data/raw/genomes "
                         "--genomes data/raw/genomes_redownload")
    ap.add_argument("--out",
                    default="results/qc_concordance/fasta_integrity_report.tsv")
    ap.add_argument("--ids")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args()

    dirs = a.genomes or ["data/raw/genomes"]
    # Later directories override earlier ones, so the audit describes the file
    # that was actually genotyped rather than the first copy retrieved.
    chosen = {}
    for d in dirs:
        for f in os.listdir(d):
            if f.endswith(".fna"):
                chosen[f[:-4]] = os.path.join(d, f)
    files = [chosen[g] for g in sorted(chosen)]
    if a.ids:
        want = {l.strip() for l in open(a.ids) if l.strip()}
        files = [f for f in files if os.path.basename(f)[:-4] in want]
    if not files:
        sys.exit(f"no .fna files under {a.genomes}")
    print(f"{len(files):,} assemblies to check", flush=True)

    with Pool(a.jobs) as pool:
        local = pool.map(local_metrics, files, chunksize=32)
    paths = {os.path.basename(f)[:-4]: f for f in files}
    meta = deposited([g for g, *_ in local])

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    truncated, no_record = [], 0
    with open(a.out, "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["genome_id", "source_dir", "local_size_bp", "local_contigs",
                     "ends_newline", "bvbrc_size_bp", "bvbrc_contigs",
                     "size_ratio", "verdict"])
        for gid, size, nrec, newline in local:
            m = meta.get(gid)
            exp = int(m.get("genome_length") or 0) if m else 0
            exp_ctg = int(m.get("contigs") or 0) if m else 0
            ratio = f"{size / exp:.4f}" if exp else ""
            # A missing deposited record is UNVERIFIED, never COMPLETE: the
            # comparison could not be made. The test is exact agreement on
            # length AND contig count, not "not shorter": a file that ends
            # cleanly can still be short, which an earlier newline-only rule
            # let through.
            if not m:
                no_record += 1
                verdict = "UNVERIFIED"
            elif size != exp or (exp_ctg and nrec != exp_ctg) or not newline:
                verdict = "INCOMPLETE"
                truncated.append(gid)
            else:
                verdict = "COMPLETE"
            wr.writerow([gid, os.path.basename(os.path.dirname(paths[gid])),
                         size, nrec, newline, exp or "", exp_ctg or "",
                         ratio, verdict])

    ids_path = os.path.join(os.path.dirname(a.out) or ".", "truncated_ids.txt")
    with open(ids_path, "w") as fh:
        fh.write("\n".join(truncated) + ("\n" if truncated else ""))

    n = len(local)
    print(f"\n  COMPLETE   : {n - len(truncated) - no_record:,}")
    print(f"  INCOMPLETE : {len(truncated):,} ({100*len(truncated)/n:.1f}%)")
    print(f"  UNVERIFIED : {no_record:,}  (no deposited record to compare against)")
    print(f"\n  report : {a.out}")
    print(f"  repair : {ids_path}")
    if truncated:
        print("\n  Next: scripts/01d_redownload_verified.py then "
              "scripts/04c_fix_truncated_genomes.sh")
    sys.exit(1 if truncated else 0)


if __name__ == "__main__":
    main()
