#!/usr/bin/env python3
"""
19a_fetch_derivation_accessions.py — NCBI accessions for the derivation set.

External validation is only external if the two sets do not share isolates.
BV-BRC genome identifiers (for example 573.12771) are not NCBI accessions, so
the overlap cannot be computed from the analysis files alone. This script
resolves every derivation genome to its NCBI assembly and BioSample accession
through the BV-BRC API, so that 19_external_validation.py can exclude them.

Without this step any apparent external performance is partly a re-measurement
of the training data, and should not be reported as validation.

USAGE
    python scripts/19a_fetch_derivation_accessions.py
    python scripts/19a_fetch_derivation_accessions.py --batch 100 --sleep 0.5
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

API = "https://www.bv-brc.org/api/genome/"
SUP = Path("results/supplementary")

# The API returns HTTP 403 to the default Python-urllib user agent. Identify the
# pipeline explicitly: it is accepted, and it tells the service operator who is
# making the request.
UA = ("kp_discordance/1.0 (research pipeline; "
      "+https://github.com/Aamiralizai/kp_discordance)")


def fetch(ids: list[str], retries: int = 3, timeout: int = 90) -> list[dict]:
    """One batch of genome_ids -> accession records.

    The query goes in the URL rather than a POST body: the RQL POST content
    type is rejected (HTTP 403) through some network paths, while the
    equivalent GET is accepted. Keep batches small enough that the URL stays
    well inside server limits.
    """
    q = ("in(genome_id,(" + ",".join(ids) + "))"
         "&select(genome_id,assembly_accession,biosample_accession,"
         "genome_name,sra_accession)"
         "&limit(" + str(len(ids) + 10) + ")"
         "&http_accept=application/json")
    url = API + "?" + q
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
    ap.add_argument("--data", type=Path, default=SUP / "S1_genome_level_data.tsv")
    ap.add_argument("--out", type=Path, default=SUP / "S1_derivation_accessions.tsv")
    ap.add_argument("--batch", type=int, default=100)
    ap.add_argument("--sleep", type=float, default=0.3,
                    help="pause between batches, to stay polite to the API")
    a = ap.parse_args()

    if not a.data.exists():
        sys.exit(f"Missing {a.data}")
    ids = (pd.read_csv(a.data, sep="\t", dtype=str, low_memory=False)
           ["BVBRC_genome_id"].dropna().astype(str).unique().tolist())
    print(f"derivation genomes : {len(ids):,}")

    rows, missing = [], []
    for i in range(0, len(ids), a.batch):
        chunk = ids[i:i + a.batch]
        got = fetch(chunk)
        rows += got
        found = {r.get("genome_id") for r in got}
        missing += [g for g in chunk if g not in found]
        print(f"  {min(i + a.batch, len(ids)):>6,}/{len(ids):,}  "
              f"resolved {len(rows):,}", end="\r", flush=True)
        time.sleep(a.sleep)
    print()

    df = pd.DataFrame(rows)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, sep="\t", index=False)

    n_asm = int(df["assembly_accession"].notna().sum()) if "assembly_accession" in df else 0
    n_bio = int(df["biosample_accession"].notna().sum()) if "biosample_accession" in df else 0
    print(f"resolved           : {len(df):,}")
    print(f"  assembly acc     : {n_asm:,}")
    print(f"  biosample acc    : {n_bio:,}")
    print(f"unresolved         : {len(missing):,}")
    if missing:
        print(f"  first few        : {', '.join(missing[:8])}")
        print("  Unresolved genomes cannot be excluded by accession. Report the")
        print("  count: it bounds the residual overlap in the validation set.")
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
