#!/usr/bin/env python
"""05_merge_kleborate.py — merge Kleborate shard outputs correctly.

Kleborate 3.x writes SEVERAL output files per run (one per module group).
Naively concatenating every .txt inflates row counts. This script:
  * groups shard files by filename (so each module's table stays separate)
  * concatenates same-named files across shards, keeping ONE header
  * de-duplicates on the strain/genome column
  * reports what it found so you can see the module tables

Run (either env, needs pandas):
    python scripts/05_merge_kleborate.py

Outputs:
    data/processed/kleborate_<filename>.tsv   one per module table
    data/processed/kleborate_results.tsv      the main KpSC table (primary)
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

SHARD_ROOT = Path("data/interim/kleborate")
OUT = Path("data/processed")
SKIP = {"kleborate.log", "failed_chunks.txt", "done.txt"}


def main():
    if not SHARD_ROOT.exists():
        sys.exit(f"ERROR: {SHARD_ROOT} not found — run 04_run_kleborate.sh first.")

    groups = defaultdict(list)
    for shard in sorted(SHARD_ROOT.glob("shard_*")):
        for f in sorted(shard.glob("*.txt")):
            if f.name in SKIP or f.stat().st_size == 0:
                continue
            groups[f.name].append(f)

    if not groups:
        sys.exit(f"ERROR: no output .txt files under {SHARD_ROOT}/shard_*/")

    print(f"Found {len(groups)} distinct output table(s) across shards:\n")
    written = []
    for name, files in sorted(groups.items()):
        frames = []
        for f in files:
            try:
                frames.append(pd.read_csv(f, sep="\t", dtype=str, low_memory=False))
            except Exception as e:
                print(f"  ! skipping unreadable {f}: {e}")
        if not frames:
            continue
        df = pd.concat(frames, ignore_index=True)

        # de-duplicate on the strain/genome identifier column
        key = next((c for c in df.columns
                    if c.lower().strip() in {"strain", "genome", "genome_id", "assembly", "name"}), None)
        n_before = len(df)
        if key:
            df = df.drop_duplicates(subset=[key], keep="first")
        else:
            df = df.drop_duplicates(keep="first")

        stem = Path(name).stem
        dest = OUT / f"kleborate_{stem}.tsv"
        df.to_csv(dest, sep="\t", index=False)
        written.append((stem, dest, len(df), key))
        print(f"  {name:<50} {n_before:>6} rows -> {len(df):>6} unique "
              f"(key: {key or 'n/a'}) [{len(files)} shard files]")

    # promote the main KpSC genotype table (one row per genome) to the canonical name.
    # NOTE: hAMRonization output is one row PER GENE HIT, not per genome — never primary.
    candidates = [(stem, dest, n, key) for stem, dest, n, key in written
                  if "hamronization" not in stem.lower() and key is not None]
    main = None
    # must be the K. PNEUMONIAE complex table — not oxytoca, not escherichia
    for stem, dest, n, key in candidates:
        s = stem.lower()
        if ("pneumo" in s or "kpsc" in s) and "oxytoca" not in s and "escherichia" not in s:
            main = (dest, n); break
    if main is None and candidates:
        main = max(((d, n) for _, d, n, _ in candidates), key=lambda x: x[1])
    if main is None:
        print("\n! No per-genome table found (none had a strain/genome key).")
        print("! Check the shard outputs manually.")

    if main:
        dest, n = main
        canonical = OUT / "kleborate_results.tsv"
        pd.read_csv(dest, sep="\t", dtype=str, low_memory=False).to_csv(
            canonical, sep="\t", index=False)
        print(f"\nPrimary table -> {canonical}  ({n} genomes)")
        print("(this is the one downstream scripts read)")


if __name__ == "__main__":
    main()
