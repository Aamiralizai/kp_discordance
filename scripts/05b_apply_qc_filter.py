#!/usr/bin/env python3
"""
05b_apply_qc_filter.py — drop QC-failed genomes from the genotype table.

Quarantining an assembly moves a file. It does not remove that genome from
data/processed/kleborate_klebsiella_pneumo_complex_output.tsv, which is what
08b_discordance_recordlevel.py actually reads. Without this step, running
00b_qc_assembly.py --apply changes nothing downstream and the analysis still
includes every genome that failed quality control.

Retyping the survivors would also work, but Kleborate's output for a genome is
a function of its assembly, so a genome that was typed once does not need
typing again. Filtering the existing table is equivalent and takes seconds
rather than days.

IDEMPOTENT. The first run copies the genotype table to
<name>.prefilter.tsv and filters from it. Every later run filters from that
same backup, so re-running after a changed QC threshold gives the same answer
as running once, rather than compounding.

USAGE
    python scripts/05b_apply_qc_filter.py
    python scripts/05b_apply_qc_filter.py --keep-warn      # default
    python scripts/05b_apply_qc_filter.py --no-keep-warn   # strict: PASS only
    python scripts/05b_apply_qc_filter.py --restore        # undo
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

KLEB = Path("data/processed/kleborate_klebsiella_pneumo_complex_output.tsv")
QC = Path("data/interim/assembly_qc_metrics.tsv")
STRAIN = "strain"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kleborate", type=Path, default=KLEB)
    ap.add_argument("--qc", type=Path, default=QC)
    ap.add_argument("--keep-warn", dest="keep_warn", action="store_true",
                    default=True,
                    help="keep assemblies with a WARN verdict (default)")
    ap.add_argument("--no-keep-warn", dest="keep_warn", action="store_false",
                    help="keep only PASS; use for the QC-stringency sensitivity")
    ap.add_argument("--restore", action="store_true",
                    help="restore the unfiltered table and exit")
    a = ap.parse_args()

    backup = a.kleborate.with_suffix(".prefilter.tsv")

    if a.restore:
        if not backup.exists():
            sys.exit(f"No backup at {backup}; nothing to restore.")
        shutil.copy2(backup, a.kleborate)
        print(f"restored {a.kleborate} from {backup}")
        return 0

    for p in (a.kleborate, a.qc):
        if not p.exists():
            sys.exit(f"Missing {p}.\nRun 00b_qc_assembly.py (and the genotyping "
                     f"steps) before this.")

    # Filter from the pristine copy every time, never from a filtered table.
    if not backup.exists():
        shutil.copy2(a.kleborate, backup)
        print(f"backed up unfiltered genotype table -> {backup}")
    source = backup

    kb = pd.read_csv(source, sep="\t", dtype=str, low_memory=False)
    if STRAIN not in kb.columns:
        sys.exit(f"{source} has no {STRAIN!r} column (found: "
                 f"{', '.join(kb.columns[:8])}...). 08b joins on this column, "
                 f"so the filter cannot be applied safely.")

    qc = pd.read_csv(a.qc, sep="\t", dtype=str, low_memory=False)
    qc["sample"] = qc["sample"].astype(str).str.strip()
    verdicts = dict(zip(qc["sample"], qc["verdict"].astype(str)))

    kb["_k"] = kb[STRAIN].astype(str).str.strip()
    keep_set = {"PASS"} | ({"WARN"} if a.keep_warn else set())

    matched = kb["_k"].isin(verdicts)
    kb["_verdict"] = kb["_k"].map(verdicts).fillna("no_qc_result")
    keep = kb["_verdict"].isin(keep_set)

    n_before = len(kb)
    n_unmatched = int((~matched).sum())
    partial = (qc["partial"].astype(str) == "yes").any() if "partial" in qc else False

    print(f"genotype table   : {n_before:,} genomes")
    print(f"QC results       : {len(verdicts):,} assemblies")
    print(f"keeping verdicts : {', '.join(sorted(keep_set))}")
    print("\nverdict breakdown in the genotype table:")
    for v, n in kb["_verdict"].value_counts().items():
        mark = "keep" if v in keep_set else "DROP"
        print(f"  {v:16s} {n:7,}   {mark}")

    if n_unmatched:
        print(f"\n{n_unmatched:,} genomes in the genotype table have no QC "
              f"result and are DROPPED.\nA genome that was typed but never "
              f"quality-controlled cannot be reported as\nhaving passed. If "
              f"this count is large, the sample names in the QC table do not\n"
              f"match the {STRAIN!r} column and the filter is mismatched — "
              f"check before proceeding.")

    out = kb[keep].drop(columns=["_k", "_verdict"])
    out.to_csv(a.kleborate, sep="\t", index=False)

    print(f"\nkept             : {len(out):,}")
    print(f"removed          : {n_before - len(out):,}")
    print(f"written          : {a.kleborate}")
    print(f"unfiltered copy  : {backup}  (restore with --restore)")

    if partial:
        print("\nNOTE: the QC verdicts are PARTIAL — completeness, contamination")
        print("or CDS count were not evaluated. Supply them before describing")
        print("this as a complete QualiBact assessment.")

    print("\nNext: bash run_all.sh")
    return 0


if __name__ == "__main__":
    sys.exit(main())
