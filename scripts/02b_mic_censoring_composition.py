#!/usr/bin/env python
"""02b_mic_censoring_composition.py — censoring composition of the MIC records.

The Methods state what share of the retrieved minimum inhibitory concentration
records are left-censored, right-censored and exact. That claim had no released
table behind it, so it is computed here and written out.

A record counts when a numeric concentration can be read from it and the unit is
not millimetres, which excludes disk-diffusion zone diameters. The sign is taken
from the sign field, or from the leading operator of the measurement string when
that field is empty.

Usage:
    python scripts/02b_mic_censoring_composition.py
"""
import argparse
import csv
import re
from pathlib import Path

RAW = Path("data/raw/phenotypes/kp_amr_anchor.tsv")
OUT = Path("results/tables/mic_censoring_composition.tsv")
OP = re.compile(r"\s*(<=|>=|<|>|=)")


def sign_of(row):
    s = (row.get("measurement_sign") or "").strip()
    if s:
        return s
    m = OP.match((row.get("measurement") or ""))
    return m.group(1) if m else ""


def value_of(row):
    v = (row.get("measurement_value") or "").strip()
    if not v:
        v = OP.sub("", (row.get("measurement") or "").strip(), count=1)
    try:
        return float(v)
    except ValueError:
        return None


def kind_of(sign):
    if sign in ("<=", "<"):
        return "left_censored"
    if sign in (">", ">="):
        return "right_censored"
    return "exact"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drugs", default="meropenem,imipenem")
    a = ap.parse_args()
    drugs = [d.strip().lower() for d in a.drugs.split(",") if d.strip()]

    with open(RAW, newline="", encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["drug", "n_records_retrieved", "n_concentration_records",
                     "censoring", "n", "pct"])
        for drug in drugs:
            got = [r for r in rows
                   if (r.get("antibiotic") or "").strip().lower() == drug]
            conc = [r for r in got
                    if value_of(r) is not None
                    and (r.get("measurement_unit") or "").strip() != "mm"]
            total = len(conc) or 1
            counts = {"left_censored": 0, "right_censored": 0, "exact": 0}
            for r in conc:
                counts[kind_of(sign_of(r))] += 1
            for kind in ("left_censored", "right_censored", "exact"):
                wr.writerow([drug, len(got), len(conc), kind, counts[kind],
                             round(100 * counts[kind] / total, 4)])
            print(f"  {drug}: {len(conc)} concentration records of {len(got)} "
                  f"retrieved — " + ", ".join(
                      f"{k.replace('_', ' ')} {100 * counts[k] / total:.1f}%"
                      for k in ("left_censored", "right_censored", "exact")))
    print(f"[out] {OUT}")


if __name__ == "__main__":
    main()
