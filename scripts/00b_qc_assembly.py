#!/usr/bin/env python3
"""
00b_qc_assembly.py — species-specific assembly QUALITY control.

This is not the same step as scripts/00_qc_fasta.sh, and the distinction is the
point. That script checks whether a file is a structurally valid FASTA, because
Kleborate aborts an entire batch on a malformed record. Structural validity says
nothing about whether an assembly is good enough to genotype reliably. This
script asks the second question, against published species-specific criteria.

Why it matters for this study specifically, rather than as general hygiene:

  * A fragmented assembly can split or truncate ompK36. The truncation is then
    read as porin disruption, which is one of the two exposures under study.
    Fragmentation therefore manufactures the exposure.
  * A fragmented or incomplete assembly can drop a carbapenemase gene entirely,
    which manufactures the other exposure in the opposite direction.
  * A contaminated assembly can contribute a carbapenemase belonging to another
    organism, which is a false positive that no downstream filter would catch.

So assembly quality is not a background covariate here. It is confounded with
both exposures, in opposite directions, and must be controlled before typing.

METRICS

Five are computed directly from the assembly FASTA:
    genome_size_bp, gc_content_pct, n50_bp, n_contigs, longest_contig_bp

Three require external tools and are evaluated only when supplied:
    completeness_pct, contamination_pct   --checkm2  CheckM2 quality_report.tsv
    n_cds                                 --cds      TSV of sample<TAB>cds_count

Any metric that cannot be evaluated is reported as "not_evaluated" per assembly
and listed in the summary. A verdict computed without the full metric set is
labelled PARTIAL so that no run can silently claim a completeness check it did
not perform.

VERDICT

Per metric: FAIL / WARN / PASS against config/qc_thresholds.csv.
Overall:   FAIL if any metric fails, else WARN if any warns, else PASS.
           Suffixed "(partial)" when some metrics were not evaluable.

USAGE

    # report only
    python scripts/00b_qc_assembly.py

    # with CheckM2 completeness and contamination
    python scripts/00b_qc_assembly.py --checkm2 data/interim/checkm2/quality_report.tsv

    # quarantine failures so they are excluded from genotyping
    python scripts/00b_qc_assembly.py --apply

    # lenient comparator, for the QC-stringency sensitivity analysis
    python scripts/00b_qc_assembly.py --scheme enterobase-v2.3 --tag _enterobase
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

GEN_DEFAULT = Path("data/raw/genomes")
THRESHOLDS = Path("config/qc_thresholds.csv")

FASTA_METRICS = ["genome_size_bp", "gc_content_pct", "n50_bp",
                 "n_contigs", "longest_contig_bp"]
TOOL_METRICS = {"completeness_pct": "checkm2",
                "contamination_pct": "checkm2",
                "n_cds": "annotation"}
ALL_METRICS = FASTA_METRICS + list(TOOL_METRICS)


# --------------------------------------------------------------- thresholds
def load_thresholds(scheme: str) -> dict:
    """Read the verified bounds for one scheme. Refuses unverified rows."""
    if not THRESHOLDS.exists():
        sys.exit(f"Missing {THRESHOLDS}. QC thresholds are not hard-coded; the "
                 f"table must be present and verified.")
    out, seen_schemes = {}, set()
    with THRESHOLDS.open() as fh:
        for row in csv.DictReader(r for r in fh if not r.lstrip().startswith("#")):
            seen_schemes.add(row["scheme"])
            if row["scheme"] != scheme:
                continue
            if row["verified"].strip().upper() != "YES":
                sys.exit(f"Threshold row for {scheme}/{row['metric']} is not "
                         f"verified. Verify it against the source document "
                         f"before analysing.")

            def num(k):
                v = row[k].strip()
                return float(v) if v else None

            out[row["metric"]] = dict(
                fail_lower=num("fail_lower"), warn_lower=num("warn_lower"),
                warn_upper=num("warn_upper"), fail_upper=num("fail_upper"),
                needs_tool=row["needs_tool"].strip() or None,
                source=row["source"].strip())
    if not out:
        sys.exit(f"Scheme {scheme!r} not found in {THRESHOLDS}. "
                 f"Available: {', '.join(sorted(seen_schemes))}")
    return out


def classify(value, t) -> str:
    """FAIL / WARN / PASS for one metric. Bounds are inclusive of the good side."""
    if value is None:
        return "not_evaluated"
    if t["fail_lower"] is not None and value < t["fail_lower"]:
        return "FAIL"
    if t["fail_upper"] is not None and value > t["fail_upper"]:
        return "FAIL"
    if t["warn_lower"] is not None and value < t["warn_lower"]:
        return "WARN"
    if t["warn_upper"] is not None and value > t["warn_upper"]:
        return "WARN"
    return "PASS"


# ------------------------------------------------------------------ metrics
def assembly_metrics(path: Path) -> dict:
    """Compute assembly statistics in a single pass over the file.

    GC is taken over unambiguous bases only, so that a run of Ns depresses the
    ambiguous fraction rather than silently dragging GC toward zero.
    """
    lengths: list[int] = []
    cur = 0
    a = c = g = t = n = other = 0
    try:
        with path.open("rb") as fh:
            for raw in fh:
                if raw.startswith(b">"):
                    if cur:
                        lengths.append(cur)
                    cur = 0
                    continue
                line = raw.strip()
                if not line:
                    continue
                cur += len(line)
                low = line.lower()
                a += low.count(b"a")
                c += low.count(b"c")
                g += low.count(b"g")
                t += low.count(b"t")
                n += low.count(b"n")
    except OSError as exc:
        return {"sample": path.stem, "error": str(exc)}
    if cur:
        lengths.append(cur)

    if not lengths:
        return {"sample": path.stem, "error": "no_sequence"}

    total = sum(lengths)
    lengths.sort(reverse=True)
    half, run, n50 = total / 2.0, 0, lengths[-1]
    for L in lengths:
        run += L
        if run >= half:
            n50 = L
            break
    acgt = a + c + g + t
    other = total - acgt - n

    return {
        "sample": path.stem,
        "error": "",
        "genome_size_bp": total,
        "n_contigs": len(lengths),
        "n50_bp": n50,
        "longest_contig_bp": lengths[0],
        "gc_content_pct": round(100.0 * (g + c) / acgt, 4) if acgt else None,
        # reported, not thresholded: high ambiguity is worth seeing even though
        # QualiBact does not bound it.
        "ambiguous_bases": n + other,
        "ambiguous_pct": round(100.0 * (n + other) / total, 4) if total else None,
    }


# ------------------------------------------------------------- external data
def load_checkm2(path: Path) -> dict:
    """completeness/contamination keyed by sample name, from CheckM2."""
    out = {}
    with path.open() as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = {k.lower(): k for k in (rd.fieldnames or [])}
        name = cols.get("name") or cols.get("bin id") or cols.get("genome")
        comp = cols.get("completeness") or cols.get("completeness_general")
        cont = cols.get("contamination")
        if not (name and comp and cont):
            sys.exit(f"{path} does not look like a CheckM2 quality_report.tsv "
                     f"(need Name/Completeness/Contamination; got "
                     f"{rd.fieldnames}).")
        for row in rd:
            try:
                out[row[name]] = (float(row[comp]), float(row[cont]))
            except (TypeError, ValueError):
                continue
    return out


def load_cds(path: Path) -> dict:
    out = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                out[parts[0]] = int(float(parts[1]))
            except ValueError:
                continue
    return out


# --------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", type=Path, default=GEN_DEFAULT)
    ap.add_argument("--glob", default="*.fna",
                    help="assembly filename pattern (default: *.fna)")
    ap.add_argument("--scheme", default="qualibact-v1.0")
    ap.add_argument("--checkm2", type=Path,
                    help="CheckM2 quality_report.tsv")
    ap.add_argument("--cds", type=Path,
                    help="TSV: sample<TAB>cds_count")
    ap.add_argument("--apply", action="store_true",
                    help="move FAIL assemblies to the quarantine directory")
    ap.add_argument("--quarantine", type=Path,
                    default=Path("data/interim/failed_assembly_qc"))
    ap.add_argument("--tag", default="",
                    help="suffix for the output filenames")
    ap.add_argument("--jobs", type=int, default=0,
                    help="worker processes (default: cpu_count)")
    ap.add_argument("--outdir", type=Path, default=Path("data/interim"))
    a = ap.parse_args()

    thr = load_thresholds(a.scheme)

    files = sorted(a.genomes.glob(a.glob))
    if not files:
        sys.exit(f"No assemblies matching {a.glob} under {a.genomes}")

    print(f"scheme    : {a.scheme}")
    print(f"assemblies: {len(files):,} under {a.genomes}")
    print(f"metrics   : {', '.join(m for m in ALL_METRICS if m in thr)}")

    jobs = a.jobs or None
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        rows = list(pool.map(assembly_metrics, files, chunksize=16))

    checkm = load_checkm2(a.checkm2) if a.checkm2 else {}
    cds = load_cds(a.cds) if a.cds else {}
    if a.checkm2:
        print(f"checkm2   : {len(checkm):,} records from {a.checkm2}")
    if a.cds:
        print(f"cds       : {len(cds):,} records from {a.cds}")

    # which thresholded metrics can we actually evaluate?
    evaluable, skipped = [], []
    for m in ALL_METRICS:
        if m not in thr:
            continue
        need = thr[m]["needs_tool"]
        if need == "checkm2" and not checkm:
            skipped.append(m)
        elif need == "annotation" and not cds:
            skipped.append(m)
        else:
            evaluable.append(m)

    a.outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = a.outdir / f"assembly_qc_metrics{a.tag}.tsv"
    summary_path = a.outdir / f"assembly_qc_summary{a.tag}.txt"

    verdicts = {"PASS": 0, "WARN": 0, "FAIL": 0, "ERROR": 0}
    fail_reasons: dict[str, int] = {}
    warn_reasons: dict[str, int] = {}
    quarantined = 0
    if a.apply:
        a.quarantine.mkdir(parents=True, exist_ok=True)

    cols = (["sample", "verdict", "partial", "error"] + ALL_METRICS
            + ["ambiguous_bases", "ambiguous_pct"]
            + [f"{m}_verdict" for m in ALL_METRICS])

    with metrics_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore", restval="")
        w.writeheader()

        for row, path in zip(rows, files):
            if row.get("error"):
                row["verdict"] = "ERROR"
                row["partial"] = "yes"
                verdicts["ERROR"] += 1
                w.writerow(row)
                continue

            sample = row["sample"]
            if sample in checkm:
                row["completeness_pct"], row["contamination_pct"] = checkm[sample]
            if sample in cds:
                row["n_cds"] = cds[sample]

            worst = "PASS"
            for m in ALL_METRICS:
                if m not in thr:
                    continue
                v = classify(row.get(m), thr[m])
                row[f"{m}_verdict"] = v
                if v == "FAIL":
                    worst = "FAIL"
                    fail_reasons[m] = fail_reasons.get(m, 0) + 1
                elif v == "WARN":
                    warn_reasons[m] = warn_reasons.get(m, 0) + 1
                    if worst == "PASS":
                        worst = "WARN"

            row["verdict"] = worst
            row["partial"] = "yes" if skipped else "no"
            verdicts[worst] += 1
            w.writerow(row)

            if worst == "FAIL" and a.apply:
                try:
                    shutil.move(str(path), str(a.quarantine / path.name))
                    quarantined += 1
                except OSError:
                    pass

    # ------------------------------------------------------------- summary
    n = len(files)
    lines = [
        "assembly quality control",
        f"scheme        : {a.scheme}",
        f"thresholds    : {THRESHOLDS}",
        f"assemblies    : {n:,}",
        "",
        "verdicts",
        f"  PASS        : {verdicts['PASS']:,} ({100*verdicts['PASS']/n:.1f}%)",
        f"  WARN        : {verdicts['WARN']:,} ({100*verdicts['WARN']/n:.1f}%)",
        f"  FAIL        : {verdicts['FAIL']:,} ({100*verdicts['FAIL']/n:.1f}%)",
        f"  unreadable  : {verdicts['ERROR']:,}",
        "",
    ]
    if fail_reasons:
        lines.append("failing metric (assemblies may fail more than one)")
        for m, k in sorted(fail_reasons.items(), key=lambda kv: -kv[1]):
            lines.append(f"  {m:22s} {k:,}")
        lines.append("")
    if warn_reasons:
        lines.append("warning metric")
        for m, k in sorted(warn_reasons.items(), key=lambda kv: -kv[1]):
            lines.append(f"  {m:22s} {k:,}")
        lines.append("")
    if skipped:
        lines += [
            "NOT EVALUATED — these thresholded metrics had no input:",
            *(f"  {m:22s} requires {thr[m]['needs_tool']}" for m in skipped),
            "",
            "Every verdict above is therefore PARTIAL. Supply the missing",
            "inputs before describing this as a complete QualiBact assessment.",
            "",
        ]
    lines += [
        f"metrics table : {metrics_path}",
        ("quarantined   : %s to %s" % (f"{quarantined:,}", a.quarantine))
        if a.apply else
        "DRY RUN — nothing moved. Re-run with --apply to quarantine failures.",
    ]
    text = "\n".join(lines)
    summary_path.write_text(text + "\n")
    print("\n" + text)
    print(f"\nsummary       : {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
