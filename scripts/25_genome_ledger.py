#!/usr/bin/env python3
"""25_genome_ledger.py — one row per genome, from retrieval to inclusion.

WHY THIS EXISTS
The repair of the truncated downloads is documented across several logs and
tables, and their totals nearly but not exactly agree: 3,493 assemblies were
re-retrieved and verified, 3,493 were submitted to Kleborate, one failed, and
3,491 reached the merged genotype table. Totals that almost reconcile are not
evidence; they are a question. This script answers it at genome-ID level and
states where each genome left the study, so every count the paper reports can
be traced to a set of identifiers rather than to a log line.

It is a reporting step. It computes nothing the analysis depends on, and it
fails loudly if the stages disagree in a way the paper does not disclose.

COLUMNS
    genome_id
    source_dir          which retrieval the genotyped file came from
    integrity           COMPLETE, UNVERIFIED or INCOMPLETE against BV-BRC
    local_size_bp / local_contigs / bvbrc_size_bp / bvbrc_contigs
    typed               Kleborate produced a row for this genome
    species             Kleborate's species call
    species_retained    the species name contains "pneumoniae"
    qc_verdict          QualiBact PASS, WARN, FAIL, or none if no metrics
    qc_retained         PASS or WARN
    meropenem_linked / imipenem_linked
    meropenem_call / imipenem_call      S, R, or blank when not in the binary set
    exit_stage          the first stage that removed this genome, or "included"

USAGE (from the project root)
    python3 scripts/25_genome_ledger.py

OUTPUT
    results/tables/genome_ledger.tsv
    results/tables/genome_ledger_summary.txt
"""
import csv
import os
import re
import sys
from collections import Counter, OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPECIES_COL = "enterobacterales__species__species"


def path(*parts):
    return os.path.join(ROOT, *parts)


def die(msg):
    sys.exit(f"25_genome_ledger: {msg}")


def read_tsv(p, key=None, required=True):
    if not os.path.exists(p):
        if required:
            die(f"missing {p}")
        return {} if key else []
    with open(p, newline="", encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r[key]: r for r in rows} if key else rows


def main():
    integrity = read_tsv(path("results", "qc_concordance",
                              "fasta_integrity_effective.tsv"), key="genome_id")
    typed = read_tsv(path("data", "processed",
                          "kleborate_klebsiella_pneumo_complex_output.prefilter.tsv"),
                     key="strain")
    qc = read_tsv(path("data", "interim", "assembly_qc_metrics.tsv"),
                  key="sample", required=False)
    # genomes Kleborate assigned to a species group outside the complex: they
    # are typed, but they leave the complex table, so without this they would
    # be indistinguishable from a genotyping failure
    other_group = {}
    for fname, label in (("kleborate_escherichia_output.tsv", "Escherichia"),
                         ("kleborate_klebsiella_oxytoca_complex_output.tsv",
                          "K. oxytoca complex")):
        for r in read_tsv(path("data", "processed", fname), required=False):
            gid = (r.get("strain") or "").strip()
            if gid:
                other_group[gid] = label

    # the retrieval attempt log for the repaired assemblies
    manifest = {}
    for r in read_tsv(path("data", "raw", "genomes_redownload", "manifest.tsv"),
                      required=False):
        gid = r.get("genome_id")
        if not gid:
            continue
        prev = manifest.get(gid)
        if prev is None or r.get("status") == "ok":
            manifest[gid] = r

    calls = {}
    for drug in ("meropenem", "imipenem"):
        p = path("results", "tables",
                 f"recordlevel_discordance_{drug}_clsi_bmd.tsv")
        calls[drug] = read_tsv(p, key="genome_id")

    ids = sorted(set(integrity) | set(typed)
                 | set(calls["meropenem"]) | set(calls["imipenem"]))
    if not ids:
        die("no genomes found in any input")

    out = path("results", "tables", "genome_ledger.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cols = ["genome_id", "source_dir", "redownload_status", "redownload_attempts",
            "integrity", "local_size_bp", "local_contigs", "bvbrc_size_bp",
            "bvbrc_contigs", "typed", "species_group", "species",
            "species_retained", "qc_verdict", "qc_retained",
            "meropenem_linked", "meropenem_call", "imipenem_linked",
            "imipenem_call", "exit_stage"]
    attempts = Counter(r.get("genome_id") for r in read_tsv(
        path("data", "raw", "genomes_redownload", "manifest.tsv"),
        required=False))
    stages = Counter()
    rows = []
    for gid in ids:
        ig = integrity.get(gid, {})
        tp = typed.get(gid)
        q = qc.get(gid, {})
        species = (tp or {}).get(SPECIES_COL, "") if tp else ""
        retained = bool(species) and "pneumoniae" in species
        verdict = (q.get("verdict") or "").strip().upper()
        qc_ret = verdict in ("PASS", "WARN")
        row = OrderedDict((c, "") for c in cols)
        row.update(
            genome_id=gid,
            source_dir=ig.get("source_dir", ""),
            integrity=ig.get("verdict", "NOT_RETRIEVED"),
            local_size_bp=ig.get("local_size_bp", ""),
            local_contigs=ig.get("local_contigs", ""),
            bvbrc_size_bp=ig.get("bvbrc_size_bp", ""),
            bvbrc_contigs=ig.get("bvbrc_contigs", ""),
            redownload_status=(manifest.get(gid) or {}).get("status", ""),
            redownload_attempts=attempts.get(gid, "") or "",
            typed="yes" if tp else "no",
            species_group=("K. pneumoniae complex" if tp
                           else other_group.get(gid, "")),
            species=species,
            species_retained="yes" if retained else "no",
            qc_verdict=verdict or "none",
            qc_retained="yes" if qc_ret else "no",
        )
        for drug, pre in (("meropenem", "meropenem"), ("imipenem", "imipenem")):
            c = calls[drug].get(gid)
            row[f"{pre}_linked"] = "yes" if c else "no"
            row[f"{pre}_call"] = (c or {}).get("call", "")

        if not re.match(r"^\d+\.\d+$", gid):
            # three rows in the Kleborate output carry no genome identifier.
            # They are fragments of neighbouring records left by a newline
            # inside a field, and counting them inflates both the typed total
            # and the species exclusions.
            stage = "malformed Kleborate row, no genome identifier"
        elif not tp and gid in other_group:
            stage = f"assigned to {other_group[gid]}, outside the complex"
        elif not tp:
            stage = ("not retrieved" if row["integrity"] == "NOT_RETRIEVED"
                     else "no Kleborate call")
        elif not retained:
            stage = "species not in the K. pneumoniae complex"
        elif not qc_ret:
            stage = ("failed assembly QC" if verdict == "FAIL"
                     else "no QC metrics available")
        elif row["meropenem_call"] in ("S", "R"):
            stage = "included"
        elif row["meropenem_linked"] == "yes":
            stage = "linked but not binary-classifiable"
        else:
            stage = "no CLSI broth-dilution meropenem record"
        row["exit_stage"] = stage
        stages[stage] += 1
        rows.append(row)

    with open(out, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)

    # reconciliation against the counts the paper and the figure state
    n_typed = sum(1 for r in rows if r["typed"] == "yes"
                  and re.match(r"^\d+\.\d+$", r["genome_id"]))
    n_retained = sum(1 for r in rows if r["species_retained"] == "yes")
    n_qc = sum(1 for r in rows if r["species_retained"] == "yes"
               and r["qc_retained"] == "yes")
    n_inc = stages["included"]
    redown = sum(1 for r in rows
                 if r["source_dir"] == "data/raw/genomes_redownload")
    redown_typed = sum(1 for r in rows
                       if r["source_dir"] == "data/raw/genomes_redownload"
                       and r["typed"] == "yes")
    redown_other = sum(1 for r in rows
                       if r["source_dir"] == "data/raw/genomes_redownload"
                       and r["typed"] == "no" and r["species_group"])
    lines = [
        "genome ledger — where every genome left the study",
        "",
        f"  rows in the ledger                      {len(rows):>7,}",
        f"  re-retrieved after the integrity audit  {redown:>7,}",
        f"    of those, in the complex genotype table {redown_typed:>5,}",
        f"    assigned to another species group     {redown_other:>7,}",
        f"    unaccounted                           "
        f"{redown - redown_typed - redown_other:>7,}",
        f"  typed, in the K. pneumoniae complex      {n_typed:>7,}",
        f"  species retained (complex)              {n_retained:>7,}",
        f"  retained after assembly QC              {n_qc:>7,}",
        f"  in the primary meropenem binary set     {n_inc:>7,}",
        "",
        "  where the rest left:",
    ]
    for stage, n in stages.most_common():
        if stage != "included":
            lines.append(f"    {stage:<44} {n:>7,}")
    lines += ["", "  integrity of the genotyped file:"]
    for v, n in Counter(r["integrity"] for r in rows).most_common():
        lines.append(f"    {v:<44} {n:>7,}")
    inc_bad = [r["genome_id"] for r in rows
               if r["exit_stage"] == "included" and r["integrity"] != "COMPLETE"]
    lines += ["",
              f"  included genomes whose file is not verified COMPLETE: "
              f"{len(inc_bad):,}"]
    summary = "\n".join(lines) + "\n"
    with open(path("results", "tables", "genome_ledger_summary.txt"), "w") as fh:
        fh.write(summary)
    print(summary)
    print(f"  wrote {out}")

    if inc_bad:
        die(f"{len(inc_bad)} analysed genomes were called from a file that is "
            f"not verified complete, e.g. {inc_bad[:5]}")


if __name__ == "__main__":
    main()
