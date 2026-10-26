#!/usr/bin/env python3
"""12b_flow_diagram_typeset.py — the typeset study-flow figure (Figure 1).

WHY A SECOND FLOW SCRIPT
12_flow_diagram.py draws the same flow with matplotlib. That version is the
machine-readable one: it is laid out by the plotting library, and because
matplotlib writes italic spans as one <text> element per character, its SVG
cannot be opened and corrected by hand. Figure 1 is the one figure in the paper
whose layout was settled by eye, so the submitted version is drawn from a
typeset template instead, where every label is a single text node in Times New
Roman and the boxes sit where an author placed them.

The template carries the layout and {{PLACEHOLDER}} marks where numbers go.
This script fills them from the same tables the rest of the paper is written
from, so the figure cannot drift away from the analysis the way a hand-edited
figure does. It refuses to write the figure if any placeholder is left over or
if an internal sum disagrees.

USAGE (from the project root)
    python3 scripts/12b_flow_diagram_typeset.py

INPUTS
    assets/Figure1_study_flow_template.svg
    config/flow_counts.yaml                      pre-genotyping stages
    data/interim/assembly_qc_metrics.tsv         QualiBact verdicts
    results/tables/recordlevel_audit_{drug}_clsi_bmd.tsv
    results/tables/recordlevel_discordance_{drug}_clsi_bmd.tsv

OUTPUT
    results/figures/Figure1_study_flow.svg
    results/figures/Figure1_study_flow.png / .pdf / .tif  (if a renderer exists)
"""
import csv
import os
import re
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(ROOT, "assets", "Figure1_study_flow_template.svg")
OUTDIR = os.path.join(ROOT, "results", "figures")
STEM = "Figure1_study_flow"


def die(msg):
    sys.exit(f"12b_flow_diagram_typeset: {msg}")


def read_yaml_ints(path):
    """flow_counts.yaml is flat 'key: int' with comments. No PyYAML needed."""
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line or ":" not in line:
                continue
            k, v = line.split(":", 1)
            try:
                out[k.strip()] = int(v.strip())
            except ValueError:
                pass
    return out


def read_audit(drug):
    path = os.path.join(ROOT, "results", "tables",
                        f"recordlevel_audit_{drug}_clsi_bmd.tsv")
    if not os.path.exists(path):
        die(f"missing {path} — run 08b_discordance_recordlevel.py first")
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["stage"]] = int(row["n"])
    return out


def read_discordance(drug):
    """(carb+ n, carb+ susceptible, carb- n, carb- resistant) on the binary set."""
    path = os.path.join(ROOT, "results", "tables",
                        f"recordlevel_discordance_{drug}_clsi_bmd.tsv")
    if not os.path.exists(path):
        die(f"missing {path}")
    pos = pos_s = neg = neg_r = 0
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            call = row["call"]
            if call not in ("S", "R"):
                continue
            if row["carb_pos"].strip().lower() in ("true", "1", "yes"):
                pos += 1
                pos_s += call == "S"
            else:
                neg += 1
                neg_r += call == "R"
    return pos, pos_s, neg, neg_r


def qc_verdicts():
    path = os.path.join(ROOT, "data", "interim", "assembly_qc_metrics.tsv")
    if not os.path.exists(path):
        die(f"missing {path} — run 00e_assembly_qc.py first")
    counts = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            v = (row.get("verdict") or "").strip().upper()
            counts[v] = counts.get(v, 0) + 1
    return counts


def genotype_table():
    """Counts read from the genotype table itself, not from a run log.

    Three rows in the Kleborate output carry no genome identifier. They are
    fragments of neighbouring records, left by embedded newlines in a field,
    and they die at the species filter because their species cell is not a
    species name. Counting raw rows therefore overstates both the number of
    genomes typed and the number excluded for species, which is how an earlier
    revision of this figure reported 7,336 typed and 91 species exclusions
    instead of 7,333 and 88. Rows are counted here by identifier.

    Returns (typed_in_complex, retained, strict, quasi, malformed).
    """
    path = os.path.join(ROOT, "data", "processed",
                        "kleborate_klebsiella_pneumo_complex_output.prefilter.tsv")
    if not os.path.exists(path):
        die(f"missing {path} — needed for the genotyping counts in Figure 1")
    col = "enterobacterales__species__species"
    typed = strict = quasi = retained = malformed = 0
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gid = (row.get("strain") or "").strip()
            if not re.match(r"^\d+\.\d+$", gid):
                malformed += 1
                continue
            typed += 1
            name = (row.get(col) or "").strip()
            if "pneumoniae" not in name:
                continue
            retained += 1
            if name == "Klebsiella pneumoniae":
                strict += 1
            elif "quasipneumoniae" in name:
                quasi += 1
            else:
                die(f"{gid} is retained but is neither K. pneumoniae nor "
                    f"K. quasipneumoniae ({name}); the two-way split is wrong")
    return typed, retained, strict, quasi, malformed


def main():
    if not os.path.exists(TEMPLATE):
        die(f"missing template {TEMPLATE}")
    flow = read_yaml_ints(os.path.join(ROOT, "config", "flow_counts.yaml"))
    for k in ("genomes_with_pheno", "assemblies", "submitted", "typed",
              "confirmed_kp"):
        if k not in flow:
            die(f"config/flow_counts.yaml has no '{k}'")

    audit = {d: read_audit(d) for d in ("meropenem", "imipenem")}
    qc = qc_verdicts()
    with_metrics = sum(qc.values())
    retained = qc.get("PASS", 0) + qc.get("WARN", 0)
    failed = qc.get("FAIL", 0)
    # Quality control runs on the species-confirmed genomes, not on everything
    # Kleborate typed: the species filter has already removed 91 genomes by
    # then, and no QC metrics were ever fetched for them. Counting the QC
    # exclusions against the typed total lands on the same retained figure, but
    # it reports those 91 a second time as "no computable metrics", so the
    # denominator here is the species-confirmed set.
    no_metrics = flow["confirmed_kp"] - with_metrics
    if no_metrics < 0:
        die("more QC rows than species-confirmed genomes — inputs disagree")

    v = {
        "RECORDS": sum(a["phenotype records for this drug"] for a in audit.values()),
        "GENOMES_PHENO": flow["genomes_with_pheno"],
        "NOT_RETRIEVABLE": flow["genomes_with_pheno"] - flow["assemblies"],
        "ASSEMBLIES": flow["assemblies"],
        "QC_FAILED": flow["assemblies"] - flow["submitted"],
        "SUBMITTED": flow["submitted"],
        "KLEB_FAIL": flow["submitted"] - flow["typed"],
        "SPECIES_MISMATCH": flow["typed"] - flow["confirmed_kp"],
        "TYPING_EXCL": flow["submitted"] - flow["confirmed_kp"],
        "TYPED": flow["typed"],
        "CONFIRMED_KP": flow["confirmed_kp"],
        "QC_PASS": retained,
        # The retained set is PASS plus WARN. Reporting it as "pass" reads as
        # though every retained assembly met every threshold, which is the
        # opposite of what the paper argues about WARN retention.
        "QC_PASS_STRICT": qc.get("PASS", 0),
        "QC_WARN": qc.get("WARN", 0),
        "QC_DENOM": flow["confirmed_kp"],
        "QC_THRESHOLD_FAIL": failed,
        "QC_NO_METRICS": no_metrics,
        "QC_EXCLUDED": failed + no_metrics,
    }
    typed_cx, retained, strict, quasi, malformed = genotype_table()
    if retained != flow["confirmed_kp"]:
        die(f"genotype table has {retained} genomes in the complex but "
            f"config/flow_counts.yaml says confirmed_kp={flow['confirmed_kp']}")
    v["TYPED"] = typed_cx
    v["KP_STRICT"], v["KQ"] = strict, quasi
    v["MALFORMED_ROWS"] = malformed
    # Genomes Kleborate assigned to another species group are typed, but they
    # leave the K. pneumoniae complex table, so they are named separately.
    v["OTHER_GROUP"] = flow.get("other_species_group", 0)
    v["CALLED"] = typed_cx + v["OTHER_GROUP"]
    v["SPECIES_MISMATCH"] = typed_cx - retained
    v["KLEB_FAIL"] = flow["submitted"] - v["CALLED"]
    v["TYPING_EXCL"] = flow["submitted"] - retained
    numeric_key = "records with a numeric value (mm and computational excluded)"
    for drug, pre in (("meropenem", "M"), ("imipenem", "I")):
        a = audit[drug]
        pos, pos_s, neg, neg_r = read_discordance(drug)
        v.update({
            f"{pre}_RECORDS": a["phenotype records for this drug"],
            f"{pre}_NUMERIC": a[numeric_key],
            f"{pre}_CLSI": a["records interpreted against CLSI"],
            f"{pre}_BMD": a["records generated by Broth dilution"],
            f"{pre}_NONCLSI": a[numeric_key] - a["records interpreted against CLSI"],
            f"{pre}_OTHERMETHOD": (a["records interpreted against CLSI"]
                                   - a["records generated by Broth dilution"]),
            f"{pre}_LINKED": a["genomes with linked genotype"],
            f"{pre}_S": a["susceptible"],
            f"{pre}_I": a["intermediate"],
            f"{pre}_R": a["resistant"],
            f"{pre}_U": a["unresolvable"],
            f"{pre}_BINARY": a["binary analysis set"],
            f"{pre}_POS": pos, f"{pre}_POS_S": pos_s,
            f"{pre}_NEG": neg, f"{pre}_NEG_R": neg_r,
            f"{pre}_POS_S_PCT": f"{100 * pos_s / pos:.1f}%" if pos else "n/a",
            f"{pre}_NEG_R_PCT": f"{100 * neg_r / neg:.1f}%" if neg else "n/a",
        })
        # The figure claims the two outcome boxes partition the binary set.
        if pos + neg != a["binary analysis set"]:
            die(f"{drug}: carbapenemase-positive {pos} + negative {neg} "
                f"!= binary analysis set {a['binary analysis set']}")
        if a["susceptible"] + a["resistant"] != a["binary analysis set"]:
            die(f"{drug}: S + R != binary analysis set")
    if v["QC_PASS"] + v["QC_EXCLUDED"] != v["CONFIRMED_KP"]:
        die("QC retained + excluded != species-confirmed genomes")

    def fmt(x):
        return f"{x:,}" if isinstance(x, int) else str(x)

    svg = open(TEMPLATE, encoding="utf-8").read()
    used = set()

    def sub(m):
        used.add(m.group(1))
        if m.group(1) not in v:
            die(f"template placeholder {{{{{m.group(1)}}}}} has no value")
        return fmt(v[m.group(1)])

    svg = re.sub(r"\{\{(\w+)\}\}", sub, svg)
    left = re.findall(r"\{\{(\w+)\}\}", svg)
    if left:
        die(f"placeholders left unfilled: {sorted(set(left))}")
    unused = sorted(set(v) - used)
    if unused:
        print(f"  note: computed but not placed in the figure: {', '.join(unused)}")

    os.makedirs(OUTDIR, exist_ok=True)
    out_svg = os.path.join(OUTDIR, f"{STEM}.svg")
    with open(out_svg, "w", encoding="utf-8") as fh:
        fh.write(svg)
    print(f"  wrote {out_svg}")
    for k in sorted(v):
        print(f"    {k:<20} {fmt(v[k])}")

    # Raster and print copies, with whichever renderer is installed. The SVG
    # is the figure of record; these are conveniences, and a machine without a
    # renderer still gets the figure.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from render_svg import render
    render(out_svg, scale=4)


if __name__ == "__main__":
    main()
