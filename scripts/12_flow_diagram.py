#!/usr/bin/env python
"""12_flow_diagram.py — Figure 1: study flow diagram.

Vertical CONSORT-style flow: retrieval -> QC -> genotyping -> phenotype
linkage -> stratification -> classification -> analysis sets.
Exclusions branch to the right at each stage.

Usage:
    python scripts/12_flow_diagram.py
"""
import argparse
from pathlib import Path

import sys
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parent))
import figstyle
from matplotlib.patches import FancyBboxPatch

FIG = Path("results/figures")
TAB = Path("results/tables")


# Assembly-flow counts are read from a small audit file so the figure can
# never drift from the analysis. Regenerate it with 08b (it prints these
# numbers) or edit config/flow_counts.yaml.
# Assembly-retrieval counts precede the record-level analysis and are recorded
# once in config/flow_counts.yaml. Every downstream count is read from the
# audit TSV written by 08b, so the figure cannot drift from the analysis.
FLOW_FILE = Path("config/flow_counts.yaml")


def _upstream():
    import yaml
    if FLOW_FILE.exists():
        return yaml.safe_load(FLOW_FILE.read_text())
    sys.exit(f"\nMissing {FLOW_FILE}. It records the assembly-retrieval counts "
             f"that precede the record-level analysis.\n")


def _audit(drug, tag):
    f = TAB / f"recordlevel_audit_{drug}{tag}.tsv"
    if not f.exists():
        sys.exit(f"\nMissing {f}.\nRun 08b_discordance_recordlevel.py first; "
                 f"Figure 1 counts are read from its audit output.\n")
    d = pd.read_csv(f, sep="\t")
    return dict(zip(d["stage"], d["n"]))


def counts_from_tables(tag="_clsi_bmd"):
    """Read the flow counts from generated tables. Hard-coding them previously
    left the figure showing superseded v1 values."""
    f = TAB / f"recordlevel_discordance_meropenem{tag}.tsv"
    if not f.exists():
        sys.exit(f"\nMissing {f}.\nRun 08b_discordance_recordlevel.py first; "
                 f"figure counts are read from its output, never hard-coded.\n")
    d = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
    pos = d["carb_pos"].astype(str).str.lower().isin({"true", "1"})
    binary = d["call"].isin(["S", "R"])
    c = dict(
        joined=len(d),
        binary=int(binary.sum()),
        intermediate=int((d["call"] == "I").sum()),
        ambiguous=int((d["call"] == "ambiguous").sum()),
        carb_pos=int((pos & binary).sum()),
        carb_pos_susc=int((pos & binary & (d["call"] == "S")).sum()),
        carb_neg=int((~pos & binary).sum()),
        carb_neg_res=int((~pos & binary & (d["call"] == "R")).sum()),
    )
    fi = TAB / f"recordlevel_discordance_imipenem{tag}.tsv"
    if fi.exists():
        di = pd.read_csv(fi, sep="\t", dtype=str, low_memory=False)
        c["imi_joined"] = len(di)
        c["imi_binary"] = int(di["call"].isin(["S", "R"]).sum())
    else:
        c["imi_joined"] = c["imi_binary"] = 0
    up = _upstream()
    am = _audit("meropenem", tag)
    ai = _audit("imipenem", tag) if (TAB / f"recordlevel_audit_imipenem{tag}.tsv").exists() else {}
    c.update(up)
    c["records_mero"] = am.get("phenotype records for this drug", 0)
    c["numeric_mero"] = am.get("records with a numeric value "
                               "(mm and computational excluded)", 0)
    c["std_mero"] = [v for k, v in am.items() if k.startswith("records interpreted")][0]
    c["meth_mero"] = [v for k, v in am.items() if k.startswith("records generated")][0]
    c["typed"] = am.get("Kleborate genomes", 0)
    c["confirmed_kp"] = am.get("confirmed K. pneumoniae", 0)
    c["kleborate_failed"] = am.get("genotyping failures", 0)
    c["non_kp"] = c["typed"] - c["confirmed_kp"]
    c["linked_mero"] = am.get("genomes with linked genotype", 0)
    # Assemblies dropped by QualiBact. The audit's "genomes with linked
    # genotype" row counts what survived QC, so the pre-QC figure is the sum.
    c["qc_failed_asm"] = am.get(
        "failed assembly quality control (QualiBact v1.0)", 0)
    c["linked_mero_preqc"] = c["linked_mero"] + c["qc_failed_asm"]
    c["linked_imi"] = ai.get("genomes with linked genotype", 0)
    c["linked_imi_preqc"] = c["linked_imi"] + ai.get(
        "failed assembly quality control (QualiBact v1.0)", 0)
    # The assembly-QC filter is applied to the genotype table by 05b before
    # run_all.sh, so the audit above describes only QC-retained genomes: it
    # reports every pre-QC exclusion as a genotyping failure and the QC step
    # as zero. Where config/flow_counts.yaml records the pre-QC stages, they
    # take precedence and the exclusions are recomputed as stage differences.
    # The record total was hand-entered in the YAML and drifted: it read 11,538
    # against the 7,289 + 3,925 = 11,214 the audit tables report for the two
    # carbapenems this figure describes. It is computed here instead.
    ai_recs = ai.get("phenotype records for this drug", 0)
    if c["records_mero"] and ai_recs:
        c["records"] = c["records_mero"] + ai_recs
    if "typed" in up:
        c["typed"] = up["typed"]
        c["confirmed_kp"] = up.get("confirmed_kp", c["confirmed_kp"])
        c["kleborate_failed"] = c["submitted"] - c["typed"]
        c["non_kp"] = c["typed"] - c["confirmed_kp"]
    for drug, key in (("mero", "linked_mero_preqc"), ("imi", "linked_imi_preqc")):
        if key in up:
            c[key] = up[key]
    c["qc_failed_asm"] = c["linked_mero_preqc"] - c["linked_mero"]
    c["binary_mero"] = am.get("binary analysis set", 0)
    c["binary_imi"] = ai.get("binary analysis set", 0)
    c["unres_mero"] = am.get("unresolvable", 0)
    c["inter_mero"] = am.get("intermediate", 0)
    return c

figstyle.apply_style()

INK   = "#12263A"
MAIN  = "#1B5E8C"
MAINF = "#EAF2F8"
EXC   = "#9B3B2E"
EXCF  = "#FBEFEC"
OUT   = "#1E6F52"
OUTF  = "#E9F5EF"
RULE  = "#B9C7D2"

# main-column steps: (title, detail lines, n-label)
STEPS = [
    ("BV-BRC retrieval",
     ["$\\it{Klebsiella\\ pneumoniae}$ isolates with laboratory-measured",
      "carbapenem susceptibility (computational predictions excluded)"],
     "__RETRIEVAL__"),
    ("Assembly retrieval",
     ["genome assemblies downloaded over FTPS"],
     "__ASSEMBLIES__"),
    # This step checks that a file is a parseable FASTA, because the genotyper
    # aborts a whole batch on a malformed record. It is NOT quality control,
    # and calling it that was the defect the reviewer identified. Assembly
    # quality is assessed separately, below.
    ("Structural FASTA check",
     ["valid multi-FASTA, ≥1 kb, no zero-length sequence records"],
     "__QC__"),
    ("Genotyping — Kleborate v3.2.4 (kpsc preset)",
     ["sequence type · acquired carbapenemase alleles",
      "OmpK35/OmpK36 mutations · ESBL determinants"],
     "__TYPED__"),
    ("Genotype–phenotype linkage",
     ["genomes with a linked laboratory-measured MIC"],
     "__LINKAGE__"),
    ("Assembly quality control — QualiBact v1.0",
     ["genome size · GC · N50 · contigs · completeness",
      "contamination · coding sequences (species-specific bounds)"],
     "__QCPASS__"),
    ("Methodological stratification",
     ["CLSI interpretive standard · broth dilution only",
      "(primary analysis; unstratified data reported as supplementary)"],
     "__STRAT2__"),
    ("Censoring-aware S/R classification",
     ["interval-censored MICs resolved where the censoring bound",
      "lies wholly on one side of the breakpoint"],
     "__CLASS2__"),
]

# Exclusions aligned to the gap BELOW each step (index -> (title, detail, count)).
#
# The detail and count fields are passed through the substitution table, so they
# must name a computed key rather than carry a literal. Earlier versions wrote
# the numbers in by hand and the figure went on displaying superseded v1 counts
# after the analysis was rerun; the whole point of counts_from_tables() is that
# no number in this figure is typed twice.
#
# Index 4 is the one remaining literal. The count of records dropped for using a
# non-CLSI standard or a method other than broth dilution is not recoverable
# from the audit keys this script already reads, so it needs a new key emitted
# by 08b_discordance_recordlevel.py. Until that exists the value is marked and
# must be checked against the audit table by hand.
EXCLUS = {
    0: ("Assembly not retrievable", "no sequence deposited at the\nBV-BRC path",
        "__EXC1N__"),
    1: ("Failed structural FASTA check",
        "zero-length sequence records\n(abort genotyping)", "__EXC2N__"),
    2: ("Genotyping failure or", "__EXC3__", "__EXC3N__"),
    3: ("No linked meropenem MIC",
        "genome typed but no laboratory\nMIC available for this drug",
        "__EXC4N__"),
    4: ("Failed assembly quality control",
        "QualiBact v1.0 FAIL on any\nmetric, or no metric computable",
        "__EXCQCN__"),
    5: ("Non-CLSI standard or", "method other than broth\ndilution",
        "__EXC5N__"),
    6: ("Censoring bound straddles", "breakpoint — S/R not\nresolvable",
        "__EXC6N__"),
}

TERMINAL = [
    ("Carbapenemase-positive", "__CP_N__", "__CP_S__"),
    ("Carbapenemase-negative", "__CN_N__", "__CN_R__"),
]


def rbox(ax, cx, cy, w, h, ec, fc, lw=1.4):
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                boxstyle="round,pad=0.008,rounding_size=0.012",
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2))


def vline(ax, x, y0, y1, color=MAIN, lw=1.5):
    ax.plot([x, x], [y0, y1], color=color, lw=lw, zorder=1, solid_capstyle="round")


def arrowhead(ax, x, y, color=MAIN, size=0.010):
    ax.plot([x - size, x, x + size], [y + 1.8 * size, y, y + 1.8 * size],
            color=color, lw=1.5, zorder=1, solid_joinstyle="miter")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="Figure1_study_flow")
    a = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    c = counts_from_tables()
    subs = {
        "__RETRIEVAL__": (f"{c['records']:,} phenotype records from "
                          f"{c['genomes_with_pheno']:,} genomes"),
        "__ASSEMBLIES__": f"{c['assemblies']:,} assemblies",
        "__QC__": f"{c['submitted']:,} assemblies",
        "__TYPED__": (f"{c['typed']:,} typed; {c['confirmed_kp']:,} confirmed "
                      f"$\\it{{K.\\,pneumoniae}}$"),
        "__LINKAGE2__": (f"meropenem {c['linked_mero']:,}   ·   "
                         f"imipenem {c['linked_imi']:,}"),
        "__STRAT2__": (f"meropenem {c['meth_mero']:,} records   ·   "
                       f"CLSI {c['std_mero']:,} of {c['numeric_mero']:,}"),
        "__CLASS2__": (f"binary {c['binary_mero']:,}  ·  intermediate "
                       f"{c['inter_mero']:,}  ·  unresolved {c['unres_mero']:,}"),
        # Linkage precedes assembly QC in the flow, so this row shows the
        # pre-QC count and the QC step below it shows what survived.
        "__LINKAGE__": (f"meropenem {c['linked_mero_preqc']:,}   ·   "
                        f"imipenem {c['linked_imi_preqc']:,}"),
        "__QCPASS__": (f"meropenem {c['linked_mero']:,}   ·   "
                       f"imipenem {c['linked_imi']:,}"),
        "__EXCQCN__": f"n = {c['qc_failed_asm']:,}",
        "__STRAT__":   f"meropenem {c['meth_mero']:,} records   ·   imipenem records",
        "__CLASSIFIED__": (f"binary {c['binary_mero']:,}  ·  intermediate "
                           f"{c['inter_mero']:,}  ·  unresolved {c['unres_mero']:,}"),
        "__EXC3__": (f"Kleborate failure (n = {c['kleborate_failed']}) and\n"
                     f"species mismatch (n = {c['non_kp']})"),
        "__EXC3N__": f"n = {c['kleborate_failed'] + c['non_kp']}",
        # Every exclusion count below is a difference between two stages that
        # are themselves read from the generated tables, so the arithmetic in
        # the figure cannot disagree with Table 1.
        "__EXC1N__": f"n = {c['genomes_with_pheno'] - c['assemblies']:,}",
        "__EXC2N__": f"n = {c['assemblies'] - c['submitted']:,}",
        "__EXC4N__": f"n = {c['confirmed_kp'] - c['linked_mero_preqc']:,}",
        "__EXC5N__": f"n = {c['numeric_mero'] - c['meth_mero']:,}",
        "__EXC6N__": f"n = {c['unres_mero']:,}",
        "__CP_N__": f"n = {c['carb_pos']:,}",
        "__CP_S__": (f"{c['carb_pos_susc']:,} "
                     f"({100*c['carb_pos_susc']/max(c['carb_pos'],1):.1f}%) "
                     f"phenotypically susceptible"),
        "__CN_N__": f"n = {c['carb_neg']:,}",
        "__CN_R__": (f"{c['carb_neg_res']:,} "
                     f"({100*c['carb_neg_res']/max(c['carb_neg'],1):.1f}%) "
                     f"phenotypically resistant"),
    }
    steps = [(a_, b_, subs.get(n_, n_)) for a_, b_, n_ in STEPS]
    terminal = [(a_, subs.get(b_, b_), subs.get(c_, c_)) for a_, b_, c_ in TERMINAL]

    n = len(steps)
    fig, ax = plt.subplots(figsize=(10.6, 12.8))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    CX, W, H = 0.365, 0.60, 0.088
    EX, EW, EH = 0.845, 0.295, 0.088
    top, bottom = 0.955, 0.175
    gap = (top - bottom) / (n - 1)

    ys = [top - i * gap for i in range(n)]

    for i, ((title, detail, nlab), y) in enumerate(zip(steps, ys)):
        last = (i == n - 1)
        ec, fc = (OUT, OUTF) if last else (MAIN, MAINF)
        rbox(ax, CX, y, W, H, ec, fc)
        ax.text(CX, y + H / 2 - 0.017, title, ha="center", va="top",
                fontsize=10.2, fontweight="bold", color=ec, zorder=3)
        ax.text(CX, y + H / 2 - 0.036, "\n".join(detail), ha="center", va="top",
                fontsize=8.3, color="#41525F", zorder=3, linespacing=1.45)
        # The count sits below the detail block, not at a fixed offset from the
        # box floor: a two-line detail used to run straight through it.
        ax.text(CX, y + H / 2 - 0.040 - 0.0145 * len(detail), nlab,
                ha="center", va="top",
                fontsize=9.3, fontweight="bold", color=INK, zorder=3)

        if not last:
            y_next_top = ys[i + 1] + H / 2
            y_bot = y - H / 2
            mid = (y_bot + y_next_top) / 2
            vline(ax, CX, y_bot, y_next_top + 0.012)
            arrowhead(ax, CX, y_next_top + 0.002)
            if i in EXCLUS:
                t1, t2, nn = EXCLUS[i]
                t2 = subs.get(t2, t2); nn = subs.get(nn, nn)
                ax.plot([CX, EX - EW / 2], [mid, mid], color=RULE, lw=1.2,
                        zorder=1, linestyle=(0, (4, 2)))
                rbox(ax, EX, mid, EW, EH, EXC, EXCF, lw=1.2)
                ax.text(EX, mid + EH / 2 - 0.012, t1, ha="center", va="top",
                        fontsize=8.5, fontweight="bold", color=EXC, zorder=3)
                ax.text(EX, mid + EH / 2 - 0.031, t2, ha="center", va="top",
                        fontsize=7.7, color="#6B4A43", zorder=3, linespacing=1.4)
                ax.text(EX, mid - EH / 2 + 0.009, nn, ha="center", va="bottom",
                        fontsize=8.6, fontweight="bold", color=EXC, zorder=3)

    # terminal analysis sets
    y_last = ys[-1] - H / 2
    y_term = 0.052
    xs = [0.175, 0.555]
    vline(ax, CX, y_last, y_last - 0.036, color=OUT)
    ax.plot([xs[0], xs[1]], [y_last - 0.036] * 2, color=OUT, lw=1.5, zorder=1)
    for cx, (t, nlab, sub) in zip(xs, terminal):
        vline(ax, cx, y_last - 0.036, y_term + 0.050, color=OUT)
        arrowhead(ax, cx, y_term + 0.040, color=OUT)
        rbox(ax, cx, y_term, 0.335, 0.086, OUT, OUTF)
        ax.text(cx, y_term + 0.030, t, ha="center", va="center",
                fontsize=9.8, fontweight="bold", color=OUT, zorder=3)
        ax.text(cx, y_term + 0.006, nlab, ha="center", va="center",
                fontsize=10.4, fontweight="bold", color=INK, zorder=3)
        ax.text(cx, y_term - 0.022, sub, ha="center", va="center",
                fontsize=8.0, color="#41525F", zorder=3)

    ax.text(0.5, y_term - 0.062,
            "Primary analysis set — meropenem, CLSI broth dilution",
            ha="center", va="center", fontsize=8.6, style="italic", color="#6A7A87")

    out = FIG / a.out
    figstyle.save(fig, out)
    print(f"[out] {out}.png / .pdf / .svg")


if __name__ == "__main__":
    main()
