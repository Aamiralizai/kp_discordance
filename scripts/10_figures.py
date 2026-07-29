#!/usr/bin/env python
"""10_figures.py — publication-quality multi-panel figures.

FIGURE 1 (main): the discordance result
   A  surveillance confusion matrix + performance metrics
   B  discordance by carbapenemase family (forest plot, Wilson CIs)
   C  mechanism: porin state within family
   D  allele-level resolution
   E  meropenem vs imipenem (drug-specific reversal)

FIGURE 2 (supp): MIC distributions by genotype group

Usage:
    python scripts/10_figures.py
    python scripts/10_figures.py --drug imipenem
"""
import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parent))
import figstyle
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

TAB = Path("results/tables")
FIG = Path("results/figures")

# --- house style: restrained, print-safe, colourblind-friendly ---
figstyle.apply_style()

C = {
    "OXA-48-like": "#E8743B",   # orange
    "KPC":         "#1F78B4",   # blue
    "NDM":         "#33A02C",   # green
    "VIM":         "#6A3D9A",   # purple
    "IMP":         "#B15928",   # brown
    "other":       "#999999",
    "grey":        "#666666",
    "light":       "#D9D9D9",
    "accent":      "#CB181D",
}
import sys
sys.path.insert(0, str(Path(__file__).parent))
from kp_defs import PORIN_ORDER, PORIN_LABEL
PORIN_C = {"none": "#8EC4E8", "ompK35_loss": "#4A90C2",
           "ompK36_reduced_expression": "#F2B880", "ompK36_loop3": "#E8743B",
           "ompK36_other_variant": "#C98B6B", "ompK36_loss": "#A8322D"}
PORIN_LBL = {k: v.replace("OmpK36 ", "OmpK36\n").replace("OmpK35 ", "OmpK35\n")
             for k, v in PORIN_LABEL.items()}


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n; z = 1.959964
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (p, max(0, c - h), min(1, c + h))


def fam_colour(f):
    s = str(f).upper()
    if any(x in s for x in ("OXA-48", "OXA-181", "OXA-232", "OXA-244", "OXA-162")):
        return C["OXA-48-like"]
    if "KPC" in s: return C["KPC"]
    if "NDM" in s: return C["NDM"]
    if "VIM" in s: return C["VIM"]
    if "IMP" in s: return C["IMP"]
    return C["other"]


def panel_label(ax, letter, dx=-0.16, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")


def load(drug, tag=""):
    """Prefer the record-level table; fall back to the legacy one."""
    f = TAB / f"recordlevel_discordance_{drug}{tag}.tsv"
    if not f.exists():
        raise SystemExit(
            f"\nMissing {f}.\nRun 08b_discordance_recordlevel.py first. "
            f"Legacy discordance_*.tsv outputs are INVALID (see "
            f"archive/deprecated/) and are no longer used as a fallback.\n")
    if not f.exists():
        raise SystemExit(f"{f} not found — run 08_discordance_analysis.py first")
    d = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
    if "mic" in d.columns:
        d["mic"] = pd.to_numeric(d["mic"], errors="coerce")
    else:
        # record-level output: interval [lo, hi]; use the finite bound for display
        lo = pd.to_numeric(d.get("lo"), errors="coerce")
        hi = pd.to_numeric(d.get("hi"), errors="coerce")
        d["mic"] = np.where(np.isfinite(hi), hi, lo)
        d["censored_display"] = ~np.isfinite(hi) | (lo == 0)
    d["carb_pos"] = d["carb_pos"].astype(str).str.lower().isin({"true", "1"})
    # Rates must be computed on the BINARY set. Including intermediate and
    # unresolved isolates in the denominator understates every percentage
    # (e.g. OXA-48-like 46/145 = 31.7% instead of 46/120 = 38.3%).
    d.attrs["n_all"] = len(d)
    d = d[d["call"].isin(["S", "R"])].copy()
    d["susceptible"] = d["call"] == "S"
    return d


# ---------------------------------------------------------------- panels
def panel_confusion(ax, d, drug):
    pos, neg = d[d["carb_pos"]], d[~d["carb_pos"]]
    tp = int((pos["call"] == "R").sum()); fp = int((pos["call"] == "S").sum())
    fn = int((neg["call"] == "R").sum()); tn = int((neg["call"] == "S").sum())
    M = np.array([[tp, fp], [fn, tn]])
    ax.imshow(np.array([[1, .35], [.35, 1]]), cmap="Blues", vmin=0, vmax=1.6, aspect="auto")
    for i in range(2):
        for j in range(2):
            disc = (i == 0 and j == 1) or (i == 1 and j == 0)
            ax.text(j, i, f"{M[i,j]:,}", ha="center", va="center",
                    fontsize=15, fontweight="bold",
                    color=C["accent"] if disc else "#123")
            if disc:
                ax.text(j, i + .30, "discordant", ha="center", va="center",
                        fontsize=7.5, color=C["accent"], style="italic")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["resistant", "susceptible"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["carbapenemase\ndetected",
                                               "no carbapenemase"])
    ax.set_xlabel(f"{drug} phenotype", labelpad=6)
    ax.set_title("Genotype vs phenotype", pad=8)
    ppv = tp / (tp + fp) if tp + fp else np.nan
    sens = tp / (tp + fn) if tp + fn else np.nan
    ax.text(0.5, -0.22, f"PPV {ppv:.1%}    sensitivity {sens:.1%}",
            transform=ax.transAxes, ha="center", fontsize=8.5, color=C["grey"],
            fontweight="bold")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)


def panel_family(ax, d, min_n=10):
    pos = d[d["carb_pos"]]
    rows = []
    for fam, g in pos.groupby("carb_family"):
        if len(g) < min_n or fam == "none":
            continue
        k = int(g["susceptible"].sum()); p, lo, hi = wilson(k, len(g))
        rows.append((fam, len(g), 100 * p, 100 * lo, 100 * hi))
    if not rows:
        ax.set_axis_off(); return
    t = pd.DataFrame(rows, columns=["fam", "n", "p", "lo", "hi"]).sort_values("p")
    y = np.arange(len(t))
    for yi, (_, r) in zip(y, t.iterrows()):
        col = fam_colour(r["fam"])
        ax.plot([r["lo"], r["hi"]], [yi, yi], color=col, lw=2.2,
                solid_capstyle="round", alpha=.85)
        ax.scatter(r["p"], yi, s=48, color=col, zorder=3,
                   edgecolor="white", linewidth=.9)
        ax.text(r["hi"] + 1.6, yi, f"n={r['n']:,}", va="center",
                fontsize=7.5, color=C["grey"])
    ax.set_yticks(y); ax.set_yticklabels(t["fam"])
    ax.set_xlabel("carbapenemase-positive isolates\nphenotypically susceptible (%)")
    ax.set_title("Discordance by enzyme family", pad=8)
    ax.set_xlim(-2, max(t["hi"]) * 1.32 + 6)
    ax.grid(axis="x", alpha=.25, lw=.6)
    ax.set_axisbelow(True)


def panel_mechanism(ax, d, fams=("OXA-48-like", "KPC"), min_n=10):
    pos = d[d["carb_pos"]]
    width = 0.8
    xt, xl, xpos = [], [], 0
    for fam in fams:
        sub = pos[pos["carb_family"] == fam]
        if len(sub) < min_n:
            continue
        start = xpos
        for pr in PORIN_ORDER:
            g = sub[sub["porin"] == pr]
            if len(g) < min_n:
                continue
            k = int(g["susceptible"].sum()); p, lo, hi = wilson(k, len(g))
            ax.bar(xpos, 100 * p, width, color=PORIN_C[pr],
                   edgecolor="white", linewidth=.8)
            ax.errorbar(xpos, 100 * p, yerr=[[100 * (p - lo)], [100 * (hi - p)]],
                        color="#333", lw=1.1, capsize=2.5, capthick=1.1)
            ax.text(xpos, 100 * hi + 1.8, f"{len(g):,}", ha="center",
                    fontsize=7, color=C["grey"])
            xpos += 1
        if xpos > start:
            xt.append((start + xpos - 1) / 2); xl.append(fam)
            xpos += 0.9
    ax.set_xticks(xt); ax.set_xticklabels(xl)
    ax.set_ylabel("phenotypically susceptible (%)")
    ax.set_title("Porin state drives discordance", pad=8)
    ax.margins(y=.16)
    ax.grid(axis="y", alpha=.25, lw=.6); ax.set_axisbelow(True)
    handles = [mpatches.Patch(facecolor=PORIN_C[p], label=PORIN_LBL[p])
               for p in PORIN_ORDER]
    ax.legend(handles=handles, frameon=False, loc="upper right",
              handlelength=1.1, borderpad=.3, labelspacing=.35)


def panel_allele(ax, d, min_n=10, top=9):
    pos = d[d["carb_pos"]]
    rows = []
    for al, g in pos.groupby("carb_allele"):
        if len(g) < min_n or ";" in str(al):
            continue
        k = int(g["susceptible"].sum()); p, lo, hi = wilson(k, len(g))
        rows.append((al, len(g), 100 * p, 100 * lo, 100 * hi))
    if not rows:
        ax.set_axis_off(); return
    t = (pd.DataFrame(rows, columns=["a", "n", "p", "lo", "hi"])
         .sort_values("p", ascending=False).head(top).iloc[::-1])
    y = np.arange(len(t))
    for yi, (_, r) in zip(y, t.iterrows()):
        col = fam_colour(r["a"])
        ax.barh(yi, r["p"], color=col, alpha=.85, height=.66,
                edgecolor="white", linewidth=.8)
        ax.errorbar(r["p"], yi, xerr=[[r["p"] - r["lo"]], [r["hi"] - r["p"]]],
                    color="#333", lw=1, capsize=2, capthick=1)
        ax.text(r["hi"] + 1.4, yi, f"n={r['n']:,}", va="center",
                fontsize=7.5, color=C["grey"])
    ax.set_yticks(y); ax.set_yticklabels(t["a"])
    ax.set_xlabel("phenotypically susceptible (%)")
    ax.set_title("Allele-level resolution", pad=8)
    ax.set_xlim(0, max(t["hi"]) * 1.3 + 5)
    ax.grid(axis="x", alpha=.25, lw=.6); ax.set_axisbelow(True)


def panel_drug_compare(ax, min_n=10, tag=""):
    """Meropenem vs imipenem discordance per family — the kinetic signature."""
    got = {}
    for drug in ("meropenem", "imipenem"):
        f = TAB / f"recordlevel_discordance_{drug}{tag}.tsv"
        if not f.exists():
            f = TAB / f"discordance_{drug}{tag}.tsv"
        if not f.exists():
            continue
        dd = load(drug, tag); p = dd[dd["carb_pos"]]
        r = {}
        for fam, g in p.groupby("carb_family"):
            if len(g) < min_n or ";" in str(fam) or "+" in str(fam):
                continue
            r[fam] = (100 * wilson(int(g["susceptible"].sum()), len(g))[0], len(g))
        got[drug] = r
    if len(got) < 2:
        ax.text(.5, .5, "imipenem table not found", ha="center", va="center",
                transform=ax.transAxes, color=C["grey"]); ax.set_axis_off(); return
    fams = [f for f in got["meropenem"] if f in got["imipenem"]]
    if not fams:
        ax.set_axis_off(); return
    x = np.arange(len(fams)); w = 0.36
    m = [got["meropenem"][f][0] for f in fams]
    i = [got["imipenem"][f][0] for f in fams]
    ax.bar(x - w / 2, m, w, label="meropenem", color="#2C7FB8",
           edgecolor="white", linewidth=.8)
    ax.bar(x + w / 2, i, w, label="imipenem", color="#F0A860",
           edgecolor="white", linewidth=.8)
    for xi, f in zip(x, fams):
        ax.text(xi - w / 2, m[fams.index(f)] + 1.2,
                f"{got['meropenem'][f][1]:,}", ha="center", fontsize=6.8,
                color=C["grey"])
        ax.text(xi + w / 2, i[fams.index(f)] + 1.2,
                f"{got['imipenem'][f][1]:,}", ha="center", fontsize=6.8,
                color=C["grey"])
    ax.set_xticks(x); ax.set_xticklabels(fams, rotation=18, ha="right")
    ax.set_ylabel("phenotypically susceptible (%)")
    ax.set_title("Substrate specificity signature", pad=8)
    ax.legend(frameon=False, loc="upper right", handlelength=1.1)
    ax.grid(axis="y", alpha=.25, lw=.6); ax.set_axisbelow(True)


def panel_mic_dist(ax, d, drug):
    pos = d[d["carb_pos"]]
    groups, labels, cols = [], [], []
    for fam in ["OXA-48-like", "KPC", "NDM"]:
        for pr in ["none", "ompK36_loop3"]:
            g = pos[(pos["carb_family"] == fam) & (pos["porin"] == pr)]
            if len(g) < 10:
                continue
            groups.append(np.log2(g["mic"].dropna()))
            labels.append(f"{fam}\n({PORIN_LBL[pr]})")
            cols.append(fam_colour(fam))
    if not groups:
        ax.set_axis_off(); return
    bp = ax.boxplot(groups, patch_artist=True, showfliers=False, widths=.6)
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(.55); patch.set_edgecolor("#333")
    for el in ("medians", "whiskers", "caps"):
        for ln in bp[el]:
            ln.set_color("#333"); ln.set_linewidth(1.1)
    for i, g in enumerate(groups, start=1):
        ax.scatter(np.random.normal(i, .055, len(g)), g, s=4,
                   color="#222", alpha=.16, zorder=1)
    ax.axhline(0, ls="--", lw=1, color=C["accent"])
    ax.text(len(groups) + .45, 0, "S breakpoint", fontsize=7,
            color=C["accent"], va="center")
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=6.8, rotation=30, ha="right",
                       rotation_mode="anchor")
    ax.set_ylabel(f"log2 {drug} MIC (mg/L)")
    ax.set_title("MIC distribution by genotype", pad=8)
    ax.grid(axis="y", alpha=.25, lw=.6); ax.set_axisbelow(True)


def panel_kinetic_contrast(ax, tag="", min_n=10):
    """OXA-48-like vs KPC discordance across meropenem/imipenem, annotated with
    the measured turnover numbers. This is the strongest mechanistic evidence:
    the drug-specific reversal follows the enzyme's substrate preference."""
    data = {}
    for drug in ("meropenem", "imipenem"):
        f = TAB / f"recordlevel_discordance_{drug}{tag}.tsv"
        if not f.exists():
            f = TAB / f"discordance_{drug}{tag}.tsv"
        if not f.exists():
            continue
        dd = load(drug, tag); pos = dd[dd["carb_pos"]]
        for fam in ("OXA-48-like", "KPC", "NDM"):
            g = pos[pos["carb_family"] == fam]
            if len(g) < min_n:
                continue
            p, lo, hi = wilson(int(g["susceptible"].sum()), len(g))
            data[(fam, drug)] = (100 * p, 100 * lo, 100 * hi, len(g))
    fams = [f for f in ("OXA-48-like", "KPC", "NDM")
            if (f, "meropenem") in data and (f, "imipenem") in data]
    if not fams:
        ax.text(.5, .5, "insufficient data", ha="center", va="center",
                transform=ax.transAxes, color=C["grey"]); ax.set_axis_off(); return
    x = np.arange(len(fams)); w = .34
    for off, drug, col in ((-w / 2, "meropenem", "#2C7FB8"),
                           (w / 2, "imipenem", "#F0A860")):
        vals = [data[(f, drug)][0] for f in fams]
        los = [data[(f, drug)][0] - data[(f, drug)][1] for f in fams]
        his = [data[(f, drug)][2] - data[(f, drug)][0] for f in fams]
        ax.bar(x + off, vals, w, label=drug, color=col,
               edgecolor="white", linewidth=.8)
        ax.errorbar(x + off, vals, yerr=[los, his], fmt="none",
                    ecolor="#333", lw=1.1, capsize=2.5)
        for xi, f in zip(x, fams):
            ax.text(xi + off, data[(f, drug)][2] + 1.4, f"{data[(f,drug)][3]:,}",
                    ha="center", fontsize=6.8, color=C["grey"])
    ax.set_xticks(x); ax.set_xticklabels(fams)
    ax.set_ylabel("phenotypically susceptible (%)")
    ax.set_title("Discordance follows substrate preference", pad=8)
    ax.legend(frameon=False, loc="upper right", handlelength=1.1)
    ax.grid(axis="y", alpha=.25, lw=.6); ax.set_axisbelow(True)
    # annotate with measured kinetics
    note = ("OXA-48 hydrolyses imipenem ~28x faster\nthan meropenem "
            "(kcat 2.8 vs ~0.1 s$^{-1}$)")
    ax.text(.02, .97, note, transform=ax.transAxes, va="top", ha="left",
            fontsize=7.2, color=C["grey"], style="italic")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="meropenem")
    ap.add_argument("--tag", default="",
                    help="table suffix, e.g. _clsi_bmd for the stratified analysis")
    a = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    d = load(a.drug, a.tag)

    # Vertical layout: 2 columns x 3 rows -> A B / C D / E F
    fig = plt.figure(figsize=(9.6, 13.2))
    gs = GridSpec(3, 2, figure=fig, hspace=.58, wspace=.34)

    axA = fig.add_subplot(gs[0, 0]); panel_confusion(axA, d, a.drug);   panel_label(axA, "A", dx=-0.34)
    axB = fig.add_subplot(gs[0, 1]); panel_family(axB, d);              panel_label(axB, "B", dx=-0.30)
    axC = fig.add_subplot(gs[1, 0]); panel_mechanism(axC, d);           panel_label(axC, "C", dx=-0.22)
    axD = fig.add_subplot(gs[1, 1]); panel_allele(axD, d);              panel_label(axD, "D", dx=-0.32)
    axE = fig.add_subplot(gs[2, 0]); panel_kinetic_contrast(axE, a.tag); panel_label(axE, "E", dx=-0.22)
    axF = fig.add_subplot(gs[2, 1]); panel_mic_dist(axF, d, a.drug);    panel_label(axF, "F", dx=-0.22)

    # No title or caption is drawn inside the figure.

    # Name the output by its role in the manuscript, not by the template that
    # produced it. All of these share the same six-panel layout, which is why
    # they were previously all prefixed "figure2_" and impossible to tell apart.
    ROLE = {
        ("meropenem", "_clsi_bmd"):        "Figure2_main_meropenem_CLSI",
        ("imipenem",  "_clsi_bmd"):        "FigureS1_imipenem_CLSI",
        ("meropenem", "_unstratified"):    "FigureS2_meropenem_unstratified",
        ("meropenem", "_clsi_bmd_dedup"):  "FigureS3_meropenem_deduplicated",
        ("meropenem", "_clsi_bmd_eucast"): "FigureS4_meropenem_EUCAST",
    }
    stem = ROLE.get((a.drug, a.tag),
                    f"figure_discordance_{a.drug}{a.tag}")
    out = FIG / stem
    figstyle.save(fig, out)
    print(f"[out] {out.name}.png / .pdf / .svg")
    print(f"      binary S/R isolates plotted: {len(d):,} "
          f"(of {d.attrs.get('n_all', len(d)):,} linked)")


if __name__ == "__main__":
    main()
