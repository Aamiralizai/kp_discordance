#!/usr/bin/env python
"""15_figure3_robustness.py — Figure 3: robustness of the findings.

Four panels showing how the reported quantities vary across every analysis
variant that was run:

  A  headline estimate (carbapenemase-positive classified susceptible)
  B  the same estimate by carbapenemase family
  C  predictive performance of carbapenemase detection
  D  the OXA-48-like porin gradient — the central mechanistic result

All values are read from results/manuscript_results.json and the record-level
mechanism tables, so the figure cannot drift from the reported numbers.
No caption is drawn inside the figure.

Usage:
    python scripts/15_figure3_robustness.py
"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parent))
import figstyle
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

RES = Path("results/manuscript_results.json")
TAB = Path("results/tables")
FIG = Path("results/figures")

figstyle.apply_style()

# analysis key -> (short label, tag, colour)
VARIANTS = [
    ("primary (CLSI, broth dilution)", "CLSI\nbroth dilution", "_clsi_bmd", "#1B5E8C"),
    ("meropenem (EUCAST breakpoints)", "EUCAST\nv16.1", "_clsi_bmd_eucast", "#E8743B"),
    ("meropenem deduplicated", "deduplicated\nby genotype", "_clsi_bmd_dedup", "#2E7D5B"),
    ("meropenem unstratified", "all methods\nand standards", "_unstratified", "#6A3D9A"),
]
FAM_COL = {"OXA-48-like": "#E8743B", "KPC": "#1F78B4",
           "NDM": "#33A02C", "VIM": "#6A3D9A"}
PORIN_ORDER = ["none", "ompK35_loss", "ompK36_loop3"]
PORIN_LBL = {"none": "intact", "ompK35_loss": "OmpK35\ntruncation",
             "ompK36_loop3": "OmpK36\nloop-3"}


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p, z = k / n, 1.959964
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (100 * p, 100 * max(0, c - h), 100 * min(1, c + h))


def label(ax, letter, dx=-0.14, dy=1.10):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")


# ------------------------------------------------------------------ panels
def panel_headline(ax, rows):
    y = np.arange(len(rows))[::-1]
    prim = rows[0]
    ax.axvspan(prim["lo"], prim["hi"], color=prim["colour"], alpha=.07, zorder=0)
    ax.axvline(prim["p"], color=prim["colour"], lw=1, ls="--", alpha=.5, zorder=1)
    for yi, r in zip(y, rows):
        ax.plot([r["lo"], r["hi"]], [yi, yi], color=r["colour"], lw=2.6,
                solid_capstyle="round", zorder=2)
        ax.scatter(r["p"], yi, s=70, color=r["colour"], zorder=3,
                   edgecolor="white", linewidth=1.1)
        ax.text(r["hi"] + 1.0, yi, f"{r['p']:.1f}%", va="center",
                fontsize=8.4, color="#2A3B47")
        ax.text(-1.5, yi, f"{r['k']}/{r['n']:,}", va="center", ha="right",
                fontsize=7.6, color="#6A7A87")
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"].replace("\n", " ") for r in rows])
    ax.set_xlabel("carbapenemase-positive classified susceptible (%)")
    ax.set_title("Headline estimate", pad=8)
    ax.set_xlim(-7, max(r["hi"] for r in rows) * 1.28)
    ax.grid(axis="x", alpha=.22, lw=.6); ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)


def panel_by_family(ax, data):
    fams = ["OXA-48-like", "KPC", "NDM"]
    x = np.arange(len(fams)); w = 0.2
    for i, (key, lab, tag, col) in enumerate(VARIANTS):
        a = data.get(key)
        if not a:
            continue
        vals, los, his, ns = [], [], [], []
        for f in fams:
            e = a["by_family"].get(f)
            if e:
                p, lo, hi = e["pct_ci"]
                vals.append(p)
                los.append(max(0.0, p - lo)); his.append(max(0.0, hi - p))
                ns.append(e["n"])
            else:
                vals.append(np.nan); los.append(0.0); his.append(0.0); ns.append(0)
        off = (i - (len(VARIANTS) - 1) / 2) * w
        ax.bar(x + off, vals, w * 0.92, color=col, label=lab.replace("\n", " "),
               edgecolor="white", linewidth=.6)
        ax.errorbar(x + off, vals, yerr=[los, his], fmt="none",
                    ecolor="#3A4A55", lw=.9, capsize=2)
    ax.set_xticks(x); ax.set_xticklabels(fams)
    ax.set_ylabel("classified susceptible (%)")
    ax.set_title("By carbapenemase family", pad=8)
    ax.legend(frameon=False, ncol=2, loc="upper right", handlelength=1.0,
              columnspacing=1.0, borderpad=.2)
    ax.grid(axis="y", alpha=.22, lw=.6); ax.set_axisbelow(True)
    ax.margins(y=.22)


def panel_metrics(ax, data):
    metrics = [("PPV", "PPV"), ("sensitivity", "sensitivity"),
               ("specificity", "specificity")]
    x = np.arange(len(metrics)); w = 0.2
    for i, (key, lab, tag, col) in enumerate(VARIANTS):
        a = data.get(key)
        if not a:
            continue
        vals = [a.get(k) for k, _ in metrics]
        off = (i - (len(VARIANTS) - 1) / 2) * w
        ax.bar(x + off, vals, w * 0.92, color=col, edgecolor="white", linewidth=.6)
    ax.set_xticks(x); ax.set_xticklabels([l for _, l in metrics])
    ax.set_ylabel("value")
    ax.set_ylim(0, 1.05)
    ax.set_title("Performance of carbapenemase detection", pad=8)
    ax.grid(axis="y", alpha=.22, lw=.6); ax.set_axisbelow(True)


def panel_porin(ax):
    xs, seen = np.arange(len(PORIN_ORDER)), False
    w = 0.2
    for i, (key, lab, tag, col) in enumerate(VARIANTS):
        f = TAB / f"recordlevel_mechanism_meropenem{tag}.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t")
        d = d[d["family"] == "OXA-48-like"]
        vals, los, his = [], [], []
        for pr in PORIN_ORDER:
            r = d[d["porin"] == pr]
            if len(r):
                p, lo, hi = r.iloc[0]["pct"], r.iloc[0]["ci_lo"], r.iloc[0]["ci_hi"]
                vals.append(p)
                los.append(max(0.0, p - lo)); his.append(max(0.0, hi - p))
            else:
                vals.append(np.nan); los.append(0); his.append(0)
        off = (i - (len(VARIANTS) - 1) / 2) * w
        ax.bar(xs + off, vals, w * 0.92, color=col, edgecolor="white", linewidth=.6)
        ax.errorbar(xs + off, vals, yerr=[los, his], fmt="none",
                    ecolor="#3A4A55", lw=.9, capsize=2)
        seen = True
    if not seen:
        ax.set_axis_off(); return
    ax.set_xticks(xs)
    ax.set_xticklabels([PORIN_LBL[p] for p in PORIN_ORDER])
    ax.set_ylabel("classified susceptible (%)")
    ax.set_title("OXA-48-like: porin gradient", pad=8)
    ax.grid(axis="y", alpha=.22, lw=.6); ax.set_axisbelow(True)
    ax.margins(y=.20)


def main():
    if not RES.exists():
        sys.exit(f"\nMissing {RES}. Run scripts/14_results_manifest.py first.\n")
    data = {a["analysis"]: a for a in json.load(RES.open())["analyses"]}

    rows = []
    for key, lab, tag, col in VARIANTS:
        a = data.get(key)
        if not a:
            continue
        p, lo, hi = a["carbapenemase_positive_susceptible_pct_ci"]
        rows.append(dict(label=lab, colour=col, p=p, lo=lo, hi=hi,
                         n=a["carbapenemase_positive_n"],
                         k=a["carbapenemase_positive_susceptible_n"]))
    if not rows:
        sys.exit("No analyses found in the manifest.")

    FIG.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(11.6, 7.4))
    gs = GridSpec(2, 2, figure=fig, hspace=.52, wspace=.34)

    axA = fig.add_subplot(gs[0, 0]); panel_headline(axA, rows);  label(axA, "A", dx=-0.34)
    axB = fig.add_subplot(gs[0, 1]); panel_by_family(axB, data); label(axB, "B")
    axC = fig.add_subplot(gs[1, 0]); panel_metrics(axC, data);   label(axC, "C", dx=-0.18)
    axD = fig.add_subplot(gs[1, 1]); panel_porin(axD);           label(axD, "D")

    out = FIG / "Figure3_robustness_meropenem"
    figstyle.save(fig, out)
    print(f"[out] {out}.png / .pdf / .svg")
    for r in rows:
        print(f"  {r['label'].replace(chr(10), ' '):<28}{r['p']:>6.1f}% "
              f"[{r['lo']:.1f}-{r['hi']:.1f}]")


if __name__ == "__main__":
    main()
