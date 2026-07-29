#!/usr/bin/env python
"""16_figure4_false_resistant.py — Figure 4: the other direction of discordance.

The paper reports discordance in BOTH directions, but Figures 2 and 3 concern
carbapenemase-POSITIVE isolates only. This figure covers the reverse: isolates
that are phenotypically resistant with no acquired carbapenemase.

  A  resistance rate by ESBL determinant x OmpK36 disruption
  B  association of each beta-lactam determinant with resistance (OR, 95% CI)
  C  how much of the resistant, carbapenemase-negative group is explained
  D  the two directions of discordance side by side

Everything is computed from the record-level table, so the figure is
self-contained and cannot drift from the reported numbers. No caption is drawn
inside the figure.

Usage:
    python scripts/16_figure4_false_resistant.py
"""
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
from scipy import stats

TAB = Path("results/tables")
FIG = Path("results/figures")

figstyle.apply_style()

C_POS = "#B04A3F"     # resistant / positive association
C_NEG = "#4A7FA5"     # susceptible / inverse association
C_NEU = "#8C99A3"
C_EXP = "#2E7D5B"
C_UNEXP = "#9B3B2E"


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p, z = k / n, 1.959964
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (100 * p, 100 * max(0, c - h), 100 * min(1, c + h))


def or_ci(a11, a12, a21, a22):
    """Odds ratio with Woolf 95% CI; Haldane correction if a cell is zero."""
    if 0 in (a11, a12, a21, a22):
        a11, a12, a21, a22 = a11 + .5, a12 + .5, a21 + .5, a22 + .5
    orr = (a22 / a21) / (a12 / a11)
    se = math.sqrt(1 / a11 + 1 / a12 + 1 / a21 + 1 / a22)
    return orr, orr * math.exp(-1.96 * se), orr * math.exp(1.96 * se)


def label(ax, letter, dx=-0.15, dy=1.10):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")


def load():
    f = TAB / "recordlevel_discordance_meropenem_clsi_bmd.tsv"
    if not f.exists():
        sys.exit(f"\nMissing {f}. Run 08b_discordance_recordlevel.py first.\n")
    d = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
    for c in ("carb_pos", "porin_disrupted", "esbl"):
        d[c] = d[c].astype(str).str.strip().str.lower().isin({"true", "1"})
    d = d[d["call"].isin(["S", "R"])].copy()
    d["res"] = d["call"] == "R"
    return d


# ------------------------------------------------------------------ panels
def panel_cross(ax, neg):
    cells = [(False, False, "neither"), (True, False, "ESBL only"),
             (False, True, "OmpK36\ndisruption only"), (True, True, "both")]
    palette = [C_NEU, C_NEG, C_POS, "#D4826A"]
    xs, vals, los, his, ns, cols = [], [], [], [], [], []
    for (e, p, lab), col in zip(cells, palette):
        g = neg[(neg["esbl"] == e) & (neg["porin_disrupted"] == p)]
        if len(g) == 0:          # empty cell: nothing to plot
            continue
        v, lo, hi = wilson(int(g["res"].sum()), len(g))
        xs.append(lab); vals.append(v); ns.append(len(g)); cols.append(col)
        los.append(max(0.0, v - lo)); his.append(max(0.0, hi - v))
    if not vals:
        ax.set_axis_off(); return
    x = np.arange(len(vals))
    ax.bar(x, vals, .62, color=cols, edgecolor="white", linewidth=.8)
    ax.errorbar(x, vals, yerr=[los, his], fmt="none", ecolor="#3A4A55",
                lw=1.0, capsize=3)
    for xi, v, hi, n in zip(x, vals, his, ns):
        ax.text(xi, v + hi + 1.8, f"n={n:,}", ha="center", fontsize=7.6,
                color="#6A7A87")
    ax.set_xticks(x); ax.set_xticklabels(xs, fontsize=8)
    ax.set_ylabel("phenotypically resistant (%)")
    ax.set_title("Resistance without a carbapenemase", pad=8)
    ax.grid(axis="y", alpha=.22, lw=.6); ax.set_axisbelow(True)
    ax.margins(y=.22)


def panel_determinants(ax, neg):
    dets = [("porin_disrupted", "OmpK36 disruption"),
            ("esbl", "acquired ESBL")]
    rows = []
    for col, name in dets:
        a11 = int((~neg[col] & ~neg["res"]).sum())
        a12 = int((~neg[col] & neg["res"]).sum())
        a21 = int((neg[col] & ~neg["res"]).sum())
        a22 = int((neg[col] & neg["res"]).sum())
        orr, lo, hi = or_ci(a11, a12, a21, a22)
        _, pv = stats.fisher_exact([[a11, a12], [a21, a22]])
        rows.append(dict(name=name, orr=orr, lo=lo, hi=hi, p=pv,
                         n=int(neg[col].sum())))
    y = np.arange(len(rows))[::-1]
    ax.axvline(1, color="#6A7A87", lw=1, ls="--", alpha=.7)
    for yi, r in zip(y, rows):
        col = C_POS if r["orr"] > 1 else C_NEG
        ax.plot([r["lo"], r["hi"]], [yi, yi], color=col, lw=2.6,
                solid_capstyle="round")
        ax.scatter(r["orr"], yi, s=72, color=col, zorder=3,
                   edgecolor="white", linewidth=1.1)
        ax.text(r["hi"] * 1.15, yi, f"{r['orr']:.2f}", va="center",
                fontsize=8.4, color="#2A3B47")
    ax.set_yticks(y); ax.set_yticklabels([r["name"] for r in rows])
    ax.set_xscale("log")
    ax.set_xlabel("odds ratio for resistance (log scale)")
    ax.set_title("Determinant associations", pad=8)
    ax.grid(axis="x", alpha=.22, lw=.6); ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.text(0.02, -0.42,
            "inverse associations are attributed to collection composition,\n"
            "not to a protective effect (see text)",
            transform=ax.transAxes, fontsize=7.4, color="#6A7A87", style="italic")


def panel_explained(ax, neg):
    fr = neg[neg["res"]]
    n = len(fr)
    porin = int(fr["porin_disrupted"].sum())
    unexp = n - porin
    vals = [porin, unexp]
    labs = [f"OmpK36 disruption\n{porin} ({100*porin/n:.1f}%)",
            f"not explained\n{unexp} ({100*unexp/n:.1f}%)"]
    ax.barh([1, 0], vals, .5, color=[C_EXP, C_UNEXP],
            edgecolor="white", linewidth=.8)
    for yi, v, lab in zip([1, 0], vals, labs):
        ax.text(v + n * 0.02, yi, lab, va="center", fontsize=8.4,
                color="#2A3B47")
    ax.set_yticks([]); ax.set_xlim(0, n * 1.5)
    ax.set_xlabel(f"resistant isolates without a carbapenemase (n = {n})")
    ax.set_title("How much is explained", pad=8)
    ax.grid(axis="x", alpha=.22, lw=.6); ax.set_axisbelow(True)
    for s in ("left",):
        ax.spines[s].set_visible(False)


def panel_directions(ax, d):
    pos, neg = d[d["carb_pos"]], d[~d["carb_pos"]]
    fs = int((pos["call"] == "S").sum()); fr = int(neg["res"].sum())
    a = wilson(fs, len(pos)); b = wilson(fr, len(neg))
    x = np.arange(2)
    vals = [a[0], b[0]]
    los = [max(0.0, a[0] - a[1]), max(0.0, b[0] - b[1])]
    his = [max(0.0, a[2] - a[0]), max(0.0, b[2] - b[0])]
    ax.bar(x, vals, .5, color=[C_NEG, C_POS], edgecolor="white", linewidth=.8)
    ax.errorbar(x, vals, yerr=[los, his], fmt="none", ecolor="#3A4A55",
                lw=1.0, capsize=3)
    for xi, v, hi, k, n in zip(x, vals, his, (fs, fr), (len(pos), len(neg))):
        ax.text(xi, v + hi + 0.5, f"{k}/{n:,}", ha="center", fontsize=7.8,
                color="#6A7A87")
    ax.set_xticks(x)
    ax.set_xticklabels(["carbapenemase +\nphenotypically S",
                        "carbapenemase −\nphenotypically R"], fontsize=8.2)
    ax.set_ylabel("% of that genotype group")
    ax.set_title("Both directions of discordance", pad=8)
    ax.grid(axis="y", alpha=.22, lw=.6); ax.set_axisbelow(True)
    ax.margins(y=.26)


def main():
    d = load()
    neg = d[~d["carb_pos"]].copy()
    FIG.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(11.6, 7.6))
    gs = GridSpec(2, 2, figure=fig, hspace=.62, wspace=.36)
    axA = fig.add_subplot(gs[0, 0]); panel_cross(axA, neg);        label(axA, "A")
    axB = fig.add_subplot(gs[0, 1]); panel_determinants(axB, neg); label(axB, "B", dx=-0.26)
    axC = fig.add_subplot(gs[1, 0]); panel_explained(axC, neg);    label(axC, "C")
    axD = fig.add_subplot(gs[1, 1]); panel_directions(axD, d);     label(axD, "D")

    out = FIG / "Figure4_carbapenemase_negative_meropenem"
    figstyle.save(fig, out)
    print(f"[out] {out}.png / .pdf / .svg")

    fr = neg[neg["res"]]
    print(f"  carbapenemase-negative: {len(neg):,}; resistant {len(fr)} "
          f"({100*len(fr)/len(neg):.1f}%)")
    print(f"  of those, OmpK36-disrupted: {int(fr['porin_disrupted'].sum())} "
          f"({100*fr['porin_disrupted'].mean():.1f}%)")


if __name__ == "__main__":
    main()
