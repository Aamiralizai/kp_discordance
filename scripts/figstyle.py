#!/usr/bin/env python
"""figstyle.py — shared figure styling and export.

Two things every figure script needs:

  * Times New Roman with metric-compatible fallbacks. Times is present on
    Windows and macOS but usually absent on Linux/WSL, where Liberation Serif
    is the metric-compatible substitute — the rendered PNG/PDF will look
    correct either way.

  * EDITABLE TEXT in the exported SVG. matplotlib converts glyphs to paths by
    default, which makes subheadings and axis labels impossible to edit. We set
    svg.fonttype="none" to emit real <text> elements, then rewrite the
    font-family declaration so the SVG explicitly asks for Times New Roman even
    when it was rendered with a substitute.
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt

SERIF_STACK = ["Times New Roman", "Liberation Serif", "Nimbus Roman",
               "FreeSerif", "DejaVu Serif"]
SVG_FONT_FAMILY = "'Times New Roman', 'Liberation Serif', Times, serif"


def apply_style(extra=None):
    """Journal-style serif defaults with editable vector text."""
    rc = {
        "font.family": "serif",
        "font.serif": SERIF_STACK,
        "mathtext.fontset": "stix",
        "svg.fonttype": "none",     # keep text as <text>, not paths
        "pdf.fonttype": 42,         # TrueType: editable in vector editors
        "ps.fonttype": 42,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    }
    if extra:
        rc.update(extra)
    plt.rcParams.update(rc)


def _fix_svg_fonts(path: Path):
    """Name the intended font in the SVG, whatever was used to render it."""
    s = path.read_text(encoding="utf-8")
    # matplotlib may emit an empty or substitute family; normalise both
    s = re.sub(r"font-family:\s*[^;\"']*;", f"font-family:{SVG_FONT_FAMILY};", s)
    s = re.sub(r'font-family="[^"]*"', f'font-family="{SVG_FONT_FAMILY}"', s)
    if "font-family" not in s:
        s = s.replace("<svg ", f'<svg style="font-family:{SVG_FONT_FAMILY}" ', 1)
    path.write_text(s, encoding="utf-8")


def save(fig, out_stem: Path, formats=(".png", ".pdf", ".svg")):
    """Save in each format; post-process the SVG so its text stays editable
    and names Times New Roman."""
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for ext in formats:
        p = out_stem.with_suffix(ext)
        fig.savefig(p)
        if ext == ".svg":
            _fix_svg_fonts(p)
        written.append(p)
    return written
