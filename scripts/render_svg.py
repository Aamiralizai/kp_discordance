#!/usr/bin/env python3
"""render_svg.py — turn an SVG into PNG, PDF and TIFF.

The figures of record in this project are SVG. Journals ask for raster and
print copies too, and different machines have different renderers installed,
so this module tries them in order of fidelity and reports which one it used
rather than failing when the first is absent:

    inkscape        best: honours the fonts and the markers exactly
    rsvg-convert    good
    headless Chromium (Playwright's copy, or any chrome on PATH)
    cairosvg        last resort: no font substitution table

Chromium needs one quirk: asked for a window the exact size of the drawing it
sometimes returns a shorter one and silently crops the bottom of the figure, so
it is asked for a larger window and the screenshot is cropped back.

USAGE
    python3 scripts/render_svg.py FIGURE.svg [--scale 4] [--formats png,pdf,tif]
    from render_svg import render; render("fig.svg", scale=4)
"""
import os
import re
import shutil
import subprocess
import sys

CHROME_CANDIDATES = [
    "/opt/pw-browsers/chromium-1194/chrome-linux/chrome",
    "/opt/pw-browsers/chromium/chrome-linux/chrome",
]


def _chrome():
    for p in CHROME_CANDIDATES:
        if os.path.exists(p):
            return p
    for name in ("chromium", "chromium-browser", "google-chrome", "chrome"):
        p = shutil.which(name)
        if p:
            return p
    for root in ("/opt/pw-browsers",):
        if os.path.isdir(root):
            for dirpath, _dirs, files in os.walk(root):
                for f in ("chrome", "headless_shell"):
                    if f in files:
                        return os.path.join(dirpath, f)
    return None


def _size(svg_text):
    m = re.search(r'viewBox="([-\d.eE\s]+)"', svg_text)
    if m:
        vb = [float(x) for x in m.group(1).split()]
        return vb[2], vb[3]
    w = re.search(r'<svg\b[^>]*\bwidth="([\d.]+)', svg_text)
    h = re.search(r'<svg\b[^>]*\bheight="([\d.]+)', svg_text)
    if w and h:
        return float(w.group(1)), float(h.group(1))
    raise ValueError("cannot determine the SVG's size")


def _chrome_png(svg_path, out_png, scale):
    chrome = _chrome()
    if not chrome:
        return False
    svg = open(svg_path, encoding="utf-8").read()
    w, h = _size(svg)
    svg = re.sub(r'(<svg\b[^>]*?)\s(?:width|height)="[^"]*"', r"\1", svg, count=2)
    svg = svg.replace("<svg", f'<svg width="{w}" height="{h}"', 1)
    html = ("<!doctype html><meta charset='utf-8'><style>"
            "html,body{margin:0;padding:0;background:#fff}"
            f"svg{{display:block;width:{w}px;height:{h}px}}</style>" + svg)
    tmp = os.path.join(os.path.dirname(os.path.abspath(out_png)),
                       f"._render_{os.path.basename(out_png)}.html")
    open(tmp, "w", encoding="utf-8").write(html)
    try:
        subprocess.run([chrome, "--headless", "--disable-gpu", "--no-sandbox",
                        "--hide-scrollbars", "--default-background-color=FFFFFFFF",
                        f"--force-device-scale-factor={scale}",
                        f"--window-size={int(w) + 200},{int(h) + 300}",
                        f"--screenshot={out_png}", f"file://{tmp}"],
                       capture_output=True, timeout=180)
    except Exception:
        return False
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    if not os.path.exists(out_png):
        return False
    try:
        from PIL import Image
        im = Image.open(out_png).convert("RGB")
        im.crop((0, 0, int(round(w * scale)), int(round(h * scale)))).save(out_png)
    except Exception:
        pass
    return True


def _chrome_pdf(svg_path, out_pdf):
    """A vector PDF: Chromium prints the page rather than photographing it."""
    chrome = _chrome()
    if not chrome:
        return False
    svg = open(svg_path, encoding="utf-8").read()
    w, h = _size(svg)
    svg = re.sub(r'(<svg\b[^>]*?)\s(?:width|height)="[^"]*"', r"\1", svg, count=2)
    svg = svg.replace("<svg", f'<svg width="{w}" height="{h}"', 1)
    # CSS px are 1/96 in, which is also how Chromium measures @page
    html = ("<!doctype html><meta charset='utf-8'><style>"
            f"@page{{size:{w / 96:.4f}in {h / 96:.4f}in;margin:0}}"
            "html,body{margin:0;padding:0;background:#fff}"
            f"svg{{display:block;width:{w}px;height:{h}px}}</style>" + svg)
    tmp = os.path.join(os.path.dirname(os.path.abspath(out_pdf)),
                       f"._print_{os.path.basename(out_pdf)}.html")
    open(tmp, "w", encoding="utf-8").write(html)
    try:
        subprocess.run([chrome, "--headless", "--disable-gpu", "--no-sandbox",
                        "--no-pdf-header-footer", f"--print-to-pdf={out_pdf}",
                        f"file://{tmp}"], capture_output=True, timeout=180)
    except Exception:
        return False
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return os.path.exists(out_pdf) and os.path.getsize(out_pdf) > 1000


def render(svg_path, scale=4, formats=("png", "pdf", "tif"), quiet=False):
    """Write the raster and print copies beside the SVG. Returns the renderer."""
    stem = os.path.splitext(svg_path)[0]
    out = {f: f"{stem}.{f}" for f in formats}
    used = None

    ink = shutil.which("inkscape")
    if ink:
        for f in ("png", "pdf"):
            if f in out:
                cmd = [ink, svg_path, f"--export-type={f}",
                       f"--export-filename={out[f]}"]
                if f == "png":
                    cmd.append(f"--export-dpi={int(96 * scale)}")
                subprocess.run(cmd, capture_output=True)
        used = "inkscape"
    elif shutil.which("rsvg-convert"):
        if "png" in out:
            subprocess.run(["rsvg-convert", "-z", str(scale), "-o", out["png"],
                            svg_path], capture_output=True)
        if "pdf" in out:
            subprocess.run(["rsvg-convert", "-f", "pdf", "-o", out["pdf"],
                            svg_path], capture_output=True)
        used = "rsvg-convert"
    elif "png" in out and _chrome_png(svg_path, out["png"], scale):
        used = "chromium"
        if "pdf" in out and not _chrome_pdf(svg_path, out["pdf"]):
            out.pop("pdf", None)
    else:
        try:
            import cairosvg
            if "png" in out:
                cairosvg.svg2png(url=svg_path, write_to=out["png"], scale=scale)
            if "pdf" in out:
                cairosvg.svg2pdf(url=svg_path, write_to=out["pdf"])
            used = "cairosvg"
        except Exception as exc:
            if not quiet:
                print(f"  no SVG renderer available ({exc}); SVG only")
            return None

    if "tif" in out and os.path.exists(out.get("png", "")):
        made = False
        if shutil.which("convert"):
            made = subprocess.run(["convert", out["png"], "-compress", "lzw",
                                   out["tif"]], capture_output=True).returncode == 0
        if not made:
            try:
                from PIL import Image
                Image.open(out["png"]).convert("RGB").save(
                    out["tif"], compression="tiff_lzw",
                    dpi=(96 * scale, 96 * scale))
            except Exception:
                out.pop("tif", None)

    if not quiet:
        for f in formats:
            p = out.get(f)
            if p and os.path.exists(p):
                print(f"    {p}  ({used})")
    return used


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    scale, formats = 4.0, ("png", "pdf", "tif")
    if "--scale" in args:
        i = args.index("--scale")
        scale = float(args[i + 1]); del args[i:i + 2]
    if "--formats" in args:
        i = args.index("--formats")
        formats = tuple(args[i + 1].split(",")); del args[i:i + 2]
    for p in args:
        render(p, scale=scale, formats=formats)
