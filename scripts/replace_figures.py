#!/usr/bin/env python3
"""replace_figures.py — swap the figures embedded in a Word document for the
regenerated ones, resizing the display box to each new aspect ratio.

The images live in word/media and are referenced by relationship id. Replacing
the bytes alone leaves the stored width and height untouched, so a figure whose
proportions changed comes out stretched. Both the drawing's wp:extent and the
picture's a:ext are rewritten here from the new image, keeping the column width
the document already uses.
"""
import shutil
import sys
from pathlib import Path

from lxml import etree
from PIL import Image

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"


def w(t):
    return f"{{{W}}}{t}"


def main(docdir, width_emu, mapping):
    docdir = Path(docdir)
    rels = {r.get("Id"): r.get("Target")
            for r in etree.parse(docdir / "word/_rels/document.xml.rels").getroot()}
    tree = etree.parse(str(docdir / "word/document.xml"))
    root = tree.getroot()

    # drawings in document order, each paired with the caption that follows it
    drawings = []
    for p in root.find(w("body")):
        if p.tag != w("p"):
            continue
        for blip in p.findall(f".//{{{A}}}blip"):
            drawings.append((p, blip))

    if len(drawings) != len(mapping):
        sys.exit(f"{docdir}: {len(drawings)} embedded images but "
                 f"{len(mapping)} replacements given")

    for (para, blip), newpng in zip(drawings, mapping):
        rid = blip.get(f"{{{R}}}embed")
        target = rels[rid]
        dest = docdir / "word" / target
        shutil.copyfile(newpng, dest)
        im = Image.open(newpng)
        cx = width_emu
        cy = int(round(width_emu * im.size[1] / im.size[0]))
        for ext in para.iter(f"{{{WP}}}extent"):
            ext.set("cx", str(cx)); ext.set("cy", str(cy))
        for ext in para.iter(f"{{{A}}}ext"):
            if ext.get("cx") is not None:
                ext.set("cx", str(cx)); ext.set("cy", str(cy))
        print(f"  {target} <- {Path(newpng).name}  "
              f"{im.size[0]}x{im.size[1]}  box {cx}x{cy} EMU")

    tree.write(str(docdir / "word/document.xml"), xml_declaration=True,
               encoding="UTF-8", standalone=True)


if __name__ == "__main__":
    docdir, width = sys.argv[1], int(sys.argv[2])
    main(docdir, width, sys.argv[3:])
