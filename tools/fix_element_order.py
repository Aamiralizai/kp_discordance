#!/usr/bin/env python3
"""fix_element_order.py — put run and paragraph properties back in schema order.

WordprocessingML fixes the order of the children of <w:rPr> and <w:pPr>. A
revision marker on a paragraph mark, <w:ins> inside <w:rPr>, has to come first,
and a formatting-revision record, <w:rPrChange> or <w:pPrChange>, has to come
last. Writing them anywhere else makes Word refuse the file with "unreadable
content", which is how this was found.

The fix is a stable sort of the children by their position in the schema
sequence. Nothing is added or removed, so the text and the revisions are
untouched.

USAGE
    python3 fix_element_order.py word/document.xml [--apply]
"""
import sys

from lxml import etree

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"


def t(e):
    return e.tag.rsplit("}", 1)[-1] if isinstance(e.tag, str) else None


RPR = ["ins", "del", "moveFrom", "moveTo", "rStyle", "rFonts", "b", "bCs", "i",
       "iCs", "caps", "smallCaps", "strike", "dstrike", "outline", "shadow",
       "emboss", "imprint", "noProof", "snapToGrid", "vanish", "webHidden",
       "color", "spacing", "w", "kern", "position", "sz", "szCs", "highlight",
       "u", "effect", "bdr", "shd", "fitText", "vertAlign", "rtl", "cs", "em",
       "lang", "eastAsianLayout", "specVanish", "oMath", "rPrChange"]

PPR = ["pStyle", "keepNext", "keepLines", "pageBreakBefore", "framePr",
       "widowControl", "numPr", "suppressLineNumbers", "pBdr", "shd", "tabs",
       "suppressAutoHyphens", "kinsoku", "wordWrap", "overflowPunct",
       "topLinePunct", "autoSpaceDE", "autoSpaceDN", "bidi", "adjustRightInd",
       "snapToGrid", "spacing", "ind", "contextualSpacing", "mirrorIndents",
       "suppressOverlap", "jc", "textDirection", "textAlignment",
       "textboxTightWrap", "outlineLvl", "divId", "cnfStyle", "rPr",
       "sectPr", "pPrChange"]


def order(parent, seq):
    idx = {n: i for i, n in enumerate(seq)}
    kids = list(parent)
    if any(t(k) not in idx for k in kids):
        return False
    want = sorted(kids, key=lambda k: idx[t(k)])
    if want == kids:
        return False
    for k in want:
        parent.append(k)
    return True


def main(path, apply=False):
    tree = etree.parse(path)
    root = tree.getroot()
    fixed = 0
    for name, seq in (("rPr", RPR), ("pPr", PPR)):
        for el in root.iter(f"{{{W}}}{name}"):
            if order(el, seq):
                fixed += 1
    print(f"  {fixed} properties element(s) reordered")
    if fixed and apply:
        tree.write(path, xml_declaration=True, encoding="UTF-8", standalone=True)
        print("  APPLIED")
    elif fixed:
        print("  dry run, pass --apply to write")
    return fixed


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], "--apply" in sys.argv)
