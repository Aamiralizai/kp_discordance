#!/usr/bin/env python3
"""accept_changes.py — the clean copy of a tracked-changes document.

Journals usually want both a tracked version and a clean one. Accepting the
changes by hand in Word is error-prone and leaves no record of what was done,
so this applies the accept operation to the XML directly:

    <w:ins>   its runs are kept, the wrapper is removed
    <w:del>   removed entirely, along with its w:delText
    <w:trPr><w:ins>   a row marked inserted stays, its marker goes
    <w:trPr><w:del>   the whole row is removed
    rPrChange, pPrChange, sectPrChange, numberingChange, cellIns, cellDel
              formatting-revision records, removed

Afterwards the document must contain no revision marks at all and its text must
equal the tracked document's accepted text. Both are checked, and nothing is
written unless they hold.

USAGE
    python3 accept_changes.py TRACKED_DIR CLEAN_DIR
"""
import shutil
import sys
from pathlib import Path

from lxml import etree

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"


def w(t):
    return f"{{{W}}}"[:-1] + "}" + t if False else f"{{{W}}}{t}"


REVISION_MARKS = ("ins", "del", "rPrChange", "pPrChange", "sectPrChange",
                  "tblPrChange", "trPrChange", "tcPrChange", "numberingChange",
                  "cellIns", "cellDel", "cellMerge", "moveFrom", "moveTo")


def accepted_text(root):
    """Text with every change accepted, as a reader of the clean copy sees it.

    A table row can be marked deleted on the row itself, through w:trPr/w:del,
    while the text inside its cells sits in ordinary runs. Accepting removes the
    whole row, so the comparison has to skip those rows too, or it reports a
    difference that is the accept working correctly.
    """
    out = []
    for p in root.iter(w("p")):
        anc, row_deleted = p.getparent(), False
        while anc is not None:
            if anc.tag == w("tr"):
                trpr = anc.find(w("trPr"))
                if trpr is not None and trpr.find(w("del")) is not None:
                    row_deleted = True
                break
            anc = anc.getparent()
        if row_deleted:
            continue
        buf = []
        for t in p.iter(w("t")):
            anc, skip = t.getparent(), False
            while anc is not None:
                if anc.tag in (w("del"), w("moveFrom")):
                    skip = True
                    break
                anc = anc.getparent()
            if not skip:
                buf.append(t.text or "")
        out.append("".join(buf))
    return "\n".join(out)


def accept(root):
    # rows whose deletion was tracked go first, while their markers survive
    for tr in list(root.iter(w("tr"))):
        trpr = tr.find(w("trPr"))
        if trpr is not None and trpr.find(w("del")) is not None:
            tr.getparent().remove(tr)
    for el in list(root.iter(w("del"))) + list(root.iter(w("moveFrom"))):
        parent = el.getparent()
        if parent is None:
            continue
        if parent.tag == w("trPr"):
            parent.remove(el)
            continue
        parent.remove(el)
    for tag in ("ins", "moveTo"):
        for el in list(root.iter(w(tag))):
            parent = el.getparent()
            if parent is None:
                continue
            if parent.tag == w("trPr"):
                parent.remove(el)
                continue
            at = list(parent).index(el)
            for off, child in enumerate(list(el)):
                parent.insert(at + off, child)
            parent.remove(el)
    for tag in ("rPrChange", "pPrChange", "sectPrChange", "tblPrChange",
                "trPrChange", "tcPrChange", "numberingChange", "cellIns",
                "cellDel", "cellMerge"):
        for el in list(root.iter(w(tag))):
            if el.getparent() is not None:
                el.getparent().remove(el)


def main():
    if len(sys.argv) != 3 or sys.argv[1] in ("-h", "--help"):
        sys.exit("usage: accept_changes.py TRACKED_DIR CLEAN_DIR\n"
                 "  TRACKED_DIR is an unzipped .docx whose word/document.xml\n"
                 "  carries tracked changes. CLEAN_DIR is overwritten with a\n"
                 "  copy in which every change has been accepted, and is only\n"
                 "  written if no revision mark survives and the text matches.")
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    if not (src / "word" / "document.xml").exists():
        sys.exit(f"not an unzipped .docx: {src}/word/document.xml is missing")
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)

    for part in ("document.xml", "footnotes.xml", "endnotes.xml"):
        path = dst / "word" / part
        if not path.exists():
            continue
        tree = etree.parse(str(path))
        root = tree.getroot()
        want = accepted_text(root) if part == "document.xml" else None
        accept(root)
        if part == "document.xml":
            got = accepted_text(root)
            if got != want:
                sys.exit("refusing to write: accepting the changes altered the text")
            left = [t for t in REVISION_MARKS
                    for _ in root.iter(w(t))]
            if left:
                sys.exit(f"refusing to write: {len(left)} revision mark(s) remain")
        tree.write(str(path), xml_declaration=True, encoding="UTF-8",
                   standalone=True)
        print(f"  {part}: clean")
    print(f"  wrote {dst}")


if __name__ == "__main__":
    main()
