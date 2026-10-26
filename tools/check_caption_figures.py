#!/usr/bin/env python3
"""check_caption_figures.py — check every figure stated in a new caption.

A caption is prose, so the cell-by-cell table check cannot see it. This round
shipped a wrong model length range in a caption for exactly that reason.

Matching a caption number against the pool of all released values does not
work and was tried first. With roughly thirty thousand released values, any
small integer matches something, so the check passed on the wrong figure too.
A pool match has a resolution limit and bare integers are below it.

What works is naming the source. Each entry below states a figure as the
caption prints it and how to derive it from the release. The two must agree,
and the figure must also be present in the built document, so a caption cannot
drift from the analysis and the gate cannot pass on a caption that was edited
away.

USAGE
    python3 tools/check_caption_figures.py [path to the built .docx]
"""
import csv
import os
import re
import sys
import zipfile

from lxml import etree

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
REL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT = ("/mnt/user-data/outputs/r17/"
           "Supplementary_Information_MG_revised_clean.docx")


def tsv(path):
    with open(os.path.join(REL, path), newline="", encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if not list(r.values())[0].startswith("#")]


def load():
    return {
        "prox": tsv("results/tables/breakpoint_proximity_meropenem_clsi_bmd.tsv"),
        "sens": tsv("results/tables/breakpoint_sensitivity_meropenem_clsi_bmd.tsv"),
        "nest": tsv("results/tables/provenance_nesting_meropenem_clsi_bmd.tsv"),
        "fam": tsv("results/tables/family_within_provenance_meropenem_clsi_bmd.tsv"),
        "floor": tsv("results/tables/allele_search_detection_floor.tsv"),
        "summ": tsv("results/tables/profile_screen_summary.tsv"),
        "unm": tsv("results/tables/profile_screen_unmatched_proteins.tsv"),
        "prf": tsv("results/allele_search/ref/betalactamase_families.tsv"),
        "adj": tsv("results/tables/table4_determinants_adjusted.tsv"),
    }


def checks(d):
    prox, sens, nest = d["prox"], d["sens"], d["nest"]
    fam, floor = d["fam"], d["floor"]
    summ, unm = d["summ"], d["unm"]
    prf = {r["family"]: r for r in d["prf"]}
    adj = {r["term"]: r for r in d["adj"]}
    allpos = next(r for r in sens if r["group"] == "all carbapenemase-positive")
    spec = [r for r in unm if r["profile_specific"] == "True"]
    broad = [r for r in unm if r["profile_specific"] == "False"]
    f15 = {r["family"]: r for r in fam if r["reporting_year"] == "2015"}

    # (table, the figure as the caption prints it, the value derived from the
    #  release, where it came from)
    out = [
        ("S20a", "124", allpos["susceptible_as_published"],
         "breakpoint_sensitivity, susceptible_as_published"),
        ("S20a", "1",
         next(r["mic_upper_mgL"] for r in prox
              if r["at_or_above_breakpoint"] == "True"),
         "breakpoint_proximity, the row at or above the breakpoint"),

        ("S21b", "ten", "10",
         "the stated minimum family size, a rule not a result"),
        ("S21b", "two", "2",
         "the stated minimum number of families, a rule"),

        ("S22", "Ninety", str(
            sum(1 for r in floor
                if abs(float(r["target_identity"]) - 0.95) < 1e-9)),
         "allele_search_detection_floor, rows at the first rung"),
        ("S22", "five", str(len({r["replicate"] for r in floor})),
         "allele_search_detection_floor, distinct replicates"),
        ("S22", "three", str(len({r["model"] for r in floor})),
         "allele_search_detection_floor, distinct divergence models"),
        ("S22", "90",
         str(round(100 * min(float(r["target_identity"]) for r in floor
                             if r["verdict"] == "novel_candidate"))),
         "the lowest rung at which any mutant is still called a new allele"),
        ("S22", "70",
         str(round(100 * min(
             float(r["target_identity"]) for r in floor
             if r["verdict"] == "divergent_homologue"))),
         "the lowest rung at which any mutant is still a divergent homologue"),
        ("S22", "42.8", max(
            r["background_ceiling_max_pct"] for r in
            tsv("results/supplementary/S19a_allele_search_summary.tsv")),
         "S19a_allele_search_summary, the highest background ceiling"),

        ("S23", "136", str(len(prf)), "betalactamase_families, profiles built"),
        ("S23", "55", str(sum(1 for v in prf.values()
                              if int(v["carbapenemase_members"]) > 0)),
         "betalactamase_families, profiles containing a carbapenemase"),
        ("S23", "100", "100", "the stated bit-score threshold, a rule"),
        ("S23", "1.25", "1.25", "the stated specificity threshold, a rule"),
        ("S23", "135", str(sum(1 for v in prf.values()
                               if v["profile_specific"] == "True")),
         "betalactamase_families, profiles flagged specific"),
        ("S23", "26", prf["blaOXA-PR"]["members"],
         "betalactamase_families, blaOXA-PR members"),
        ("S23", "265", prf["blaOXA-PR"]["member_aa_min"],
         "betalactamase_families, blaOXA-PR shortest member"),
        ("S23", "610", prf["blaOXA-PR"]["member_aa_max"],
         "betalactamase_families, blaOXA-PR longest member"),
        ("S23", "44", str(len(broad)),
         "profile_screen_unmatched_proteins, hits from a broad profile"),
        ("S23", "73", str(len(unm)),
         "profile_screen_unmatched_proteins, all unaccounted hits"),
        ("S23", "29", str(len(spec)),
         "profile_screen_unmatched_proteins, hits from a specific profile"),
        ("S23", "298", str(sum(int(r["genomes"]) for r in summ)),
         "profile_screen_summary, genomes screened"),
        ("S23", "63", str(min(int(r["aa_length"]) for r in spec)),
         "the shortest specific-profile unaccounted protein"),
        ("S23", "168", str(max(int(r["aa_length"]) for r in spec)),
         "the longest specific-profile unaccounted protein"),
        ("S23", "266", str(min(int(prf[r["profile_family"]]["model_length"])
                               for r in spec)),
         "the shortest model behind a specific-profile unaccounted hit"),
        ("S23", "381", str(max(int(prf[r["profile_family"]]["model_length"])
                               for r in spec)),
         "the longest model behind a specific-profile unaccounted hit"),
        ("S23", "three", "3", "the stated minimum family size, a rule"),

        ("S24", "Table 4", "Table 4", "a cross-reference, not a result"),
    ]
    return out


def doc_text(path):
    if path.endswith(".docx"):
        xml = zipfile.ZipFile(path).read("word/document.xml")
    else:
        xml = open(path, "rb").read()
    root = etree.fromstring(xml)
    body = root.find(W + "body")
    # Only body-level paragraphs. A paragraph inside a table cell is also a
    # w:p, and walking every w:p would mix table contents into the caption.
    seq = []
    for el in body:
        if el.tag == W + "p":
            seq.append(("p", " ".join(
                "".join(e.text or "" for e in el.iter(W + "t")).split())))
        elif el.tag == W + "tbl":
            seq.append(("tbl", ""))
    caps = {}
    for i, (kind, t) in enumerate(seq):
        if kind != "p":
            continue
        m = re.match(r"Table (S\d+[ab]?)\.", t)
        if not m:
            continue
        # A journal caption is often a short title followed by one or more
        # prose paragraphs carrying the explanation, labelled Notes or not.
        # Everything up to the table itself, or up to the next caption, is
        # part of the caption. Without this the gate reports every figure
        # absent on a document whose captions are split that way.
        parts = [t]
        j = i + 1
        while j < len(seq) and seq[j][0] == "p":
            nxt = seq[j][1]
            if re.match(r"Table S\d+[ab]?\.", nxt) or not nxt:
                break
            parts.append(nxt)
            j += 1
        caps[m.group(1)] = " ".join(parts)
    return caps


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    caps = doc_text(path)
    d = load()
    rows = checks(d)
    bad = []
    missing = []
    for table, printed, derived, where in rows:
        if str(printed).replace(",", "") != str(derived).replace(",", ""):
            # a word form stands for its number
            words = {"Ninety": "90", "five": "5", "three": "3",
                     "two": "2", "ten": "10"}
            if words.get(str(printed)) != str(derived).replace(",", ""):
                bad.append((table, printed, derived, where))
                continue
        cap = caps.get(table)
        if cap is None:
            missing.append((table, printed, "the caption is not in the document"))
        elif str(printed) not in cap:
            missing.append((table, printed, "the figure is not in the caption"))
    print(f"  {len(rows)} figure(s) checked across {len({r[0] for r in rows})} "
          f"caption(s)")
    for table, printed, derived, where in bad:
        print(f"    MISMATCH Table {table}: caption says {printed!r}, "
              f"the release gives {derived!r} ({where})")
    for table, printed, why in missing:
        print(f"    ABSENT   Table {table}: {printed!r} {why}")
    if bad or missing:
        sys.exit(f"{len(bad)} mismatch(es) and {len(missing)} absent figure(s)")
    print("  every figure stated in a new caption matches its named source")


if __name__ == "__main__":
    main()
