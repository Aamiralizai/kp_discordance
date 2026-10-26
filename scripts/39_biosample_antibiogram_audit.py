#!/usr/bin/env python3
"""39_biosample_antibiogram_audit.py — audit the external cohort against primary records.

Pathogen Detection is a secondary aggregator. It restates susceptibility results
that submitters deposited on the BioSample record, so the BioSample antibiogram
is the primary source for every external isolate. This script retrieves that
antibiogram for every external isolate and compares every drug it reports with
the call Pathogen Detection carries for the same isolate and drug.

The comparison is made on the deposited phenotype alone, across all drugs
reported. It never reads the carbapenemase genotype, so it cannot be tuned by
the genotype-phenotype agreement this study measures, and a submission whose
record is sound is left untouched.

Where the two sources disagree the BioSample call is the one retained, because
it is the record the submitting laboratory deposited.

Inputs. The deposited BioSample records are cached under data/biosample_xml, so
that side of the comparison needs no network. The Pathogen Detection snapshot is
too large to redistribute and is fetched by scripts/19b_fetch_ncbi_pathogen.sh.

Output
    results/tables/external_biosample_audit.tsv          one row per isolate
    results/tables/external_biosample_concordance.tsv    per BioProject
    results/tables/external_biosample_audit_summary.tsv  cohort totals

USAGE
    python3 scripts/39_biosample_antibiogram_audit.py
"""
from __future__ import annotations

import collections
import csv
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

REL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CACHE = os.path.join(REL, "data", "biosample_xml")
BATCH = 40
PAUSE = 0.4
RETRY = 4

# A submission is treated as carrying an inverted AST record only on this
# evidence: every comparable call disagrees, over at least this many calls and
# this many drugs. One drug or a handful of calls would not distinguish a
# systematic fault from ordinary disagreement between a submitter and an
# aggregator, so the bar is set where chance is not a competing explanation.
INVERSION_MIN_CALLS = 100
INVERSION_MIN_DRUGS = 5


# Resistant and susceptible MIC breakpoints in mg/L for Enterobacterales, used
# only to check the deposited numeric MIC against the deposited phenotype label.
# They are not used to classify anything in the study. EUCAST 2024 and CLSI
# M100 agree on the direction for both drugs, which is all this check needs.
BREAKPOINTS = {
    ("meropenem", "EUCAST"): (2.0, 8.0),     # S <= 2, R > 8
    ("meropenem", "CLSI"): (1.0, 4.0),       # S <= 1, R >= 4
    ("ampicillin", "EUCAST"): (8.0, 8.0),    # S <= 8, R > 8
    ("ampicillin", "CLSI"): (8.0, 32.0),     # S <= 8, R >= 32
}

# The antibiogram table columns this script reads, by the name BioSample uses.
# Asserting them catches a field-mapping error, where the right value is read
# from the wrong column because a submitter ordered the table differently.
REQUIRED_COLUMNS = ("antibiotic", "resistance phenotype", "measurement sign",
                    "measurement", "measurement units")


def mic_call(drug: str, sign: str, value: str, units: str,
             standard: str) -> str:
    """S or R from the deposited MIC alone, or "" when it cannot be decided.

    This reads the numeric measurement, which is a different field from the
    phenotype label and is parsed by different code. If the label parsing were
    at fault, the two would not agree.

    A censored measurement is treated as the interval it denotes, and a call is
    returned only when the whole interval falls on one side of a breakpoint.
    """
    if (units or "").strip().lower() not in ("mg/l", "ug/ml", "mcg/ml", ""):
        return ""
    bp = BREAKPOINTS.get((drug, (standard or "EUCAST").strip().upper()))
    if bp is None:
        bp = BREAKPOINTS.get((drug, "EUCAST"))
    if bp is None:
        return ""
    s_max, r_min = bp
    raw = (value or "").strip()
    # a combination such as "32/2" is reported as the first component
    raw = raw.split("/")[0].strip()
    try:
        x = float(raw)
    except ValueError:
        return ""
    sg = (sign or "==").strip()
    if sg in (">", ">="):
        # the true value is at or above x, so only a high x decides
        return "R" if x >= r_min else ""
    if sg in ("<", "<="):
        return "S" if x <= s_max else ""
    if x >= r_min:
        return "R"
    if x <= s_max:
        return "S"
    return ""


def norm(s: str) -> str:
    """S, I or R from a BioSample resistance-phenotype string."""
    s = (s or "").strip().lower()
    if s.startswith("susceptible"):
        return "I" if "dose" in s else "S"
    if s.startswith("intermediate"):
        return "I"
    if s.startswith("resistant") or s.startswith("nonsusceptible"):
        return "R"
    return ""


def fetch(ids: list[str]) -> str:
    """One efetch call, cached on disk so a rerun is offline and identical."""
    os.makedirs(CACHE, exist_ok=True)
    key = os.path.join(CACHE, f"{ids[0]}_{len(ids)}.xml")
    if os.path.exists(key) and os.path.getsize(key) > 0:
        with open(key, encoding="utf-8") as fh:
            return fh.read()
    data = urllib.parse.urlencode(
        {"db": "biosample", "id": ",".join(ids),
         "rettype": "full", "retmode": "xml"}).encode()
    last = None
    for attempt in range(RETRY):
        try:
            with urllib.request.urlopen(
                    urllib.request.Request(EUTILS, data=data), timeout=180) as r:
                txt = r.read().decode("utf-8", "replace")
            if "<BioSampleSet" not in txt:
                raise ValueError("no BioSampleSet in the reply")
            with open(key, "w", encoding="utf-8") as fh:
                fh.write(txt)
            return txt
        except (urllib.error.URLError, ValueError, TimeoutError) as exc:
            last = exc
            time.sleep(2 ** attempt)
    sys.exit(f"could not retrieve {len(ids)} BioSample record(s): {last}")


def parse(xml: str) -> dict[str, dict[str, tuple]]:
    """accession -> {drug: (call, sign, value, units, standard)}"""
    out: dict[str, dict[str, tuple]] = {}
    for bs in ET.fromstring(xml).findall("BioSample"):
        acc = bs.get("accession") or ""
        rec: dict[str, tuple] = {}
        for tbl in bs.iter("Table"):
            if not (tbl.get("class") or "").startswith("Antibiogram"):
                continue
            header = tbl.find("Header")
            body = tbl.find("Body")
            if header is None or body is None:
                continue
            head = [(c.text or "").strip().lower()
                    for c in header.findall("Cell")]
            missing = [c for c in REQUIRED_COLUMNS if c not in head]
            if missing:
                sys.exit(
                    f"{acc}: its antibiogram table does not name the "
                    f"column(s) {missing}. The header is {head}. Reading it "
                    "positionally could map a value to the wrong field, so "
                    "the audit stops rather than guess.")
            widths = {len(r.findall("Cell")) for r in body.findall("Row")}
            if widths - {len(head)}:
                sys.exit(
                    f"{acc}: its antibiogram rows have widths {sorted(widths)} "
                    f"against a header of {len(head)} column(s), so the "
                    "column mapping is not reliable.")
            for row in body.findall("Row"):
                d = dict(zip(head, [(c.text or "").strip()
                                    for c in row.findall("Cell")]))
                drug = (d.get("antibiotic") or "").strip().lower()
                call = norm(d.get("resistance phenotype", ""))
                if not drug or not call:
                    continue
                rec[drug] = (call,
                             d.get("measurement sign", ""),
                             d.get("measurement", ""),
                             d.get("measurement units", ""),
                             d.get("testing standard", ""))
        out[acc] = rec
    return out


def tsv(path: str) -> list[dict]:
    with open(os.path.join(REL, path), newline="", encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if not list(r.values())[0].startswith("#")]


PD_SNAPSHOT = "data/external/ncbi_kp_ast.tsv"


def require_snapshot() -> str:
    """Checked before any retrieval, so a missing input fails in a second."""
    path = os.path.join(REL, PD_SNAPSHOT)
    if not os.path.exists(path):
        sys.exit(
            f"The Pathogen Detection snapshot is absent:\n  {path}\n"
            "It is too large to redistribute, so the repository does not ship "
            "it. Retrieve it first with\n"
            "  bash scripts/19b_fetch_ncbi_pathogen.sh\n"
            "The BioSample side of the audit is cached under data/biosample_xml "
            "and needs no network.")
    return path


def pd_calls(path: str) -> dict[str, dict[str, str]]:
    """accession -> {drug: S|I|R} as Pathogen Detection carries it."""
    csv.field_size_limit(10 ** 7)
    out: dict[str, dict[str, str]] = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            acc = rec.get("biosample_acc", "")
            if not acc:
                continue
            d: dict[str, str] = {}
            for tok in (rec.get("AST_phenotypes") or "").split(","):
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    c = v.strip().upper()[:1]
                    if c in ("S", "I", "R"):
                        d[k.strip().lower()] = c
            if d:
                out[acc] = d
    return out


def write(path: str, rows: list[dict]) -> str:
    full = os.path.join(REL, path)
    with open(full, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return full


def main() -> None:
    pd_path = require_snapshot()
    iso = tsv("results/tables/external_isolates.tsv")
    accs = sorted({r["biosample_acc"] for r in iso if r["biosample_acc"]})
    print(f"external cohort: {len(iso)} isolate(s), "
          f"{len(accs)} distinct BioSample accession(s)")

    bs: dict[str, dict[str, tuple]] = {}
    for i in range(0, len(accs), BATCH):
        bs.update(parse(fetch(accs[i:i + BATCH])))
        done = min(i + BATCH, len(accs))
        if done % 200 == 0 or done == len(accs):
            print(f"  retrieved {done}/{len(accs)}")
        time.sleep(PAUSE)
    pdc = pd_calls(pd_path)

    # ---- per-BioProject concordance over every drug both sources report ----
    per = collections.defaultdict(
        lambda: {"n_iso": 0, "n_calls": 0, "n_conc": 0,
                 "drugs": set(), "iso_with_ab": 0})
    for r in iso:
        p = per[r["bioproject"]]
        p["n_iso"] += 1
        acc = r["biosample_acc"]
        bd = bs.get(acc) or {}
        if bd:
            p["iso_with_ab"] += 1
        pdd = pdc.get(acc) or {}
        for drug, tup in bd.items():
            got = pdd.get(drug)
            if got in ("S", "R") and tup[0] in ("S", "R"):
                p["n_calls"] += 1
                p["drugs"].add(drug)
                if got == tup[0]:
                    p["n_conc"] += 1

    conc_rows, inverted = [], set()
    for bp, p in sorted(per.items(), key=lambda kv: -kv[1]["n_calls"]):
        n, c, nd = p["n_calls"], p["n_conc"], len(p["drugs"])
        flag = (n >= INVERSION_MIN_CALLS and nd >= INVERSION_MIN_DRUGS
                and c == 0)
        if flag:
            inverted.add(bp)
        conc_rows.append({
            "bioproject": bp,
            "n_isolates": p["n_iso"],
            "n_with_biosample_antibiogram": p["iso_with_ab"],
            "n_comparable_calls": n,
            "n_drugs_compared": nd,
            "n_concordant": c,
            "pct_concordant": f"{100.0 * c / n:.1f}" if n else "",
            "record_inverted": int(flag),
        })

    # ---- one row per isolate, with the phenotype taken from the primary record
    rows, changed = [], 0
    for r in iso:
        acc = r["biosample_acc"]
        bd = bs.get(acc) or {}
        mem = bd.get("meropenem", ("", "", "", "", ""))
        amp = bd.get("ampicillin", ("", "", "", "", ""))
        # Script 19 overwrites `call` with the corrected phenotype once this
        # audit exists. Comparing against that column would compare the
        # correction with itself and report perfect agreement, so the
        # as-deposited call is the one read whenever 19 has recorded it.
        pd_mem = (r.get("call_as_deposited") or r.get("call") or "").strip()
        primary = mem[0] if mem[0] in ("S", "I", "R") else pd_mem
        if primary != pd_mem:
            changed += 1
        rows.append({
            "biosample_acc": acc,
            "asm_acc": r["asm_acc"],
            "bioproject": r["bioproject"],
            "pd_meropenem": pd_mem,
            "biosample_meropenem": mem[0],
            "biosample_meropenem_mic": (mem[1] + mem[2]).strip(),
            "biosample_meropenem_units": mem[3],
            "pd_ampicillin": r.get("intrinsic_drug_call", ""),
            "biosample_ampicillin": amp[0],
            "biosample_ampicillin_mic": (amp[1] + amp[2]).strip(),
            "testing_standard": mem[4] or amp[4],
            "has_biosample_antibiogram": int(bool(bd)),
            "record_inverted_submission": int(r["bioproject"] in inverted),
            "meropenem_primary": primary,
            "meropenem_source": ("biosample" if mem[0] in ("S", "I", "R")
                                 else "pathogen_detection"),
            "meropenem_disagrees": int(bool(mem[0]) and bool(pd_mem)
                                       and mem[0] != pd_mem),
        })

    # ---------------- parsing and field-mapping control ----------------
    # The inversion is asserted on the phenotype label. If that label were
    # being parsed or mapped wrongly, the deposited numeric MIC would not
    # support it. So each call is re-derived from the MIC alone, by different
    # code reading a different field, and compared with both sources.
    ctrl_rows = []
    for r in iso:
        acc = r["biosample_acc"]
        bd = bs.get(acc) or {}
        pdd = pdc.get(acc) or {}
        for drug in ("meropenem", "ampicillin"):
            tup = bd.get(drug)
            if not tup:
                continue
            label, sign, value, units, std = tup
            derived = mic_call(drug, sign, value, units, std)
            got = pdd.get(drug)
            if not derived or label not in ("S", "R") or got not in ("S", "R"):
                continue
            ctrl_rows.append({
                "biosample_acc": acc,
                "bioproject": r["bioproject"],
                "drug": drug,
                "mic_sign": sign,
                "mic_value": value,
                "mic_units": units,
                "testing_standard": std,
                "mic_derived_call": derived,
                "biosample_label": label,
                "pathogen_detection_label": got,
                "mic_supports_biosample": int(derived == label),
                "mic_supports_pathogen_detection": int(derived == got),
                "record_inverted_submission": int(r["bioproject"] in inverted),
            })

    def control(rows, label):
        n = len(rows)
        a = sum(x["mic_supports_biosample"] for x in rows)
        b = sum(x["mic_supports_pathogen_detection"] for x in rows)
        return {
            "cohort": label, "n_calls_with_decidable_mic": n,
            "mic_agrees_with_biosample": a,
            "pct_agrees_with_biosample": f"{100.0 * a / n:.1f}" if n else "",
            "mic_agrees_with_pathogen_detection": b,
            "pct_agrees_with_pathogen_detection": (
                f"{100.0 * b / n:.1f}" if n else ""),
        }

    ctrl_summary = [
        control(ctrl_rows, "all external"),
        control([x for x in ctrl_rows if x["record_inverted_submission"]],
                "inverted submission"),
        control([x for x in ctrl_rows if not x["record_inverted_submission"]],
                "all other submissions"),
    ]
    if ctrl_rows:
        write("results/tables/external_biosample_parsing_control.tsv", ctrl_rows)
        write("results/tables/external_biosample_parsing_control_summary.tsv",
              ctrl_summary)

    tot_calls = sum(r["n_comparable_calls"] for r in conc_rows)
    tot_conc = sum(r["n_concordant"] for r in conc_rows)
    with_ab = sum(r["has_biosample_antibiogram"] for r in rows)
    mem_cmp = [r for r in rows if r["biosample_meropenem"] and r["pd_meropenem"]]
    mem_bad = [r for r in mem_cmp if r["meropenem_disagrees"]]
    inv_iso = [r for r in rows if r["record_inverted_submission"]]

    def pct(a, b):
        return f"{100.0 * a / b:.1f}" if b else ""

    summary = [
        ("external isolates", len(rows), ""),
        ("with a BioSample antibiogram", with_ab, pct(with_ab, len(rows))),
        ("comparable susceptibility calls, all drugs", tot_calls, ""),
        ("concordant with Pathogen Detection", tot_conc,
         pct(tot_conc, tot_calls)),
        ("meropenem comparable", len(mem_cmp), pct(len(mem_cmp), len(rows))),
        ("meropenem disagrees", len(mem_bad), pct(len(mem_bad), len(mem_cmp))),
        ("submissions with a wholly inverted record", len(inverted), ""),
        ("isolates in those submissions", len(inv_iso),
         pct(len(inv_iso), len(rows))),
        ("meropenem calls taken from the BioSample record", changed,
         pct(changed, len(rows))),
    ]
    sp = os.path.join(REL,
                      "results/tables/external_biosample_audit_summary.tsv")
    with open(sp, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["quantity", "n", "pct"])
        for a, b, c in summary:
            w.writerow([a, b, c])

    # The corrected records on their own, so the evidence behind the
    # correction can be read without filtering the full cohort.
    corrected_rows = [r for r in rows if r["meropenem_disagrees"]]
    if corrected_rows:
        for r in corrected_rows:
            r = r
        write("results/tables/external_biosample_corrected_records.tsv",
              [{**r,
                "reason": ("the deposited antibiogram contradicts the "
                           "aggregated record, in a submission whose every "
                           "comparable call is contradicted")}
               for r in corrected_rows])

    a1 = write("results/tables/external_biosample_audit.tsv", rows)
    a2 = write("results/tables/external_biosample_concordance.tsv", conc_rows)

    for a, b, c in summary:
        print(f"  {a:46s} {b:>6}  {c}")
    print("\nparsing and field-mapping control, the call re-derived from the "
          "deposited MIC alone")
    for c in ctrl_summary:
        print(f"  {c['cohort']:24s} n={c['n_calls_with_decidable_mic']:5d}  "
              f"agrees with BioSample "
              f"{c['pct_agrees_with_biosample']:>6s}%  "
              f"with Pathogen Detection "
              f"{c['pct_agrees_with_pathogen_detection']:>6s}%")
    print(f"\n  submissions flagged as inverted: "
          f"{', '.join(sorted(inverted)) or 'none'}")
    print("\nper-BioProject concordance (largest first)")
    for r in conc_rows[:12]:
        print(f"  {r['bioproject']:16s} iso={r['n_isolates']:4d} "
              f"calls={r['n_comparable_calls']:5d} "
              f"drugs={r['n_drugs_compared']:3d} "
              f"concordant={r['pct_concordant']:>6s}%"
              f"{'   INVERTED' if r['record_inverted'] else ''}")
    for p in (a1, a2, sp):
        print(f"wrote {p}")


if __name__ == "__main__":
    main()
