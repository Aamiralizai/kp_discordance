#!/usr/bin/env python
"""kp_defs.py — single source of truth for genotype definitions.

Imported by every analysis script so that carbapenemase family assignment and
porin classification cannot drift between analyses (a defect in earlier
versions, where OXA-162 was included in the OXA-48-like family by one script
and omitted by others).

Import as:
    sys.path.insert(0, str(Path(__file__).parent))
    from kp_defs import carb_family, carb_allele, porin_status, has_carb
"""
import re

# ---------------------------------------------------------------- families
# Curated allele -> family dictionary. Matching is on the cleaned allele token.
# Extend deliberately; do not add ad-hoc substring rules in analysis scripts.
OXA48_LIKE = {
    "OXA-48", "OXA-162", "OXA-181", "OXA-204", "OXA-232", "OXA-244",
    "OXA-245", "OXA-247", "OXA-370", "OXA-405", "OXA-484", "OXA-517",
}

FAMILY_PREFIX = [
    ("KPC", re.compile(r"^KPC[-_]", re.I)),
    ("NDM", re.compile(r"^NDM[-_]", re.I)),
    ("VIM", re.compile(r"^VIM[-_]", re.I)),
    ("IMP", re.compile(r"^IMP[-_]", re.I)),
    ("GES", re.compile(r"^GES[-_]", re.I)),
    ("SPM", re.compile(r"^SPM[-_]", re.I)),
]

EMPTY = {"", "-", "nan", "NA", "None", "none"}

# Kleborate appends confidence/context flags to allele strings; strip them.
_FLAGS = re.compile(r"[\*\?\^\$\~\!]+")


def has_carb(v) -> bool:
    """Any acquired carbapenemase reported?"""
    return isinstance(v, str) and v.strip() not in EMPTY


def _tokens(v):
    if not has_carb(v):
        return []
    out = []
    for tok in re.split(r"[;,]", v):
        tok = _FLAGS.sub("", tok).strip()
        if tok and tok not in EMPTY:
            out.append(tok)
    return out


def carb_allele(v) -> str:
    """Cleaned, sorted, semicolon-joined allele string ('none' if absent)."""
    toks = _tokens(v)
    return ";".join(sorted(set(toks))) if toks else "none"


def allele_family(tok: str) -> str:
    """Family for a single cleaned allele token."""
    t = tok.strip().upper()
    if t in {a.upper() for a in OXA48_LIKE}:
        return "OXA-48-like"
    for fam, rx in FAMILY_PREFIX:
        if rx.match(t):
            return fam
    # OXA alleles not in the curated OXA-48-like set are reported separately
    # rather than silently folded into that family.
    if t.startswith("OXA-"):
        return "OXA-other"
    return "other"


def carb_family(v) -> str:
    """Family label for an isolate. Isolates carrying alleles of more than one
    family receive a combined label (e.g. 'KPC+NDM') and are analysed apart."""
    toks = _tokens(v)
    if not toks:
        return "none"
    fams = sorted({allele_family(t) for t in toks})
    return "+".join(fams)


# ---------------------------------------------------------------- porins
# Recognised loop-3 di-amino-acid insertion positions in OmpK36 (HGVS protein
# numbering as emitted by Kleborate v3). Insertions elsewhere in the protein
# are NOT loop-3 events and are classified separately.
LOOP3_POSITIONS = {(134, 135), (135, 136), (136, 137)}
LOOP3_RESIDUES = {"GLYASP", "THRASP", "ASPTHR", "ASPGLY"}   # GD / TD

_TRUNCATING = re.compile(r"(fs|Ter|\*)", re.I)          # frameshift / nonsense
_INFRAME_DEL = re.compile(r"\bdel\b|del[A-Za-z]{3}", re.I)
_INS = re.compile(r"p\.(\d+)_(\d+)ins([A-Za-z]+)", re.I)
_NT_DEL = re.compile(r"c\.[A-Z]?\d+del", re.I)
_C25CT = re.compile(r"c\.25C>T", re.I)

PORIN_ORDER = ["none", "ompK35_loss", "ompK36_reduced_expression",
               "ompK36_loop3", "ompK36_other_variant", "ompK36_loss"]

PORIN_LABEL = {
    "none": "intact",
    "ompK35_loss": "OmpK35 truncation",
    "ompK36_reduced_expression": "OmpK36 c.25C>T (reduced expression)",
    "ompK36_loop3": "OmpK36 loop-3 insertion",
    "ompK36_other_variant": "OmpK36 other in-frame variant",
    "ompK36_loss": "OmpK36 truncation",
}


def porin_status(v, collapse=False):
    """Classify OmpK35/OmpK36 changes from Kleborate v3 HGVS notation.

    Categories (ordered by expected impact on carbapenem influx):
      none
      ompK35_loss                 OmpK35 frameshift / nonsense / nt deletion
      ompK36_reduced_expression   OmpK36 c.25C>T — SYNONYMOUS; reduces OmpK36
                                  abundance post-transcriptionally via an
                                  inhibitory mRNA secondary structure. This is
                                  NOT a premature stop codon.
      ompK36_loop3                di-amino-acid insertion at a recognised
                                  loop-3 position (GD / TD)
      ompK36_other_variant        in-frame insertion/deletion elsewhere in
                                  OmpK36 — impact not established
      ompK36_loss                 OmpK36 frameshift / nonsense

    OmpK36 categories take precedence over OmpK35. Within OmpK36, truncation >
    loop-3 insertion > reduced expression > other variant.

    collapse=True merges ompK36_reduced_expression and ompK36_other_variant
    into ompK36_loss and ompK36_loop3 respectively, reproducing the coarser
    scheme used in earlier analyses (for comparison only).
    """
    if not isinstance(v, str) or v.strip() in EMPTY:
        return "none"

    k36_trunc = k36_loop3 = k36_other = k36_lowexp = k35_trunc = False

    for tok in re.split(r"[;,]", v):
        tok = tok.strip()
        m = re.match(r"(OmpK3[56])\s*:\s*(.+)$", tok, flags=re.I)
        if not m:
            continue
        gene, chg = m.group(1).upper(), m.group(2).strip()

        truncating = bool(_TRUNCATING.search(chg)) or bool(_NT_DEL.search(chg))

        if gene == "OMPK36":
            if _C25CT.search(chg):
                k36_lowexp = True
                continue
            if truncating:
                k36_trunc = True
                continue
            mi = _INS.search(chg)
            if mi:
                pos = (int(mi.group(1)), int(mi.group(2)))
                res = mi.group(3).upper()
                if pos in LOOP3_POSITIONS and res in LOOP3_RESIDUES:
                    k36_loop3 = True
                else:
                    k36_other = True
                continue
            if _INFRAME_DEL.search(chg):
                k36_other = True
                continue
            k36_other = True
        elif gene == "OMPK35":
            if truncating:
                k35_trunc = True

    if k36_trunc:
        cls = "ompK36_loss"
    elif k36_loop3:
        cls = "ompK36_loop3"
    elif k36_lowexp:
        cls = "ompK36_reduced_expression"
    elif k36_other:
        cls = "ompK36_other_variant"
    elif k35_trunc:
        cls = "ompK35_loss"
    else:
        cls = "none"

    if collapse:
        cls = {"ompK36_reduced_expression": "ompK36_loss",
               "ompK36_other_variant": "ompK36_loop3"}.get(cls, cls)
    return cls


def porin_disrupted(cls: str) -> bool:
    """OmpK36-level disruption (the mechanistically potent category)."""
    return cls in {"ompK36_loop3", "ompK36_loss",
                   "ompK36_reduced_expression", "ompK36_other_variant"}


# ---------------------------------------------------------------- breakpoints
# Breakpoints are loaded from config/breakpoints.csv and are USED ONLY IF the
# row is marked verified=YES. This is a hard gate: an unverified breakpoint
# previously produced a fabricated discordance estimate and an apparent
# reversal of a central finding. Guessing is not permitted here.
import csv as _csv
from pathlib import Path as _Path

_BP_FILE = _Path(__file__).resolve().parent.parent / "config" / "breakpoints.csv"


def _load_breakpoints():
    table, unverified = {}, []
    if not _BP_FILE.exists():
        return table, unverified
    with open(_BP_FILE) as fh:
        for row in _csv.DictReader(
                r for r in fh if not r.lstrip().startswith("#")):
            if not row.get("drug"):
                continue
            drug = row["drug"].strip().lower()
            std = row["standard"].strip().upper()
            entry = dict(S=float(row["S_value_mgL"]), R=float(row["R_value_mgL"]),
                         S_op=row.get("S_operator", "<=").strip() or "<=",
                         R_op=row.get("R_operator", ">=").strip() or ">=",
                         source=row.get("source_document", "").strip(),
                         notes=row.get("notes", "").strip())
            if row.get("verified", "").strip().upper() == "YES":
                table.setdefault(std, {})[drug] = entry
            else:
                unverified.append((drug, std))
    return table, unverified


BREAKPOINTS, UNVERIFIED_BREAKPOINTS = _load_breakpoints()


def get_breakpoints(drug, standard):
    """Return verified breakpoints, or raise with an actionable message."""
    d, s = drug.strip().lower(), standard.strip().upper()
    if s in BREAKPOINTS and d in BREAKPOINTS[s]:
        return BREAKPOINTS[s][d]
    if (d, s) in UNVERIFIED_BREAKPOINTS:
        raise SystemExit(
            f"\nREFUSING TO ANALYSE {drug} under {standard}.\n"
            f"The breakpoint row exists in {_BP_FILE} but is marked verified=NO.\n"
            f"Open the {standard} document you intend to cite, confirm the\n"
            f"Enterobacterales breakpoints for {drug}, record the edition in\n"
            f"source_document, and set verified=YES.\n"
            f"Analysing with an unchecked breakpoint has already produced one\n"
            f"fabricated result in this project.\n")
    raise SystemExit(
        f"\nNo breakpoint row for {drug} under {standard} in {_BP_FILE}.\n"
        f"Add one, verify it against the source document, and set verified=YES.\n")


def available_drugs(standard="CLSI"):
    """Drugs with VERIFIED breakpoints under a given standard."""
    return sorted(BREAKPOINTS.get(standard.upper(), {}))
