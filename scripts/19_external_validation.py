#!/usr/bin/env python3
"""
19_external_validation.py — test the derived rules on an independent collection.

The rules in results/tables/rules_decision_list*.tsv are fitted on BV-BRC
isolates. This script applies them, unchanged, to NCBI Pathogen Detection
isolates that are not in that set, and reports the predictive values obtained.
Nothing is refitted here. A rule that is re-estimated on the validation data is
not being validated.

DATA

NCBI Pathogen Detection publishes, per isolate, submitted antimicrobial
susceptibility results (AST_phenotypes), AMRFinderPlus genotypes
(AMR_genotypes) including ompK35/ompK36 protein changes, and precomputed
assembly statistics. That combination allows the family and porin tiers to be
evaluated without retyping 165,000 assemblies.

    scripts/19b_fetch_ncbi_pathogen.sh      retrieves and filters the metadata
    scripts/19a_fetch_derivation_accessions.py   accessions used for exclusion

FOUR DIFFERENCES FROM THE DERIVATION SET, ALL OF WHICH MATTER

1. Genotyping tool. Carbapenemase and porin calls come from AMRFinderPlus, not
   Kleborate. The carbapenemase family dictionary is reused from kp_defs.py so
   family assignment is identical, but porin classification must be mapped
   between two different notations, and the mapping is approximate. In
   particular Kleborate's ompK36 c.25C>T reduced-expression call is a
   nucleotide-level variant that AMRFinderPlus does not report, so that
   category cannot be evaluated externally and is reported as unavailable
   rather than merged into another class.

2. Susceptibility testing. The derivation set was restricted to CLSI broth
   dilution. NCBI AST submissions do not carry a uniform method or standard, so
   this set is method-heterogeneous. An earlier version of this note claimed
   that heterogeneity biases towards the null. It does not follow: a different
   method or standard can place a borderline isolate on either side of the
   breakpoint, and whether that helps or hurts a given rule depends on where
   its errors already sit. The direction cannot be determined from these data,
   so external and derivation estimates are not directly comparable in either
   direction.

3. Only interpreted calls are available, not MIC values, so the
   censoring-aware classification used in the derivation cannot be reproduced.
   Calls of ND and DD are dropped; only S and R enter the binary set.

4. Incomplete gene calls. AMRFinderPlus flags calls as PARTIAL,
   PARTIAL_END_OF_CONTIG or MISTRANSLATION. By default a flagged call does not
   establish gene presence, because a partial hit at a contig boundary is
   precisely the fragmentation artefact under investigation. --count-partial
   reruns the analysis counting them, and both results are reported: the
   difference is an estimate of how much the genotype depends on assembly
   contiguity.

USAGE

    python scripts/19_external_validation.py
    # NCBI publishes 3 of the 7 QualiBact metrics, so --qc-filter drops only
    # the isolates failing those 3. It is not a QualiBact PASS filter and no
    # external isolate holds a full QualiBact PASS.
    python scripts/19_external_validation.py --qc-filter
    python scripts/19_external_validation.py --count-partial
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from kp_defs import OXA48_LIKE, allele_family  # noqa: E402

TAB = Path("results/tables")
SUP = Path("results/supplementary")

# Acquired carbapenemases outside the OXA-48-like group. This list exists
# because the OXA family is overwhelmingly NOT carbapenemase: blaOXA-1,
# blaOXA-2, blaOXA-9 and blaOXA-10 together appear tens of thousands of times
# in this collection and are narrow- or extended-spectrum penicillinases.
# Counting any blaOXA as a carbapenemase would invent a carbapenemase for a
# large fraction of the validation set.
CARBAPENEMASE_OXA_OTHER = {
    "OXA-23", "OXA-24", "OXA-40", "OXA-25", "OXA-26", "OXA-27", "OXA-49",
    "OXA-58", "OXA-72", "OXA-96", "OXA-143", "OXA-182", "OXA-235", "OXA-236",
    "OXA-237", "OXA-253", "OXA-420",
}
CARB_PREFIX = re.compile(r"^bla(KPC|NDM|VIM|IMP|GES|SPM)[-_]", re.I)
BLA_OXA = re.compile(r"^blaOXA[-_]", re.I)

UNTRUSTED_FLAGS = {"PARTIAL", "PARTIAL_END_OF_CONTIG", "MISTRANSLATION"}

# AMRFinderPlus protein-change notation for the porins.
_OMPK35 = re.compile(r"^ompK35_", re.I)
_OMPK36 = re.compile(r"^ompK36_", re.I)
_TRUNC = re.compile(r"(fsTer|Ter\d*$|fs$)", re.I)
# Loop-3 di-residue insertion: a residue at 134-137 replaced by itself plus
# GD / TD / DT, e.g. D135DGD, T136TDT. This is the classic OmpK36 loop-3 event.
_LOOP3 = re.compile(r"^ompK36_([A-Z])(13[4-7])\1(GD|TD|DT|DG)$", re.I)
_SINGLE_INS = re.compile(r"^ompK36_([A-Z])(13[4-7])\1([A-Z])$", re.I)

# QualiBact v1.0 Klebsiella pneumoniae FAIL bounds for the three metrics NCBI
# precomputes. GC, completeness, contamination and CDS count are unavailable
# here, so this is a partial application of the scheme and is labelled as such.
QC_FAIL = {"genome_size_bp": (4_900_000, 6_400_000),
           "n50_bp": (58_000, None),
           "n_contigs": (None, 350)}


def wilson(k: int, n: int, z: float = 1.959963984540054):
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = (z / d) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (round(100 * p, 1), round(100 * max(0.0, c - h), 1),
            round(100 * min(1.0, c + h), 1))



def bioproject_bootstrap(frame, preds_by_tier, tiers, n_boot=2000, seed=20261013):
    """Uncertainty for the gain from the porin rule, resampling BioProjects.

    The external cohort is a stack of submissions. Isolates within one
    submission share a laboratory, a testing method and often a clone, so an
    interval that treats them as independent is too narrow, and the quantity
    that matters is paired: the same isolates scored by two rules. This
    resamples whole BioProjects with replacement and reports percentile
    intervals for the difference each rule makes, which is the comparison the
    text draws.
    """
    a_tier, b_tier = tiers
    y = frame.y.to_numpy()
    pa = np.asarray(preds_by_tier[a_tier])
    pb = np.asarray(preds_by_tier[b_tier])
    groups = frame.bioproject.fillna("<none>").to_numpy()
    idx_by_group = {g: np.flatnonzero(groups == g) for g in np.unique(groups)}
    keys = list(idx_by_group)
    rng = np.random.default_rng(seed)

    def metrics(sel):
        yy, aa, bb = y[sel], pa[sel], pb[sel]
        out = {}
        for name, pp in (("a", aa), ("b", bb)):
            tp = int(((pp == 1) & (yy == 1)).sum()); fp = int(((pp == 1) & (yy == 0)).sum())
            tn = int(((pp == 0) & (yy == 0)).sum())
            out[f"{name}_acc"] = 100.0 * (tp + tn) / len(yy) if len(yy) else float("nan")
            out[f"{name}_ppv"] = 100.0 * tp / (tp + fp) if tp + fp else float("nan")
        return out

    obs = metrics(np.arange(len(y)))
    diffs = {"accuracy": [], "PPV": []}
    for _ in range(n_boot):
        draw = rng.choice(len(keys), size=len(keys), replace=True)
        sel = np.concatenate([idx_by_group[keys[i]] for i in draw])
        m = metrics(sel)
        if np.isnan(m["a_ppv"]) or np.isnan(m["b_ppv"]):
            continue
        diffs["accuracy"].append(m["b_acc"] - m["a_acc"])
        diffs["PPV"].append(m["b_ppv"] - m["a_ppv"])
    rows = []
    for metric, key in (("accuracy", "acc"), ("PPV", "ppv")):
        v = np.asarray(diffs[metric])
        rows.append(dict(
            metric=metric, rule_a=a_tier, rule_b=b_tier,
            value_a=round(obs[f"a_{key}"], 1), value_b=round(obs[f"b_{key}"], 1),
            difference=round(obs[f"b_{key}"] - obs[f"a_{key}"], 1),
            ci_lo=round(float(np.percentile(v, 2.5)), 1),
            ci_hi=round(float(np.percentile(v, 97.5)), 1),
            n_bioprojects=len(keys), n_isolates=int(len(y)),
            n_resamples=int(len(v)), method="BioProject cluster bootstrap"))
    return rows

def predictive_values(y_true, y_pred) -> dict:
    tp = sum(1 for t, p in zip(y_true, y_pred) if p == 1 and t == 1)
    fp = sum(1 for t, p in zip(y_true, y_pred) if p == 1 and t == 0)
    tn = sum(1 for t, p in zip(y_true, y_pred) if p == 0 and t == 0)
    fn = sum(1 for t, p in zip(y_true, y_pred) if p == 0 and t == 1)
    out = dict(n=len(y_true), TP=tp, FP=fp, TN=tn, FN=fn)
    for name, k, d in (("PPV", tp, tp + fp), ("NPV", tn, tn + fn),
                       ("sensitivity", tp, tp + fn),
                       ("specificity", tn, tn + fp),
                       ("accuracy", tp + tn, len(y_true))):
        v, lo, hi = wilson(k, d)
        out[name], out[f"{name}_lo"], out[f"{name}_hi"] = v, lo, hi
    return out


# ------------------------------------------------------------------ parsing
def parse_ast(field: str, drug: str) -> str | None:
    for kv in field.split(","):
        if "=" not in kv:
            continue
        d, call = kv.split("=", 1)
        if d.strip().lower() == drug:
            c = call.strip().upper()
            return c if c in {"S", "R", "I"} else None
    return None


def split_tok(tok: str) -> tuple[str, str]:
    base, _, flag = tok.strip().partition("=")
    return base.strip(), flag.strip().upper()


def carb_alleles(field: str, count_partial: bool) -> tuple[list[str], int]:
    """Acquired-carbapenemase alleles. Returns (alleles, n_flagged_skipped)."""
    out, skipped = [], 0
    for tok in field.split(","):
        base, flag = split_tok(tok)
        if not base:
            continue
        is_carb = False
        if CARB_PREFIX.match(base):
            is_carb = True
        elif BLA_OXA.match(base):
            allele = base[3:].upper().replace("_", "-")
            is_carb = (allele in {a.upper() for a in OXA48_LIKE}
                       or allele in {a.upper() for a in CARBAPENEMASE_OXA_OTHER})
        if not is_carb:
            continue
        if flag in UNTRUSTED_FLAGS and not count_partial:
            skipped += 1
            continue
        out.append(base[3:].replace("_", "-"))
    return out, skipped


def porin_class(field: str) -> str:
    """Map AMRFinderPlus ompK changes onto the study's porin categories."""
    k36_trunc = k36_loop3 = k36_other = k35_loss = False
    for tok in field.split(","):
        base, flag = split_tok(tok)
        if not base:
            continue
        if _OMPK35.match(base):
            if _TRUNC.search(base):
                k35_loss = True
        elif _OMPK36.match(base):
            if _TRUNC.search(base):
                k36_trunc = True
            elif _LOOP3.match(base):
                k36_loop3 = True
            elif _SINGLE_INS.match(base):
                k36_other = True
            else:
                k36_other = True
    if k36_trunc:
        return "ompK36_loss"
    if k36_loop3:
        return "ompK36_loop3"
    if k36_other:
        return "ompK36_other_variant"
    if k35_loss:
        return "ompK35_loss"
    return "none"


def porin_level(v: str) -> str:
    if v in {"ompK36_loop3", "ompK36_loss", "ompK36_reduced_expression"}:
        return "ompK36_disrupted"
    if v == "ompK35_loss":
        return "ompK35_only"
    if v == "none":
        return "intact"
    return "unknown"


def family_level(alleles: list[str]) -> str:
    if not alleles:
        return "none"
    fams = sorted({allele_family(a) for a in alleles})
    if len(fams) > 1:
        return "multiple"
    f = fams[0]
    return f if f in {"KPC", "OXA-48-like", "NDM", "VIM"} else "other"


def qc_verdict(row: dict) -> str:
    """Partial QualiBact verdict from NCBI's precomputed assembly statistics."""
    for metric, (lo, hi) in QC_FAIL.items():
        try:
            v = float(row[metric])
        except (TypeError, ValueError, KeyError):
            return "unknown"
        if lo is not None and v < lo:
            return "FAIL"
        if hi is not None and v > hi:
            return "FAIL"
    return "PASS_partial"


# --------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ncbi", type=Path,
                    default=Path("data/external/ncbi_kp_ast.tsv"))
    ap.add_argument("--derivation-accessions", type=Path,
                    default=SUP / "S1_derivation_accessions.tsv")
    ap.add_argument("--rules", type=Path,
                    default=TAB / "rules_decision_list.tsv")
    ap.add_argument("--drug", default="meropenem",
                    choices=["meropenem", "imipenem"])
    ap.add_argument("--count-partial", action="store_true",
                    help="count PARTIAL/MISTRANSLATION gene calls as present")
    ap.add_argument("--qc-filter", action="store_true",
                    help="keep only assemblies passing the partial QualiBact check")
    ap.add_argument("--min-support", type=int, default=10)
    ap.add_argument("--coherence-pct", type=float, default=20.0,
                    help="flag a BioProject if fewer than this %% of its "
                         "carbapenemase-positive isolates are reported resistant")
    ap.add_argument("--allow-no-overlap", action="store_true",
                    help="proceed even if the derivation-accession file is "
                         "absent; the result is then NOT a valid external "
                         "validation and the run is for inspection only")
    ap.add_argument("--tag", default="")
    ap.add_argument("--outdir", type=Path, default=TAB)
    a = ap.parse_args()

    for p in (a.ncbi, a.rules):
        if not p.exists():
            sys.exit(f"Missing {p}")

    # ---- derivation accessions, for exclusion
    excl_asm: set[str] = set()
    excl_bio: set[str] = set()
    if a.derivation_accessions.exists():
        acc = pd.read_csv(a.derivation_accessions, sep="\t", dtype=str)
        for c, s in (("assembly_accession", excl_asm),
                     ("biosample_accession", excl_bio)):
            if c in acc.columns:
                s.update(acc[c].dropna().astype(str))
        print(f"derivation accessions : {len(excl_asm):,} assemblies, "
              f"{len(excl_bio):,} biosamples")
    elif a.allow_no_overlap:
        print("WARNING: no derivation accession file. Overlap cannot be "
              "excluded and these results are NOT a valid external "
              "validation (proceeding because --allow-no-overlap was set).",
              file=sys.stderr)
    else:
        sys.exit(
            f"Missing derivation-accession file {a.derivation_accessions}.\n"
            "Overlap with the derivation set cannot be excluded, so the run "
            "would not be a valid external validation. Generate it with "
            "scripts/19a_fetch_derivation_accessions.py, or pass "
            "--allow-no-overlap to run anyway for inspection only.")

    # ---- read and reduce the NCBI table
    rows, n_total = [], 0
    n_no_ast = n_overlap = n_not_sr = 0
    n_flagged_skipped = 0
    with a.ncbi.open() as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            n_total += 1
            call = parse_ast(rec.get("AST_phenotypes", ""), a.drug)
            if call is None:
                n_no_ast += 1
                continue
            asm = (rec.get("asm_acc") or "").strip()
            bio = (rec.get("biosample_acc") or "").strip()
            if (asm and asm in excl_asm) or (bio and bio in excl_bio):
                n_overlap += 1
                continue
            if call not in {"S", "R"}:
                n_not_sr += 1
                continue
            geno = rec.get("AMR_genotypes", "") or ""
            alleles, skipped = carb_alleles(geno, a.count_partial)
            n_flagged_skipped += skipped
            pc = porin_class(geno)
            rows.append(dict(
                label=rec.get("#label", ""), asm_acc=asm, biosample_acc=bio,
                call=call, y=1 if call == "R" else 0,
                carb_pos=1 if alleles else 0,
                alleles=";".join(sorted(set(alleles))) or "none",
                family=family_level(alleles),
                porin_raw=pc, porin=porin_level(pc),
                genome_size_bp=rec.get("asm_stats_length_bp", ""),
                n50_bp=rec.get("asm_stats_contig_n50", ""),
                n_contigs=rec.get("asm_stats_n_contig", ""),
                collection_date=rec.get("collection_date", ""),
                geo_loc_name=rec.get("geo_loc_name", ""),
                bioproject=rec.get("bioproject_acc", "") or "unknown",
                sequenced_by=rec.get("sequenced_by", "") or "unknown",
                epi_type=rec.get("epi_type", "") or "unknown"))

    d = pd.DataFrame(rows)
    # One BioSample can be represented by more than one assembly in NCBI
    # Pathogen Detection. Key on BioSample and keep the first assembly by
    # accession, so each sampled isolate is counted once.
    n_before_dedup = len(d)
    if not d.empty and "biosample_acc" in d.columns:
        has_bio = d.biosample_acc.astype(str).str.strip() != ""
        dd = d[has_bio].sort_values("asm_acc").drop_duplicates("biosample_acc",
                                                               keep="first")
        d = pd.concat([dd, d[~has_bio]], ignore_index=True)
    n_dup_removed = n_before_dedup - len(d)
    print(f"NCBI isolates scanned : {n_total:,}")
    print(f"  no {a.drug} AST      : {n_no_ast:,}")
    print(f"  in derivation set    : {n_overlap:,}  (excluded)")
    print(f"  call not S or R      : {n_not_sr:,}  (excluded)")
    print(f"  duplicate BioSample  : {n_dup_removed:,}  (collapsed, first "
          f"assembly by accession kept)")
    if not a.count_partial:
        print(f"  flagged gene calls not counted as present: {n_flagged_skipped:,}")
    if d.empty:
        sys.exit("No external isolates remain. Nothing to validate.")

    # Archive the exclusion ledger so the external cohort is reconstructible.
    a.outdir.mkdir(parents=True, exist_ok=True)
    with (a.outdir / f"external_exclusions{a.tag}.tsv").open("w") as fh:
        fh.write("stage\tn\n")
        fh.write(f"ncbi_isolates_scanned\t{n_total}\n")
        fh.write(f"no_{a.drug}_ast\t{n_no_ast}\n")
        fh.write(f"overlap_with_derivation_excluded\t{n_overlap}\n")
        fh.write(f"call_not_S_or_R_excluded\t{n_not_sr}\n")
        fh.write(f"duplicate_biosample_collapsed\t{n_dup_removed}\n")
        fh.write(f"external_cohort_final\t{len(d)}\n")

    d["qc"] = [qc_verdict(r) for r in d.to_dict("records")]
    print(f"\nexternal set          : {len(d):,} isolates "
          f"({int(d.y.sum()):,} resistant, {len(d) - int(d.y.sum()):,} susceptible)")
    print(f"  partial QualiBact    : "
          + ", ".join(f"{k} {v:,}" for k, v in d.qc.value_counts().items()))
    if a.qc_filter:
        before = len(d)
        d = d[d.qc == "PASS_partial"].copy()
        print(f"  --qc-filter applied  : {before:,} -> {len(d):,}")
        if d.empty:
            sys.exit("No isolates pass QC. Nothing to validate.")

    a.outdir.mkdir(parents=True, exist_ok=True)
    Path("data/external").mkdir(parents=True, exist_ok=True)
    d.to_csv(a.outdir / f"external_isolates{a.tag}.tsv", sep="\t", index=False)

    # ------------------------------------------- submission-batch screening
    # A public collection is an aggregate of submissions, and a submission can
    # carry systematically inconsistent AST metadata. If one batch is large
    # relative to the validation set, it can determine the result on its own.
    # Each BioProject is therefore screened on an objective internal-coherence
    # criterion before any pooled figure is reported.
    proj_rows = []
    for bp, g in d.groupby("bioproject"):
        gp = g[g.carb_pos == 1]
        n_pos = len(gp)
        r_pos = int(gp.y.sum()) if n_pos else 0
        pct_pos = round(100 * r_pos / n_pos, 1) if n_pos else float("nan")
        # A carbapenemase-positive isolate reported susceptible is the
        # phenomenon under study, so a low rate is not by itself suspicious.
        # A batch in which almost every carbapenemase-positive isolate is
        # susceptible, however, disagrees with every published series and is
        # far more likely to be a metadata error than a discovery.
        incoherent = bool(n_pos >= a.min_support and pct_pos < a.coherence_pct)
        proj_rows.append(dict(
            bioproject=bp, sequenced_by=g.sequenced_by.mode().iat[0],
            n=len(g), n_resistant=int(g.y.sum()),
            resistance_pct=round(100 * g.y.mean(), 1),
            n_carbapenemase_positive=n_pos,
            carb_pos_resistant_pct=pct_pos,
            pct_of_external_set=round(100 * len(g) / len(d), 1),
            flagged_incoherent="yes" if incoherent else "no"))
    proj = pd.DataFrame(proj_rows).sort_values("n", ascending=False)
    proj.to_csv(a.outdir / f"external_bioprojects{a.tag}.tsv",
                sep="\t", index=False)

    flagged = set(proj.loc[proj.flagged_incoherent == "yes", "bioproject"])
    print("\n" + "=" * 90)
    print("EXTERNAL SET COMPOSITION BY SUBMISSION")
    print("=" * 90)
    print(f"{'bioproject':16s}{'sequenced_by':30s}{'n':>5s}{'%set':>6s}"
          f"{'R%':>7s}{'carb+':>7s}{'carb+ R%':>10s}  flag")
    for _, r in proj.head(12).iterrows():
        cp = ("  n/a" if pd.isna(r.carb_pos_resistant_pct)
              else f"{r.carb_pos_resistant_pct:5.1f}")
        print(f"{r.bioproject:16s}{str(r.sequenced_by)[:29]:30s}{r.n:5d}"
              f"{r.pct_of_external_set:6.1f}{r.resistance_pct:7.1f}"
              f"{r.n_carbapenemase_positive:7d}{cp:>10s}"
              f"  {'FLAG' if r.flagged_incoherent == 'yes' else ''}")
    if flagged:
        nf = int(proj.loc[proj.flagged_incoherent == "yes", "n"].sum())
        print(f"\n{len(flagged)} BioProject(s) flagged: among their "
              f"carbapenemase-positive isolates fewer than\n{a.coherence_pct:.0f}% "
              f"are reported resistant. Together they are {nf:,} isolates "
              f"({100 * nf / len(d):.1f}% of\nthe external set). Results are "
              f"reported with and without them; neither is\nsuppressed.")
    else:
        print("\nNo BioProject met the incoherence criterion.")

    # ---- load the rules fitted on the derivation set
    rules = pd.read_csv(a.rules, sep="\t")
    lut: dict[str, dict] = {}
    for _, r in rules.iterrows():
        key = tuple(part.split("=", 1)[1] for part in str(r.rule).split(" AND "))
        lut.setdefault(r.tier, {})[key] = 1 if r.then_predict == "RESISTANT" else 0

    TIER_KEYS = {"T0_gene_only": ["carb_pos"],
                 "T1_family": ["family"],
                 "T2_family_porin": ["family", "porin"]}
    ORDER = ["T2_family_porin", "T1_family", "T0_gene_only"]

    def apply_rules(frame: pd.DataFrame, tier: str) -> tuple[list[int], int]:
        chain = ORDER[ORDER.index(tier):]
        preds, unseen = [], 0
        for rec in frame.to_dict("records"):
            p = None
            for t in chain:
                key = tuple(str(rec[k]) for k in TIER_KEYS[t])
                if t in lut and key in lut[t]:
                    p = lut[t][key]
                    break
            if p is None:
                unseen += 1
                p = 0
            preds.append(p)
        return preds, unseen

    def perf_for(frame: pd.DataFrame, label: str) -> list[dict]:
        rows = []
        for tier in TIER_KEYS:
            preds, unseen = apply_rules(frame, tier)
            rows.append(dict(subset=label, tier=tier, drug=a.drug,
                             unseen_strata=unseen,
                             count_partial=a.count_partial,
                             qc_filtered=a.qc_filter,
                             **predictive_values(frame.y.tolist(), preds)))
        return rows

    perf_rows, strat_rows = [], []
    perf_rows += perf_for(d, "all_external")
    d_clean = d[~d.bioproject.isin(flagged)].copy() if flagged else d
    if flagged:
        perf_rows += perf_for(d_clean, "excluding_flagged_bioprojects")

    # The derivation set passed a full QualiBact filter; this one cannot,
    # because NCBI Pathogen Detection publishes only three of the seven
    # metrics. Isolates failing on those three are labelled FAIL and the rest
    # PASS_partial, which is not a QualiBact PASS. Dropping the failures says
    # how much of the external result rests on them.
    if "qc" in d.columns and (d.qc == "FAIL").any():
        perf_rows += perf_for(d[d.qc != "FAIL"].copy(),
                              "excluding_partial_qc_failures")

    # How much the porin rule adds on this cohort, with an interval that
    # respects the fact that the cohort is a stack of submissions.
    boot_rows = []
    for label, frame in (("all_external", d),
                         ("excluding_flagged_bioprojects", d_clean)):
        if frame is d_clean and not flagged:
            continue
        preds = {t: apply_rules(frame, t)[0] for t in
                 ("T0_gene_only", "T2_family_porin")}
        for r in bioproject_bootstrap(frame, preds,
                                      ("T0_gene_only", "T2_family_porin")):
            r["subset"] = label
            boot_rows.append(r)
    if boot_rows:
        pd.DataFrame(boot_rows).to_csv(
            a.outdir / f"external_porin_gain_bootstrap{a.tag}.tsv",
            sep="\t", index=False)
        print("\nGain from the family-specific porin rule, BioProject cluster "
              "bootstrap:")
        for r in boot_rows:
            print(f"  {r['subset']:<32} {r['metric']:<9} "
                  f"{r['value_a']:5.1f} -> {r['value_b']:5.1f}  "
                  f"difference {r['difference']:+5.1f} "
                  f"(95% CI {r['ci_lo']:+.1f} to {r['ci_hi']:+.1f}, "
                  f"{r['n_bioprojects']} BioProjects)")

    # ---- leave-one-BioProject-out, so no single submission drives the result
    loo_rows = []
    big = proj[proj.n >= a.min_support].bioproject.tolist()
    for bp in big:
        sub = d_clean[d_clean.bioproject != bp]
        if len(sub) < 50 or sub.y.nunique() < 2:
            continue
        for r in perf_for(sub, f"leave_out_{bp}"):
            r["left_out"] = bp
            # Counted against the full eligible cohort, not against d_clean.
            # A flagged project is already out of d_clean, so counting there
            # reported 0 isolates removed for the very project whose removal
            # the row is about.
            r["left_out_n"] = int((d.bioproject == bp).sum())
            r["left_out_n_in_source"] = int((d_clean.bioproject == bp).sum())
            r["source_cohort"] = ("excluding_flagged_bioprojects" if flagged
                                  else "all_external")
            loo_rows.append(r)
    loo = pd.DataFrame(loo_rows)
    if not loo.empty:
        loo.to_csv(a.outdir / f"external_leave_one_project_out{a.tag}.tsv",
                   sep="\t", index=False)

    # ---- stratum rates, external vs derivation
    # Computed on the screened set: comparing a stratum rate against a
    # derivation rate is only meaningful if the external rate is not itself
    # dominated by a batch with incoherent metadata.
    for key, g in d_clean.groupby(["family", "porin"]):
        n, nr = len(g), int(g.y.sum())
        pct, lo, hi = wilson(nr, n)
        strat_rows.append(dict(
            stratum=f"family={key[0]} | porin={key[1]}",
            family=key[0], porin=key[1], n=n, n_resistant=nr,
            resistance_pct=pct, ci_lower=lo, ci_upper=hi,
            low_support="yes" if n < a.min_support else "no"))
    ext = pd.DataFrame(strat_rows)

    der_path = a.outdir / "rules_strata.tsv"
    if der_path.exists():
        der = pd.read_csv(der_path, sep="\t")
        der = der[der.tier == "T2_family_porin"][
            ["stratum", "n", "n_resistant", "resistance_pct",
             "ci_lower", "ci_upper"]]
        cmp = der.merge(ext, on="stratum", how="outer",
                        suffixes=("_derivation", "_external"))
        cmp["delta_pct"] = (cmp.resistance_pct_external
                            - cmp.resistance_pct_derivation).round(1)
        cmp["ci_overlap"] = ~((cmp.ci_lower_external > cmp.ci_upper_derivation)
                              | (cmp.ci_upper_external < cmp.ci_lower_derivation))
        cmp.to_csv(a.outdir / f"external_stratum_comparison{a.tag}.tsv",
                   sep="\t", index=False)

        print("\n" + "=" * 90)
        print("STRATUM RESISTANCE RATE: DERIVATION vs EXTERNAL")
        print("=" * 90)
        print(f"{'stratum':44s}{'derivation':>22s}{'external':>22s}  agree")
        both = cmp[(cmp.n_derivation.notna()) & (cmp.n_external.notna())
                   & (cmp.n_external >= a.min_support)
                   & (cmp.n_derivation >= a.min_support)]
        for _, r in both.sort_values("resistance_pct_derivation",
                                     ascending=False).iterrows():
            print(f"{r.stratum:44s}"
                  f"{r.resistance_pct_derivation:7.1f} "
                  f"[{r.ci_lower_derivation:5.1f}-{r.ci_upper_derivation:5.1f}] "
                  f"n={int(r.n_derivation):<5d}"
                  f"{r.resistance_pct_external:7.1f} "
                  f"[{r.ci_lower_external:5.1f}-{r.ci_upper_external:5.1f}] "
                  f"n={int(r.n_external):<5d}"
                  f"  {'yes' if r.ci_overlap else 'NO'}")
        n_ok = int(both.ci_overlap.sum())
        print(f"\nconfidence intervals overlap in {n_ok} of {len(both)} strata "
              f"with n >= {a.min_support} in both sets.")

    perf = pd.DataFrame(perf_rows + loo_rows)
    perf.to_csv(a.outdir / f"external_performance{a.tag}.tsv",
                sep="\t", index=False)
    ext.to_csv(a.outdir / f"external_strata{a.tag}.tsv", sep="\t", index=False)

    print("\n" + "=" * 90)
    print(f"EXTERNAL PREDICTIVE VALUES  —  {a.drug}, rules fitted on BV-BRC, "
          f"applied unchanged")
    print("=" * 90)
    for label in perf.subset.unique():
        if label.startswith("leave_out_"):
            continue
        sub = perf[perf.subset == label]
        print(f"\n{label}   (n = {int(sub.iloc[0].n):,})")
        for _, r in sub.iterrows():
            print(f"  {r.tier:18s} PPV {r.PPV:5.1f} "
                  f"[{r.PPV_lo:4.1f}-{r.PPV_hi:4.1f}]   "
                  f"NPV {r.NPV:5.1f} [{r.NPV_lo:4.1f}-{r.NPV_hi:4.1f}]   "
                  f"sens {r.sensitivity:5.1f}  spec {r.specificity:5.1f}  "
                  f"acc {r.accuracy:5.1f}")

    if not loo.empty:
        print("\n" + "=" * 90)
        print("LEAVE-ONE-BIOPROJECT-OUT  —  range of each statistic across "
              "dropped submissions")
        print("=" * 90)
        print(f"{'tier':20s}{'PPV min-max':>18s}{'NPV min-max':>18s}"
              f"{'accuracy min-max':>20s}   projects")
        for tier in TIER_KEYS:
            s = loo[loo.tier == tier]
            if s.empty:
                continue
            print(f"{tier:20s}"
                  f"{s.PPV.min():8.1f}-{s.PPV.max():<9.1f}"
                  f"{s.NPV.min():8.1f}-{s.NPV.max():<9.1f}"
                  f"{s.accuracy.min():10.1f}-{s.accuracy.max():<9.1f}"
                  f"   {s.left_out.nunique()}")
        worst = loo.loc[loo.groupby("tier").accuracy.idxmin()]
        print("\nsubmission whose removal lowers accuracy most (i.e. the one "
              "propping it up):")
        for _, r in worst.iterrows():
            print(f"  {r.tier:20s} {r.left_out:16s} "
                  f"(n = {int(r.left_out_n):,})  accuracy without it "
                  f"{r.accuracy:.1f}")

    print(f"\n-> {a.outdir}/external_performance{a.tag}.tsv")
    print(f"-> {a.outdir}/external_bioprojects{a.tag}.tsv")
    print("\nNOTE: this set is method-heterogeneous and genotyped with "
          "AMRFinderPlus rather\nthan Kleborate, its susceptibility testing is "
          "not standardised, and only three\nof the seven QualiBact metrics "
          "are published for it. These differences limit\ndirect comparison "
          "with the derivation estimates. Their net effect on measured\n"
          "performance could run in either direction, so this is not a bound "
          "on\ntransferability in either sense.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
