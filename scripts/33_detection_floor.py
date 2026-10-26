#!/usr/bin/env python3
"""33_detection_floor.py — how divergent can a carbapenemase be and still be found?

The allele search reported no uncatalogued carbapenemase. A negative result is
only as strong as the sensitivity behind it, and the positive control in
Table S19a demonstrates sensitivity to alleles that are already in the
catalogue, which is the easy case. This measures the other boundary. It takes
the carbapenemases the cohort actually carries, mutates each one to a ladder of
amino-acid identities, and puts every mutant through the two stages of the
search unchanged.

Stage 1 is the gate that matters. Every predicted protein is aligned against one
representative per beta-lactamase family and discarded unless it clears
E <= 1e-10. A protein that fails there is never seen again, whatever it is.
Stage 2 then finds the closest reference allele and the classifier reads its
identity and coverage.

Three divergence models are run, because one of them alone would be arguable.
"conservative" substitutes within a residue's own physicochemical group four
times out of five, which is what real divergence mostly looks like and is the
easiest case for an aligner. "uniform" substitutes at random over the other
nineteen residues, which is the hardest case and is biologically unrealistic.
"indel" applies the conservative model and then deletes and inserts short
stretches totalling five percent of the length, which is the only model that
also disturbs coverage. Reporting all three brackets the answer rather than
asserting it. Gene calling is not simulated
because it does not depend on similarity to anything: pyrodigal finds an open
reading frame from its own coding model, so a divergent gene is called exactly
as a conserved one is. What identity does change is the alignment, which is
what this measures.

USAGE
    python3 scripts/33_detection_floor.py [--replicates 5] [--seed 0]
"""
import argparse
import csv
import json
import os
import random
import sys

import pyswrd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from allele_classify import classify  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "results", "allele_search")
REF = os.path.join(WORK, "ref", "AMRProt.fa")
OUT = os.path.join(ROOT, "results", "tables", "allele_search_detection_floor.tsv")
SUMMARY = os.path.join(ROOT, "results", "tables",
                       "allele_search_detection_floor_summary.tsv")

# the carbapenemase families the two discordant cohorts actually carry, plus the
# two further families present in the wider collection
FAMILIES = ["blaKPC", "blaNDM", "blaOXA-48_fam", "blaVIM", "blaIMP", "blaGES"]
LADDER = [0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50,
          0.45, 0.40, 0.35, 0.30]
AA = "ACDEFGHIKLMNPQRSTVWY"

# BLOSUM62 off-diagonal rows, used only to weight which residue replaces which.
# A substitution model matters here because replacing every residue uniformly
# at random understates detectability: real divergence is biased towards
# conservative changes, which an alignment still scores.
B62 = {
    "A": "ACDEFGHIKLMNPQRSTVWY", "C": "ACDEFGHIKLMNPQRSTVWY",
    "D": "ACDEFGHIKLMNPQRSTVWY", "E": "ACDEFGHIKLMNPQRSTVWY",
    "F": "ACDEFGHIKLMNPQRSTVWY", "G": "ACDEFGHIKLMNPQRSTVWY",
    "H": "ACDEFGHIKLMNPQRSTVWY", "I": "ACDEFGHIKLMNPQRSTVWY",
    "K": "ACDEFGHIKLMNPQRSTVWY", "L": "ACDEFGHIKLMNPQRSTVWY",
    "M": "ACDEFGHIKLMNPQRSTVWY", "N": "ACDEFGHIKLMNPQRSTVWY",
    "P": "ACDEFGHIKLMNPQRSTVWY", "Q": "ACDEFGHIKLMNPQRSTVWY",
    "R": "ACDEFGHIKLMNPQRSTVWY", "S": "ACDEFGHIKLMNPQRSTVWY",
    "T": "ACDEFGHIKLMNPQRSTVWY", "V": "ACDEFGHIKLMNPQRSTVWY",
    "W": "ACDEFGHIKLMNPQRSTVWY", "Y": "ACDEFGHIKLMNPQRSTVWY",
}
# conservative groups, used as the weighting. A replacement is drawn from the
# residue's own group four times out of five and from the rest otherwise.
GROUPS = ["AGSTP", "ILVM", "FYW", "KRH", "DENQ", "C"]
GROUP_OF = {a: g for g in GROUPS for a in g}


def load_reference(path):
    recs, head, seq = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if head:
                recs.append((head, "".join(seq)))
            head, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if head:
        recs.append((head, "".join(seq)))
    refs = []
    for h, s in recs:
        f = h.split("|")
        if len(f) > 8 and f[8] == "BETA-LACTAM":
            refs.append({"acc": f[0], "symbol": f[3], "family": f[4],
                         "subclass": f[7], "seq": s})
    best = {}
    for i, r in enumerate(refs):
        cur = best.get(r["family"])
        if cur is None or len(r["seq"]) > len(refs[cur]["seq"]):
            best[r["family"]] = i
    return refs, sorted(best.values())


def mutate(seq, target_identity, rng, model="conservative"):
    """Diverge a sequence to a target identity under one of three models."""
    n = len(seq)
    n_mut = int(round((1.0 - target_identity) * n))
    positions = rng.sample(range(n), min(n_mut, n))
    out = list(seq)
    for i in positions:
        a = out[i]
        alt = [x for x in AA if x != a]
        if model == "uniform" or a not in GROUP_OF:
            out[i] = rng.choice(alt)
            continue
        pool = GROUP_OF[a].replace(a, "")
        out[i] = rng.choice(pool) if (pool and rng.random() < 0.8) else rng.choice(alt)
    mutated = "".join(out)
    same = sum(1 for a, b in zip(seq, mutated) if a == b)
    identity = same / n

    if model == "indel":
        # five percent of the length removed and the same amount added, in
        # short stretches, which is what disturbs coverage rather than identity
        budget = max(1, int(round(0.05 * n)))
        chars = list(mutated)
        spent = 0
        while spent < budget and len(chars) > 60:
            k = rng.randint(1, 6)
            i = rng.randrange(0, max(1, len(chars) - k))
            del chars[i:i + k]
            spent += k
        spent = 0
        while spent < budget:
            k = rng.randint(1, 6)
            i = rng.randrange(0, len(chars))
            chars[i:i] = [rng.choice(AA) for _ in range(k)]
            spent += k
        mutated = "".join(chars)
    return mutated, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--models", default="conservative,uniform,indel")
    a = ap.parse_args()

    refs, rep_idx = load_reference(REF)
    rep_seqs = [refs[i]["seq"] for i in rep_idx]
    all_seqs = [r["seq"] for r in refs]

    donors = []
    for fam in FAMILIES:
        members = [r for r in refs
                   if r["family"] == fam and r["subclass"] == "CARBAPENEM"]
        if not members:
            print(f"  no CARBAPENEM reference for {fam}, skipped")
            continue
        # the longest member stands for the family, as in the search itself
        donors.append(max(members, key=lambda r: len(r["seq"])))
    print(f"  {len(donors)} donor enzyme(s): "
          + ", ".join(f'{d["symbol"]} ({d["family"]})' for d in donors))

    models = [m.strip() for m in a.models.split(",") if m.strip()]
    rng = random.Random(a.seed)
    queries, meta = [], []
    for model in models:
        for d in donors:
            for p in LADDER:
                for rep in range(a.replicates):
                    mut, realised = mutate(d["seq"], p, rng, model)
                    queries.append(mut)
                    meta.append({"model": model,
                                 "donor_symbol": d["symbol"],
                                 "donor_family": d["family"],
                                 "donor_acc": d["acc"],
                                 "donor_aa_length": len(d["seq"]),
                                 "mutant_aa_length": len(mut),
                                 "target_identity": p,
                                 "realised_identity_to_donor": round(realised, 4),
                                 "replicate": rep})
    print(f"  {len(queries)} mutant protein(s) to screen")

    # stage 1, exactly as the search runs it
    passed = set()
    for hit in pyswrd.search(queries, rep_seqs, max_alignments=1,
                             max_evalue=1e-10):
        passed.add(hit.query_index)
    for m in meta:
        m["stage1_screened_in"] = False
    for i in passed:
        meta[i]["stage1_screened_in"] = True
    print(f"  stage 1 kept {len(passed)} of {len(queries)}")

    # stage 2, on whatever stage 1 kept
    sel = sorted(passed)
    best = {}
    if sel:
        for hit in pyswrd.search([queries[i] for i in sel], all_seqs,
                                 max_alignments=5, max_evalue=1e-5):
            r = hit.result
            ident = r.identity()
            ref = refs[hit.target_index]
            cov = (r.target_end - r.target_start + 1) / len(ref["seq"])
            cur = best.get(hit.query_index)
            if cur is None or ident > cur["identity"]:
                best[hit.query_index] = {
                    "nearest_symbol": ref["symbol"],
                    "nearest_family": ref["family"],
                    "subclass": ref["subclass"],
                    "identity": round(ident, 4),
                    "ref_coverage": round(cov, 4),
                }
    for j, i in enumerate(sel):
        b = best.get(j)
        if b is None:
            meta[i].update(stage2_hit=False, nearest_symbol="", nearest_family="",
                           subclass="", identity="", ref_coverage="", verdict="")
            continue
        meta[i].update(stage2_hit=True, **b)
        meta[i]["verdict"] = classify(b["identity"], b["ref_coverage"],
                                      b["subclass"])
    for m in meta:
        m.setdefault("stage2_hit", False)
        for k in ("nearest_symbol", "nearest_family", "subclass", "identity",
                  "ref_coverage", "verdict"):
            m.setdefault(k, "")

    cols = ["model", "donor_symbol", "donor_family", "donor_acc",
            "donor_aa_length", "mutant_aa_length",
            "target_identity", "realised_identity_to_donor", "replicate",
            "stage1_screened_in", "stage2_hit", "nearest_symbol",
            "nearest_family", "subclass", "identity", "ref_coverage", "verdict"]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for m in meta:
            w.writerow({k: m.get(k, "") for k in cols})
    print(f"[out] {OUT}")

    # one row per donor and rung, which is what the paper quotes
    rows = []
    for model in models:
      for d in donors:
        for p in LADDER:
            sub = [m for m in meta if m["model"] == model
                   and m["donor_symbol"] == d["symbol"]
                   and m["target_identity"] == p]
            n = len(sub)
            s1 = sum(1 for m in sub if m["stage1_screened_in"])
            carb = sum(1 for m in sub if m.get("subclass") == "CARBAPENEM")
            flagged = sum(1 for m in sub
                          if m.get("verdict") in ("novel_candidate",
                                                  "divergent_homologue",
                                                  "partial_divergent"))
            rows.append({
                "model": model,
                "donor_symbol": d["symbol"], "donor_family": d["family"],
                "target_identity": p, "replicates": n,
                "screened_in_stage1": s1,
                "nearest_is_carbapenemase": carb,
                "reported_as_divergent_or_novel": flagged,
                "detected_pct": round(100 * s1 / n, 1) if n else "",
            })
    with open(SUMMARY, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"[out] {SUMMARY}")

    print("\n  screened in at stage 1, and the verdict the classifier gives,")
    print("  pooled over the six donor enzymes, per divergence model")
    import collections as _c
    for model in models:
        print(f"\n    model: {model}")
        print("      {:<10s}{:>12s}  {}".format("identity", "screened in",
                                                "verdicts"))
        for p in LADDER:
            sub = [m for m in meta if m["model"] == model
                   and m["target_identity"] == p]
            s1 = sum(1 for m in sub if m["stage1_screened_in"])
            vs = _c.Counter(m.get("verdict") or "none" for m in sub)
            top = ", ".join(f"{k} {v}" for k, v in vs.most_common(3))
            print(f"      {p:<10.2f}{s1:>7d}/{len(sub):<4d}  {top}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
