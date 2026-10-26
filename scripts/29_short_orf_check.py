#!/usr/bin/env python3
"""29_short_orf_check.py — the band the main search left out.

28_allele_search.py screens predicted proteins between 180 and 500 amino acids,
because an intact beta-lactamase is about 250 to 350. That window cannot see a
badly truncated gene, which is exactly what a disrupted carbapenemase would
look like. This covers the missing band, 60 to 180 amino acids, against the
carbapenemase proteins only, so the earlier null result is not an artefact of
the length filter.

USAGE
    python3 29_short_orf_check.py COHORT [--seconds 150]
"""
import json
import os
import sys
import time

import pyrodigal
import pyswrd

# the scripts live in scripts/ and their data in results/allele_search/, so
# both are resolved from the repository root rather than from the file's own
# directory, which lets the script be moved without silently writing elsewhere
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "results", "allele_search")
OUT = os.path.join(WORK, "short_orfs.jsonl")
LO, HI = 60, 180


def load_carbapenemase_refs(path):
    recs, h, s = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if h:
                recs.append((h, "".join(s)))
            h, s = line[1:].strip(), []
        else:
            s.append(line.strip())
    if h:
        recs.append((h, "".join(s)))
    refs = []
    for head, seq in recs:
        f = head.split("|")
        if len(f) > 8 and f[7] == "CARBAPENEM":
            refs.append({"symbol": f[3], "family": f[4], "seq": seq})
    reps = {}
    for i, r in enumerate(refs):
        cur = reps.get(r["family"])
        if cur is None or len(r["seq"]) > len(refs[cur]["seq"]):
            reps[r["family"]] = i
    return refs, sorted(reps.values())


def main():
    cohort = sys.argv[1]
    secs = int(sys.argv[sys.argv.index("--seconds") + 1]) if "--seconds" in sys.argv else 150
    refs, rep_idx = load_carbapenemase_refs(os.path.join(WORK, "ref", "AMRProt.fa"))
    reps = [refs[i]["seq"] for i in rep_idx]
    allr = [r["seq"] for r in refs]
    ids = json.load(open(os.path.join(WORK, "cohorts.json")))[cohort]
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            try:
                r = json.loads(line)
                if r.get("cohort") == cohort:
                    done.add(r["genome_id"])
            except Exception:
                pass
    todo = [g for g in ids if g not in done]
    print(f"  {cohort}: {len(done)} done, {len(todo)} remaining, "
          f"{len(refs)} carbapenemase refs in {len(reps)} families")
    gf = pyrodigal.GeneFinder(meta=True)
    fh = open(OUT, "a")
    t0 = time.time()
    n = 0
    for gid in todo:
        if time.time() - t0 > secs:
            break
        path = None
        for d in ("data/raw/genomes_redownload", "data/raw/genomes"):
            p = os.path.join(ROOT, d, gid + ".fna")
            if os.path.exists(p):
                path = p
                break
        contigs, seq = [], []
        for line in open(path):
            if line.startswith(">"):
                if seq:
                    contigs.append("".join(seq))
                    seq = []
            else:
                seq.append(line.strip())
        if seq:
            contigs.append("".join(seq))
        prots = []
        for c in contigs:
            if len(c) < 100:
                continue
            for g in gf.find_genes(c):
                core = g.translate().rstrip("*")
                if LO <= len(core) <= HI:
                    prots.append(core)
        keep = set()
        if prots:
            for hit in pyswrd.search(prots, reps, max_alignments=1, max_evalue=1e-5):
                keep.add(hit.query_index)
        found = []
        if keep:
            sel = [prots[i] for i in sorted(keep)]
            best = {}
            for hit in pyswrd.search(sel, allr, max_alignments=3, max_evalue=1e-3):
                r = hit.result
                ident = r.identity()
                ref = refs[hit.target_index]
                cur = best.get(hit.query_index)
                if cur is None or ident > cur["identity"]:
                    best[hit.query_index] = {
                        "aa_length": len(sel[hit.query_index]),
                        "ref_symbol": ref["symbol"],
                        "ref_aa_length": len(ref["seq"]),
                        "identity": round(ident, 4),
                        "ref_coverage": round(
                            (r.target_end - r.target_start + 1) / len(ref["seq"]), 4),
                    }
            found = list(best.values())
        fh.write(json.dumps({"cohort": cohort, "genome_id": gid,
                             "short_orfs_screened": len(prots),
                             "carbapenemase_like": found}) + "\n")
        fh.flush()
        n += 1
    fh.close()
    print(f"  processed {n} genome(s) in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
