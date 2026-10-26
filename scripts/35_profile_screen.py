#!/usr/bin/env python3
"""35_profile_screen.py — screen the same proteins with family profiles.

Complements 28_allele_search.py. That script asks whether a predicted protein
matches a catalogued allele closely enough to be called one. This asks the
prior question: is the protein a beta-lactamase at all, judged by a profile
built from the whole family rather than by similarity to any single member.

A protein can fail the first question and pass this one. That is exactly the
case the editor's question is about, a beta-lactamase too divergent from every
catalogued allele to be named but still recognisable as one of the family.

The control cohort is screened on the same terms, because every Klebsiella
genome carries chromosomal proteins of the beta-lactamase fold and the only
meaningful reading of a profile hit is whether the discordant cohorts carry
more of them, or different ones, than the controls.

Protein names reproduce 28_allele_search.py exactly, so a hit here can be
joined to a hit there. The length window spans both windows that script and
29_short_orf_check.py used, 60 to 500 residues.

The run is resumable and time-boxed in the same way.

USAGE
    python3 scripts/35_profile_screen.py COHORT [--seconds 150] [--limit N]
"""
import argparse
import json
import os
import sys
import time

import pyhmmer
import pyrodigal

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "results", "allele_search")
HMMS = os.path.join(WORK, "ref", "betalactamase_families.hmm")
FAMMETA = os.path.join(WORK, "ref", "betalactamase_families.tsv")
OUT = os.path.join(WORK, "profile_hits.jsonl")
GENOME_DIRS = ["data/raw/genomes_redownload", "data/raw/genomes"]
MIN_AA, MAX_AA = 60, 500
MAX_EVALUE = 1e-5


def proteins(path):
    """Identical to 28_allele_search.py, so the protein names line up."""
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
    gf = pyrodigal.GeneFinder(meta=True)
    out = []
    for i, c in enumerate(contigs):
        if len(c) < 200:
            continue
        for j, g in enumerate(gf.find_genes(c)):
            out.append((f"c{i}_g{j}", g.translate().rstrip("*")))
    return out, len(contigs)


def find_file(gid):
    for d in GENOME_DIRS:
        p = os.path.join(ROOT, d, gid + ".fna")
        if os.path.exists(p):
            return p, d
    return None, None


def _name(x):
    """pyhmmer returns names as bytes or str depending on the build."""
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


def load_family_meta():
    meta = {}
    with open(FAMMETA, encoding="utf-8") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            v = dict(zip(cols, line.rstrip("\n").split("\t")))
            meta[v["family"]] = v
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort")
    ap.add_argument("--seconds", type=int, default=150)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(HMMS) as fh:
        hmms = list(fh)
    fam = load_family_meta()
    print(f"  {len(hmms)} profile(s) loaded")

    cohorts = json.load(open(os.path.join(WORK, "cohorts.json")))
    ids = cohorts[a.cohort]
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            try:
                rec = json.loads(line)
                if rec.get("cohort") == a.cohort:
                    done.add(rec["genome_id"])
            except Exception:
                pass
    todo = [g for g in ids if g not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"  {a.cohort}: {len(done)} done, {len(todo)} remaining")

    t0 = time.time()
    fh_out = open(OUT, "a")
    n = 0
    for gid in todo:
        if time.time() - t0 > a.seconds:
            break
        path, src = find_file(gid)
        if path is None:
            fh_out.write(json.dumps({"cohort": a.cohort, "genome_id": gid,
                                     "error": "assembly not found"}) + "\n")
            fh_out.flush()
            n += 1
            continue
        prots, n_contigs = proteins(path)
        cand = [(name, s) for name, s in prots if MIN_AA <= len(s) <= MAX_AA]
        seqs = [pyhmmer.easel.TextSequence(name=name.encode(), sequence=s)
                .digitize(alphabet) for name, s in cand]
        lengths = {name: len(s) for name, s in cand}

        best = {}
        if seqs:
            for top in pyhmmer.hmmsearch(hmms, seqs, E=MAX_EVALUE):
                family = _name(top.query.name)
                for hit in top:
                    if hit.evalue > MAX_EVALUE:
                        continue
                    name = _name(hit.name)
                    cur = best.get(name)
                    if cur is None or hit.score > cur["bit_score"]:
                        f = fam.get(family, {})
                        best[name] = {
                            "protein": name,
                            "aa_length": lengths.get(name),
                            "family": family,
                            "bit_score": round(float(hit.score), 1),
                            "evalue": float(hit.evalue),
                            "family_members": int(f.get("members", 0)),
                            "carbapenemase_members": int(
                                f.get("carbapenemase_members", 0)),
                            "dominant_subclass": f.get("dominant_subclass", ""),
                        }
        hits = sorted(best.values(), key=lambda h: -h["bit_score"])
        fh_out.write(json.dumps({
            "cohort": a.cohort, "genome_id": gid, "source_dir": src,
            "contigs": n_contigs, "proteins": len(prots), "screened": len(cand),
            "profile_hits": len(hits),
            "carbapenemase_family_hits": sum(
                1 for h in hits if h["carbapenemase_members"] > 0),
            "hits": hits}) + "\n")
        fh_out.flush()
        n += 1
    fh_out.close()
    print(f"  processed {n} genome(s) in {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
