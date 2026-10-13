#!/usr/bin/env python3
"""28_allele_search.py — look for beta-lactamase alleles the databases do not call.

The editor asked whether the discordant isolates carry putative carbapenemase
alleles, or other mutations, that current AMR gene databases do not capture.
Kleborate reports only what matches its catalogue, so a divergent or disrupted
allele is invisible to it. This searches the assemblies directly.

Two stages, because a full search of every predicted protein against every
reference allele takes about two minutes per genome.

  stage 1  every predicted protein of plausible length is searched against one
           representative per beta-lactamase gene family, which finds the few
           proteins worth looking at properly
  stage 2  those proteins are searched against the full reference set, giving
           the closest allele with its identity and its coverage

Each protein is then classified.

  exact        identical to a catalogued allele over its whole length
  novel        at least 90 percent identical to a carbapenemase but not equal
               to any catalogued allele, which is a putative new allele
  disrupted    at least 95 percent identical but covering less than 90 percent
               of the reference, which is a truncated or interrupted gene
  divergent    between 70 and 90 percent identical, a relative rather than an
               allele

The run is resumable. Each genome is appended to a JSONL checkpoint as it
finishes, and a second invocation skips what is already there, so the work can
be done in short sessions.

USAGE
    python3 28_allele_search.py COHORT [--seconds 150] [--limit N]
"""
import json
import os
import sys
import time

import pyrodigal
import pyswrd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from allele_classify import classify  # noqa: E402

# the scripts live in scripts/ and their data in results/allele_search/, so
# both are resolved from the repository root rather than from the file's own
# directory, which lets the script be moved without silently writing elsewhere
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "results", "allele_search")
REF = os.path.join(WORK, "ref", "AMRProt.fa")
OUT = os.path.join(WORK, "hits.jsonl")
GENOME_DIRS = ["data/raw/genomes_redownload", "data/raw/genomes"]
MIN_AA, MAX_AA = 180, 500


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
                         "subclass": f[7], "desc": f[9] if len(f) > 9 else "",
                         "seq": s})
    best = {}
    for i, r in enumerate(refs):
        cur = best.get(r["family"])
        if cur is None or len(r["seq"]) > len(refs[cur]["seq"]):
            best[r["family"]] = i
    return refs, sorted(best.values())


def proteins(path):
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
            out.append((f"c{i}_g{j}", g.translate().rstrip("*"), len(c)))
    return out, len(contigs)


def find_file(gid):
    for d in GENOME_DIRS:
        p = os.path.join(ROOT, d, gid + ".fna")
        if os.path.exists(p):
            return p, d
    return None, None




def main():
    cohort = sys.argv[1]
    secs = 150
    limit = None
    if "--seconds" in sys.argv:
        secs = int(sys.argv[sys.argv.index("--seconds") + 1])
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])

    refs, rep_idx = load_reference(REF)
    rep_seqs = [refs[i]["seq"] for i in rep_idx]
    all_seqs = [r["seq"] for r in refs]

    cohorts = json.load(open(os.path.join(WORK, "cohorts.json")))
    ids = cohorts[cohort]

    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            try:
                rec = json.loads(line)
                if rec.get("cohort") == cohort:
                    done.add(rec["genome_id"])
            except Exception:
                pass
    todo = [g for g in ids if g not in done]
    if limit:
        todo = todo[:limit]
    print(f"  {cohort}: {len(done)} done, {len(todo)} remaining")

    t_start = time.time()
    fh = open(OUT, "a")
    n = 0
    for gid in todo:
        if time.time() - t_start > secs:
            break
        path, src = find_file(gid)
        if path is None:
            fh.write(json.dumps({"cohort": cohort, "genome_id": gid,
                                 "error": "assembly not found"}) + "\n")
            fh.flush()
            n += 1
            continue
        prots, n_contigs = proteins(path)
        cand = [(name, s) for name, s, _ in prots if MIN_AA <= len(s) <= MAX_AA]
        keep = set()
        if cand:
            for hit in pyswrd.search([s for _, s in cand], rep_seqs,
                                     max_alignments=1, max_evalue=1e-10):
                keep.add(hit.query_index)
        found = []
        if keep:
            sel = [cand[i] for i in sorted(keep)]
            best = {}
            for hit in pyswrd.search([s for _, s in sel], all_seqs,
                                     max_alignments=5, max_evalue=1e-5):
                r = hit.result
                ident = r.identity()
                ref = refs[hit.target_index]
                cov = (r.target_end - r.target_start + 1) / len(ref["seq"])
                cur = best.get(hit.query_index)
                if cur is None or ident > cur["identity"]:
                    best[hit.query_index] = {
                        "protein": sel[hit.query_index][0],
                        "aa_length": len(sel[hit.query_index][1]),
                        "ref_symbol": ref["symbol"], "ref_family": ref["family"],
                        "ref_acc": ref["acc"], "subclass": ref["subclass"],
                        "ref_aa_length": len(ref["seq"]),
                        "identity": round(ident, 4),
                        "ref_coverage": round(cov, 4),
                    }
            for v in best.values():
                v["verdict"] = classify(v["identity"], v["ref_coverage"],
                                        v["subclass"])
                found.append(v)
        fh.write(json.dumps({"cohort": cohort, "genome_id": gid,
                             "source_dir": src, "contigs": n_contigs,
                             "proteins": len(prots), "screened": len(cand),
                             "betalactamase_like": len(found),
                             "hits": found}) + "\n")
        fh.flush()
        n += 1
    fh.close()
    print(f"  processed {n} genome(s) in {time.time() - t_start:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
