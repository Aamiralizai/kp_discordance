#!/usr/bin/env python3
"""Locate each carbapenemase-like gene fragment on its contig.

A fragment that abuts a contig boundary is an assembly break, not a
biological truncation. This script reports, for every flagged genome, the
contig length, the ORF coordinates, the partial flags pyrodigal assigns and
the distance from the ORF to each contig end, so the two cases can be told
apart from the output alone.

Resumable: already-written genomes in the output file are skipped.
"""
import json, os, sys

import pyrodigal
import pyswrd

# the scripts live in scripts/ and their data in results/allele_search/, so
# both are resolved from the repository root rather than from the file's own
# directory, which lets the script be moved without silently writing elsewhere
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "results", "allele_search")
GENOME_DIRS = ["data/raw/genomes_redownload", "data/raw/genomes"]
MIN_IDENTITY = 0.70
MIN_ALIGNED_AA = 60  # a local match shorter than this says nothing about the gene
OUT = os.path.join(WORK, "fragment_loci.jsonl")
REFS = os.path.join(WORK, "ref", "AMRProt.fa")  # symbol is field 3 of the pipe-delimited header

# genome -> reference symbols flagged by 28/29_*.py
TARGETS = {
    "573.13799": ["blaKPC-2"],
    "573.14223": ["blaKPC-211", "blaKPC-241"],
    "573.21572": ["blaNDM-1"],
    "573.21677": ["blaNDM-1"],
    "72407.559": ["blaCMH-3"],
    "573.20501": ["blaOXA-58"],
    "573.21958": ["blaKPC-285"],
    "573.21973": ["blaVIM-1"],
    "72407.519": ["blaOXA-233"],
}


def load_refs(wanted):
    seqs, name = {}, None
    buf = []
    with open(REFS) as fh:
        for line in fh:
            if line.startswith(">"):
                if name and name in wanted:
                    seqs[name] = "".join(buf)
                parts = line[1:].strip().split("|")
                name = parts[3] if len(parts) > 3 else parts[0]
                buf = []
            else:
                buf.append(line.strip())
    if name and name in wanted:
        seqs[name] = "".join(buf)
    return seqs


def contigs(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def find(genome, symbols, refs):
    path = None
    for d in GENOME_DIRS:
        p = os.path.join(ROOT, d, genome + ".fna")
        if os.path.exists(p):
            path = p
            break
    if path is None:
        return []
    finder = pyrodigal.GeneFinder(meta=True)
    queries, meta = [], []
    for cname, seq in contigs(path):
        if len(seq) < 60:
            continue
        for g in finder.find_genes(seq):
            aa = g.translate().rstrip("*")
            if not 50 <= len(aa) <= 600:
                continue
            queries.append(aa)
            meta.append((cname, len(seq), g.begin, g.end, g.strand,
                         bool(g.partial_begin), bool(g.partial_end)))
    targets = [refs[s] for s in symbols if s in refs]
    names = [s for s in symbols if s in refs]
    out = []
    for hit in pyswrd.search(queries, targets, gap_open=11, gap_extend=1,
                             scorer_name="BLOSUM62"):
        r = hit.result
        ident = r.identity()
        aligned = r.target_end - r.target_start + 1
        if ident < MIN_IDENTITY or aligned < MIN_ALIGNED_AA:
            continue
        cname, clen, begin, end, strand, pb, pe = meta[hit.query_index]
        out.append(dict(
            genome=genome, ref=names[hit.target_index],
            identity=round(ident, 4), aligned_aa=aligned,
            contig=cname, contig_len=clen,
            orf_start=begin, orf_end=end, strand=strand,
            partial_begin=pb, partial_end=pe,
            bases_to_contig_start=begin - 1,
            bases_to_contig_end=clen - end,
            ref_residues=f"{r.target_start + 1}-{r.target_end + 1}",
            ref_len=r.target_length,
            assembly_break=bool(pb or pe or begin - 1 <= 3 or clen - end <= 3),
        ))
    return out


def main():
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            done.add(json.loads(line)["genome"])
    refs = load_refs({s for v in TARGETS.values() for s in v})
    missing = {s for v in TARGETS.values() for s in v} - set(refs)
    if missing:
        sys.stderr.write("references not found: %s\n" % sorted(missing))
    todo = [g for g in TARGETS if g not in done]
    budget = int(os.environ.get("BUDGET", "3"))
    with open(OUT, "a") as fh:
        for genome in todo[:budget]:
            rows = find(genome, TARGETS[genome], refs)
            if not rows:
                rows = [dict(genome=genome, ref=None,
                             note="no hit passing %.2f identity over %d aa"
                                  % (MIN_IDENTITY, MIN_ALIGNED_AA))]
            for r in rows:
                fh.write(json.dumps(r) + "\n")
            fh.flush()
            print("done", genome, len(rows), flush=True)
    left = len(todo) - min(budget, len(todo))
    print("remaining:", left)
    return 0


if __name__ == "__main__":
    sys.exit(main())
