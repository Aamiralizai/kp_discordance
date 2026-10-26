#!/usr/bin/env python3
"""34_build_betalactamase_profiles.py — family profiles from the reference set.

The allele search compares each predicted protein with reference sequences one
pair at a time. Pairwise alignment is the weaker of the two methods a catalogue
search normally uses. AMRFinderPlus itself combines it with profile hidden
Markov models, and a profile built from a whole family detects a remote member
that no single pairwise comparison would reach, because it knows which
positions a family conserves and which it tolerates.

This builds one profile per beta-lactamase family in the reference set that has
at least three members, using the reference set already shipped with the
release. Nothing is downloaded, so the profiles carry the same provenance and
the same checksum as the search they complement.

USAGE
    python3 scripts/34_build_betalactamase_profiles.py
"""
import collections
import os
import sys

import pyfamsa
import pyhmmer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "results", "allele_search")
REF = os.path.join(WORK, "ref", "AMRProt.fa")
OUT = os.path.join(WORK, "ref", "betalactamase_families.hmm")
META = os.path.join(WORK, "ref", "betalactamase_families.tsv")
MIN_MEMBERS = 3
# the widest member length spread a family may show and still give a profile
# specific to one enzyme family. One of the 136 families exceeds it.
SPECIFIC_MAX_RATIO = 1.25


def load():
    recs, head, seq = [], None, []
    for line in open(REF):
        if line.startswith(">"):
            if head:
                recs.append((head, "".join(seq)))
            head, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if head:
        recs.append((head, "".join(seq)))
    out = []
    for h, s in recs:
        f = h.split("|")
        if len(f) > 8 and f[8] == "BETA-LACTAM":
            out.append({"acc": f[0], "symbol": f[3], "family": f[4],
                        "subclass": f[7], "seq": s})
    return out


def main():
    refs = load()
    fams = collections.defaultdict(list)
    for r in refs:
        fams[r["family"]].append(r)
    print(f"  {len(refs)} beta-lactamase reference(s) in {len(fams)} family(ies)")

    # one thread, so the alignment and therefore the profile are the same
    # on any machine. The .hmm file still carries a build date in its
    # header, so reproducibility is checked against the TSV beside it.
    aligner = pyfamsa.Aligner(guide_tree="upgma", threads=1)
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)

    profiles, meta, skipped = [], [], 0
    for family, members in sorted(fams.items()):
        if len(members) < MIN_MEMBERS:
            skipped += 1
            continue
        seqs = [pyfamsa.Sequence(m["acc"].encode(), m["seq"].encode())
                for m in members]
        msa = aligner.align(seqs)
        rows_aln = [s.sequence.decode() for s in msa]
        # the alignment width is taken from the aligner's own output. Reading it
        # back from the TextMSA reported a different number, which put a wrong
        # alignment width in this table and in the Supplementary caption that
        # cites it. The model length below comes from the built profile, so the
        # two are now measured independently and can be compared.
        width = len(rows_aln[0]) if rows_aln else 0
        if any(len(r) != width for r in rows_aln):
            sys.exit(f"{family}: the aligner returned rows of unequal length")
        text = pyhmmer.easel.TextMSA(
            name=family.encode(),
            sequences=[pyhmmer.easel.TextSequence(name=s.id, sequence=s.sequence.decode())
                       for s in msa],
        )
        hmm, _, _ = builder.build_msa(text.digitize(alphabet), background)
        if hmm.M > width:
            sys.exit(f"{family}: model length {hmm.M} exceeds the alignment "
                     f"width {width}, which cannot happen")
        profiles.append(hmm)
        subclasses = collections.Counter(m["subclass"] for m in members)
        lens = [len(m["seq"]) for m in members]
        ratio = max(lens) / min(lens)
        meta.append({
            "family": family,
            "members": len(members),
            "aligned_columns": width,
            "model_length": hmm.M,
            "member_aa_min": min(lens),
            "member_aa_max": max(lens),
            "member_length_ratio": round(ratio, 3),
            # A family whose members differ widely in length does not align
            # cleanly, so its profile matches on a shared fold rather than on
            # a specific enzyme and its hits carry little weight. The flag is
            # recorded rather than used to drop the profile, so the screen
            # covers every family and the summary can report the two kinds of
            # hit separately.
            "profile_specific": ratio <= SPECIFIC_MAX_RATIO,
            "dominant_subclass": subclasses.most_common(1)[0][0],
            "carbapenemase_members": subclasses.get("CARBAPENEM", 0),
            "example_symbol": members[0]["symbol"],
        })

    with open(OUT, "wb") as fh:
        for hmm in profiles:
            hmm.write(fh)
    with open(META, "w", encoding="utf-8") as fh:
        cols = list(meta[0])
        fh.write("\t".join(cols) + "\n")
        for m in meta:
            fh.write("\t".join(str(m[c]) for c in cols) + "\n")

    carb = sum(1 for m in meta if m["carbapenemase_members"] > 0)
    print(f"  {len(profiles)} profile(s) built, {skipped} family(ies) below "
          f"{MIN_MEMBERS} members and skipped")
    print(f"  {carb} profile(s) contain at least one carbapenemase reference")
    print(f"[out] {OUT}")
    print(f"[out] {META}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
