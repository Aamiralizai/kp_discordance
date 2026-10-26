#!/usr/bin/env python3
"""Single definition of what a protein-level hit to a reference allele means.

The first version of this classifier tested coverage before identity, so a
protein matching a catalogued allele at 100 percent identity over 95 percent
of its length was labelled "novel". Two such hits existed, both of them the
genome's own already-reported OXA-48 and KPC-3. Identity decides which allele
a protein is, and coverage decides whether the copy in the assembly is
complete, so identity is tested first here and coverage only afterwards.
"""

EXACT_IDENTITY = 0.999   # at or above this, the protein is the catalogued allele
EXACT_COVERAGE = 0.99    # at or above this, the copy in the assembly is complete
VARIANT_IDENTITY = 0.90  # at or above this, a plausible variant of the allele
HOMOLOGUE_IDENTITY = 0.70


def classify(identity, ref_cov, subclass):
    """Return a verdict for one protein-against-reference alignment.

    exact                  the catalogued allele, complete copy
    partial_catalogued     the catalogued allele, incomplete copy in this assembly
    novel_candidate        high coverage but identity below the allele itself
    partial_divergent      neither complete nor identical to a catalogued allele
    divergent_homologue    a related enzyme, not the same gene
    distant                background similarity
    """
    carb = subclass == "CARBAPENEM"
    if identity >= EXACT_IDENTITY:
        return "exact" if ref_cov >= EXACT_COVERAGE else "partial_catalogued"
    if identity >= VARIANT_IDENTITY:
        if ref_cov >= VARIANT_IDENTITY:
            return "novel_candidate" if carb else "novel_candidate_other_betalactamase"
        return "partial_divergent"
    if identity >= HOMOLOGUE_IDENTITY:
        return "divergent_homologue"
    return "distant"


# verdicts that mean "this genome carries the catalogued carbapenemase"
CATALOGUED = ("exact", "partial_catalogued")
# the only verdict that would answer the editor's question in the affirmative
NOVEL = ("novel_candidate",)
