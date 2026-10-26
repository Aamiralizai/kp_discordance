Frozen inputs for Supplementary Table S17

Table S17 compares every headline estimate of the primary meropenem stratum
computed two ways. One column uses genotypes called on the assembly files as
first retrieved, 3,453 of which were incomplete transfers. The other uses
genotypes called on the re-retrieved and verified assemblies.

The second column comes from the current release. The first column comes from a
superseded results directory produced before the re-retrieval, and
scripts/22_compare_truncation_fix.py needs exactly one file from it:

    manuscript_results.json    the headline estimates of that earlier run

That file is held here so the comparison is regenerable from released files
alone, rather than depending on an archive held outside the release. Verify
with:

    python scripts/22_compare_truncation_fix.py \
        results/results_before_download_repair results

which rewrites results/supplementary/S17_truncation_fix_comparison.tsv. It
reproduces the shipped table byte for byte.

QC_RERUN_NOTE.txt is that run's own note, kept unchanged.

Nothing here is used by any other analysis. The superseded genotype calls are
not part of the study's results and are retained only as the comparison's
input.
