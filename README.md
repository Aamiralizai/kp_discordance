# Carbapenemase genotype–phenotype discordance in *Klebsiella pneumoniae*

Analysis code for a study of discordance between carbapenemase detection and
measured carbapenem susceptibility in publicly available *K. pneumoniae*
genomes.

## What this repository contains

A record-level pipeline that links BV-BRC phenotype records to Kleborate
genotypes, classifies susceptibility from interval-censored MICs, and
quantifies discordance in both directions.

It carries the code, the frozen inputs the pipeline reads, and the result
files behind every number in the paper. `results/RELEASE_MANIFEST.sha256`
holds a SHA-256 for each of them, so anything here can be checked without
rerunning the analysis:

    cd <repository root>
    sha256sum -c results/RELEASE_MANIFEST.sha256

`results/RELEASE_MANIFEST.txt` records the software environment, the outcome
of every pipeline stage, and the size and digest of every tracked file. Two
retrieved inputs are not redistributed because of their size, the assembly
collection and the filtered NCBI Pathogen Detection table. The scripts named
below fetch both from pinned sources, so a replay from sequence to result
needs network access while everything downstream of the frozen tables does
not.

    bash run_all.sh           
    KEEP=1 bash run_all.sh   # keep existing results

The pipeline runs the regression tests first and aborts if any fail, empties
`results/`, and regenerates everything in dependency order. It writes
`results/logs/latest_summary.txt` (a short digest) and
`results/manuscript_results.json`.


## Design rules


0. **A downloaded assembly is not trustworthy until it is checked against the
   deposited record.** This project learned it the hard way. 3,453 of 7,360
   assemblies (46.9%) were truncated FTPS transfers: the download stopped
   mid-sequence and the client still exited 0. Each file was a parseable FASTA,
   so the structural check passed, and the deposited QualiBact metrics passed
   too, because those describe the assembly as submitted to the archive and not
   the copy on disk. Kleborate typed the incomplete files and missed the genes
   on the contigs that never arrived, which understated carbapenemase and ESBL
   detection and inverted one reported association.

   `00f_check_download_integrity.py` compares every local FASTA with the
   BV-BRC genome_length and contig count and flags a missing final newline.
   `01d_redownload_verified.py` re-retrieves what it finds and accepts a file
   only on an exact match. `04c_fix_truncated_genomes.sh` re-types the repaired
   assemblies with the same genotyper call. Run 00f before believing any
   downstream number. The evidence for this project is in
   `results/qc_concordance/`.

1. **Structural validity and assembly quality are two different checks.**
   `00_qc_fasta.sh` asks whether a file is a parseable FASTA, because the
   genotyper aborts a whole batch on a malformed record. `00b_qc_assembly.py`
   asks whether the assembly is good enough to genotype, against the
   species-specific QualiBact v1.0 thresholds in `config/qc_thresholds.csv`.
   Only the second bears on the study's validity: a fragmented assembly can
   truncate *ompK36*, which is read as porin disruption, or drop a
   carbapenemase gene entirely. Assembly quality is confounded with both
   exposures, in opposite directions, so it is controlled before typing and
   tested as a competing explanation afterwards.

   The committed QC table (results_regenerated/qc/assembly_qc_metrics.tsv)
   was produced by 00e_fetch_bvbrc_qc.py, i.e. the CheckM/PATRIC quality
   metrics and assembly statistics deposited in BV-BRC, scored against the
   QualiBact thresholds. It is not an independent local CheckM2/Prodigal run.
   00b_qc_assembly.py computes the same metrics from local FASTAs and CheckM2
   for anyone who wishes to recompute them from sequence.

2. **Filter phenotype records before aggregating per genome.** Applying
   method or standard restrictions after collapsing to one row per genome
   produces strata that are not what they claim to be.
3. **Reconcile multiple records by interval intersection**, carrying bound
   inclusivity. A record of "> 8" describes (8, ∞) and is resistant under an
   "R > 8" breakpoint. Contradictory records are excluded and reported.
4. **Breakpoints are gated on verification.** `config/breakpoints.csv` stores
   value *and* comparator *and* source document; analysis refuses to run on an
   unverified row. CLSI defines resistance as MIC ≥ R, EUCAST as MIC > R.
5. **One contrast throughout.** Fisher tests, regression and FDR all compare
   intact porins with OmpK36-level disruption, excluding OmpK35-only isolates
   from both arms.
6. **Rates use the binary S/R set.** Intermediate and unresolvable isolates are
   excluded from denominators and their counts reported.
7. **Check convergence.** Non-converged models are refitted by Firth penalised
   regression or reported as not estimable — never printed as if valid.

## Layout

    config/breakpoints.csv    verified breakpoints (value, comparator, source)
    config/qc_thresholds.csv  assembly QC thresholds (QualiBact, EnteroBase)
    config/flow_counts.yaml   assembly-retrieval counts for Figure 1
    assets/Figure1_study_flow_template.svg
                              Figure 1's layout, with {{PLACEHOLDER}} marks
                              where the counts go
    scripts/                  pipeline, numbered in execution order
    scripts/00b_qc_assembly.py      assembly QC from local FASTA
    scripts/00e_fetch_bvbrc_qc.py   assembly QC from BV-BRC (no download)
    scripts/05b_apply_qc_filter.py  removes QC-failed genomes from genotypes
    scripts/kp_defs.py        single source of genotype definitions
    scripts/12_flow_diagram.py       Figure 1 drawn by matplotlib
    scripts/12b_flow_diagram_typeset.py
                              Figure 1 again, filled from the layout template.
                              Both read the same tables, so the typeset figure
                              cannot drift from the analysis; it is the one
                              submitted, because matplotlib writes italic text
                              one character per element and cannot be corrected
                              by hand afterwards
    scripts/render_svg.py     SVG to PNG, PDF and TIFF, with whichever renderer
                              the machine has
    tests/                    regression tests, run before any rebuild
    tools/accept_changes.py   accepts the tracked changes in a manuscript file.
                              Runs no analysis
    data/                     the frozen inputs the pipeline reads
    results/                  the outputs behind every number in the paper
    results/figures/          the submitted figures, SVG with editable text,
                              PDF and PNG. The only copy, so it cannot go stale
    references/               the reference library, BibTeX and RIS
    results/allele_search/    the uncatalogued-allele search: its reference
                              set, its per-genome hits, and PROVENANCE.txt,
                              which records that its gene-calling stages ran
                              outside the pinned environment and why that does
                              not affect the tables

## Predictive rules and external validation

    scripts/17_rule_derivation.py        nested rule tiers, predictive values,
                                         clone-aware cross-validation
    scripts/18_discordance_investigation.py  what the discordant isolates are
    scripts/19_external_validation.py    the rules applied, unchanged, to NCBI
                                         Pathogen Detection isolates

Two conventions matter when reading the outputs.

Cross-validated performance is reported with folds split by sequence type, not
at random. Public genome collections are heavily clonal, and a rule scored on
near-identical isolates of a clone it was fitted on reports its own training
data back. The ST-grouped figure is the one to quote; the gap between it and
the random-fold figure measures clonal redundancy.

The external set is screened by submission before anything is pooled. A single
BioProject can be a quarter of the validation set and can carry AST metadata
that contradicts its own genotypes. Each BioProject is therefore checked for
internal coherence, results are reported with and without any that fail, and a
leave-one-BioProject-out range shows that no single submission determines the
answer.

## Search for uncatalogued carbapenemase alleles

    scripts/allele_classify.py           what a hit to a reference allele means
    scripts/28_allele_search.py          proteins of 180 to 500 residues
    scripts/29_short_orf_check.py        open reading frames of 60 to 180
    scripts/30_fragment_loci.py          where each hit sits on its contig
    scripts/31_allele_search_summary.py  Supplementary Table S19a
    scripts/32_allele_search_tables.py   Supplementary Tables S19b and S19c

Two length windows are screened, and the second one is the point. A complete
carbapenemase is 265 to 300 residues, so a window of 180 to 500 finds an
uncatalogued enzyme but excludes every gene an assembly has broken, which is
one of the two things the search exists to detect. The first version of this
analysis ran only that window and reported nothing. The 60 to 180 window found
four genomes carrying a carbapenemase at 94.9 to 100 percent identity to a
catalogued allele, each fragment abutting a contig boundary.

A hit is therefore not interpreted from identity alone. Every alignment is
located on its contig, and an open reading frame that abuts a boundary or that
the gene caller marks incomplete is recorded as an assembly break rather than a
biological truncation. The distinction is the whole result, so the coordinates
are published in Supplementary Table S19b rather than summarised.

Identity is tested before coverage, in scripts/allele_classify.py. Testing them
together, as the first version did, labels a protein matching a catalogued
allele at 100 percent identity over 95 percent of its length as a novel allele,
which is the one verdict this analysis must not get wrong. Identity decides
which allele a protein is. Coverage decides whether the copy in the assembly is
complete.

The table scripts reclassify from the stored identity and coverage rather than
from the verdict written during the run, so correcting the classifier does not
require repeating a search that takes hours.

## Verified breakpoints

| Standard | Meropenem | Imipenem | Source |
|---|---|---|---|
| CLSI | S ≤ 1, R ≥ 4 | S ≤ 1, R ≥ 4 | M100 36th ed. (2026), Table 2A-1 |
| EUCAST | S ≤ 2, R > 8 | S ≤ 2, R > 4 | v16.1, valid 24 June 2026 |

## Environment

    conda env create -f environment.yml
    conda activate kp-discordance
    # genotyping runs separately:
    conda create -n kleborate -c bioconda -c conda-forge kleborate=3.2.4

## Prerequisites

`run_all.sh` assumes the slow upstream stages have already been run:

    scripts/01_download_bvbrc.sh    phenotype records and assemblies
    scripts/00_qc_fasta.sh         assembly quality control
    scripts/04_run_kleborate.sh    genotyping
    scripts/05_merge_kleborate.py  merge shards
