# Carbapenemase genotype–phenotype discordance in *Klebsiella pneumoniae*

Analysis code for a study of discordance between carbapenemase detection and
measured carbapenem susceptibility in publicly available *K. pneumoniae*
genomes.

## What this repository contains

A record-level pipeline that links BV-BRC phenotype records to Kleborate
genotypes, classifies susceptibility from interval-censored MICs, and
quantifies discordance in both directions.

    bash run_all.sh          # clean rebuild of every table, figure and table
    KEEP=1 bash run_all.sh   # keep existing results

The rebuild runs the regression tests first and aborts if any fail, empties
`results/`, and regenerates everything in dependency order. It writes
`results/logs/latest_summary.txt` (a short digest) and
`results/manuscript_results.json`.


## Design rules


1. **Filter phenotype records before aggregating per genome.** Applying
   method or standard restrictions after collapsing to one row per genome
   produces strata that are not what they claim to be.
2. **Reconcile multiple records by interval intersection**, carrying bound
   inclusivity. A record of "> 8" describes (8, ∞) and is resistant under an
   "R > 8" breakpoint. Contradictory records are excluded and reported.
3. **Breakpoints are gated on verification.** `config/breakpoints.csv` stores
   value *and* comparator *and* source document; analysis refuses to run on an
   unverified row. CLSI defines resistance as MIC ≥ R, EUCAST as MIC > R.
4. **One contrast throughout.** Fisher tests, regression and FDR all compare
   intact porins with OmpK36-level disruption, excluding OmpK35-only isolates
   from both arms.
5. **Rates use the binary S/R set.** Intermediate and unresolvable isolates are
   excluded from denominators and their counts reported.
6. **Check convergence.** Non-converged models are refitted by Firth penalised
   regression or reported as not estimable — never printed as if valid.

## Layout

    config/breakpoints.csv    verified breakpoints (value, comparator, source)
    config/flow_counts.yaml   assembly-retrieval counts for Figure 1
    scripts/                  pipeline, numbered in execution order
    scripts/kp_defs.py        single source of genotype definitions
    tests/                    regression tests — run before any rebuild
    results/                  generated; not committed
    manuscript/               manuscript, supplementary, reference library
    archive/deprecated/       superseded code retained for provenance only

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
