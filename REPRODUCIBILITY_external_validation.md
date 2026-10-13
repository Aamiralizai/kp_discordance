# External validation — reproducibility evidence trail

> Every path below is relative to the pipeline root, the directory holding
> `run_all.sh`. In the submission package that is `05_pipeline_github/`, and the
> same files also appear under `04_results_QC/` for readers who open the results
> without the code. Two inputs are regenerated rather than shipped, for size, and
> section 2 says how.

This note lets a reviewer reconstruct the external validation from source.

## 1. Snapshot

- Source: NCBI Pathogen Detection, Klebsiella organism group.
- Snapshot identifier: **PDG000000012.2548**, retrieved 2026-10-12T05:57:55Z.
- Recorded in `results/external/SNAPSHOT.txt`.
- NCBI replaces the `latest` pointer continuously, so the identifier above is
  the only reproducible handle. A later snapshot will not reproduce these exact
  counts.

## 2. Retrieval

- Script: `scripts/19b_fetch_ncbi_pathogen.sh`.
- It resolves the snapshot, downloads
  `Metadata/PDG000000012.2548.metadata.tsv` and keeps only K. pneumoniae rows
  that carry an AST result. The full table is ~190 MB, which is why only the
  filtered file `data/external/ncbi_kp_ast.tsv` is retained rather than
  committed. Re-run the script to regenerate it.
- Field positions are resolved from the header, not hard-coded, because NCBI
  has added columns to this table in the past.

## 3. Overlap exclusion with the derivation set

- The derivation accessions are listed in
  `results/external/overlap_derivation_accessions.tsv`
  (3,347 genomes, with assembly and BioSample accessions).
- `scripts/19_external_validation.py` removes any external isolate whose
  assembly accession or BioSample accession appears in that list.
- The script now **fails** if this file is absent (pass `--allow-no-overlap`
  only to inspect, never to validate), so an external run cannot silently
  skip overlap exclusion.

## 4. Deduplication policy

- Isolates are keyed on **BioSample**. Where one BioSample is represented by
  more than one assembly, the first assembly by accession is kept and the rest
  are collapsed. In this snapshot exactly one BioSample (SAMN02142019) was
  represented by two assemblies, both carbapenemase-negative and susceptible,
  so the choice changes no categorical result.

## 5. Exclusion ledger

`results/tables/external_exclusions.tsv`:

| stage | n |
|---|---:|
| NCBI isolates scanned | 166,142 |
| no meropenem AST | 164,879 |
| overlap with derivation (excluded) | 245 |
| call not S or R (excluded) | 31 |
| duplicate BioSample (collapsed) | 1 |
| **external cohort (final)** | **986** |

## 6. BioProject internal-consistency screen

- `results/tables/external_bioprojects.tsv` lists every BioProject
  with its reported resistance and the fraction of its carbapenemase-positive
  isolates reported resistant.
- A BioProject is flagged when it has at least `--min-support` (10)
  carbapenemase-positive isolates and fewer than `--coherence-pct` (20 percent)
  of them are reported resistant. This is an **outcome-dependent screen**, not
  independent proof that any record is wrong. It is used only to define a
  sensitivity analysis, never to establish error. The full cohort of 986 is the
  primary analysis; the screen flags one BioProject (PRJNA1322038).

## 7. Rule application

- Rules are fitted on the derivation (BV-BRC) set in
  `scripts/17_rule_derivation.py` and written to
  `results/tables/rules_decision_list.tsv`.
- `scripts/19_external_validation.py` applies them unchanged. The external set
  is genotyped with AMRFinderPlus rather than Kleborate, so apparent
  performance is an independent check, and heterogeneity can bias apparent agreement in either direction.
- Outputs, all under `results/tables/`: `external_performance.tsv` (T0/T1/T2,
  full and flagged-excluded, plus leave-one-BioProject-out),
  `external_isolates.tsv` (per isolate), `external_strata.tsv`,
  `external_stratum_comparison.tsv`.

## 8. Clean-run log

`results/logs/external_validation_rerun_*.log` is the stdout of the
rerun that produced the committed tables.

## 9. Versions

Python analysis stack pinned in `requirements.txt` / `environment.yml`
(numpy 2.4.4, scipy 1.15.3, pandas 2.3.3, scikit-learn 1.8.0,
matplotlib 3.10.8, statsmodels 0.14.4). Genotyping: Kleborate 3.2.4
(derivation). External genotypes are AMRFinderPlus as published in the NCBI
snapshot. QC: the completeness, contamination and CDS metrics deposited in BV-BRC
(CheckM-based completeness/contamination, PATRIC CDS count), scored against the
QualiBact thresholds. An independent local CheckM2/Prodigal run via
scripts/00b_qc_assembly.py is an optional alternative route, not the procedure
that generated the committed results.

## 10. One-command reproduction

```
bash scripts/19b_fetch_ncbi_pathogen.sh                 # fetch snapshot
python scripts/19a_fetch_derivation_accessions.py       # overlap manifest
python scripts/19_external_validation.py --drug meropenem
```
