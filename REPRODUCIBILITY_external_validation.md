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

## 6. Record audit against primary sources

NCBI Pathogen Detection is a secondary aggregator. It restates susceptibility
results that submitters deposited on the BioSample record, so the BioSample
antibiogram is the primary source and the aggregated `AST_phenotypes` field is a
copy of it.

`scripts/39_biosample_antibiogram_audit.py` retrieves the deposited antibiogram
for every external isolate and compares every drug both sources report. The
comparison reads the deposited phenotype only. It never reads the carbapenemase
genotype, so it cannot be tuned by the genotype-phenotype agreement this study
measures.

Result, in `results/tables/external_biosample_concordance.tsv`:

| submissions | isolates | comparable calls | drugs | concordant |
| --- | --- | --- | --- | --- |
| 123 | 742 | 12,698 | up to 41 | 12,698 (100 percent) |
| 1 (PRJNA1322038) | 244 | 5,908 | 27 | 0 (0 percent) |

Every comparable call in that one submission is inverted, on all 27 drugs it
reports. Its meropenem phenotype is therefore taken from the deposited record,
which changes 244 of the 986 calls. **No isolate is excluded.** The primary
analysis is the whole cohort of 986.

### Excluding a parsing or field-mapping error

The inversion is asserted on the phenotype label, so three things establish
that the label is not simply being read wrongly.

1. **The same code agrees everywhere else.** It produces complete agreement on
   12,698 calls from 123 submissions. A parser cannot be selectively wrong in
   one submission.
2. **Columns are read by name, not by position.** The antibiogram table's
   columns are looked up by the names BioSample gives them, and the audit stops
   with an error if an expected name is absent or if any row's width disagrees
   with the header. A value therefore cannot be mapped to the wrong field.
3. **Each call is re-derived from the deposited MIC.** The numeric measurement
   is a different field, read by different code that never sees the phenotype
   string, and interpreted against the breakpoint standard the record itself
   declares. A call is returned only where the censoring interval falls wholly
   on one side of a breakpoint.

| cohort | calls with a decidable MIC | MIC agrees with deposited label | MIC agrees with aggregated label |
| --- | --- | --- | --- |
| PRJNA1322038 | 488 | 488 (100.0%) | 0 (0.0%) |
| all other submissions | 1,115 | 1,102 (98.8%) | 1,102 (98.8%) |

In the audited submission the deposited MIC supports the deposited label on
every decidable call and the aggregated label on none. Elsewhere it supports
both equally, because there the two sources agree. The 98.8 percent rather than
100 percent is the expected rate of records sitting at a breakpoint, and it is
worth noting because it shows the control is capable of disagreeing.

Per-record evidence, all under `results/tables/`:

| file | contents |
| --- | --- |
| `external_biosample_corrected_records.tsv` | the 244 corrected records only, with both calls, the deposited MIC and the reason |
| `external_biosample_audit.tsv` | all 986 isolates, with accessions, both calls, MICs, testing standard and decision |
| `external_biosample_concordance.tsv` | one row per submission |
| `external_biosample_parsing_control.tsv` | one row per call with a decidable MIC |

The retrieved XML is cached under `data/biosample_xml/`, so the BioSample side
of the comparison needs no network and reproduces byte-identically. The Pathogen
Detection snapshot is the other side, and it is too large to redistribute, so a
fresh clone must fetch it first with `bash scripts/19b_fetch_ncbi_pathogen.sh`.
The script says so rather than failing obscurely if it is absent.

### What flagged the submission, and why that is not the exclusion criterion

Reported susceptibility to ampicillin identified the submission as worth
auditing, after an unusual pooled result had been traced to it. K. pneumoniae is
intrinsically resistant to ampicillin through its chromosomal SHV
beta-lactamase. EUCAST classes this as an **expected resistant phenotype**,
defined as one in which at least 90 percent of isolates are expected to be
resistant, so the reading is grounds to investigate a record and not proof that
it is wrong. It excludes nothing here, and
`results/tables/external_record_integrity.tsv` reports it as a flag. 256
isolates carry it: 243 in the audited submission and 13 in two submissions the
audit found wholly concordant. Those 13 are retained.

The audit followed the pooled result rather than preceding it, and the
manuscript says so. What makes the correction defensible is not the order in
which it was found but that it rests on 5,908 verified contradictions against
primary records rather than on the quantity under test.

### Sensitivity analyses

| subset | n | basis |
| --- | --- | --- |
| `all_external` | 986 | primary, phenotype from the primary record |
| `all_external_as_deposited` | 986 | aggregated calls unchanged, uncorrected |
| `excluding_audited_submission` | 742 | the audited submission dropped outright |
| `excluding_partial_qc_failures` | 976 | partial QualiBact failures dropped |

Correcting the record and discarding it give the same conclusion, so the result
does not turn on which was done.

### Outcome-dependent screen, retained for reference

`results/tables/external_bioprojects.tsv` also carries the older screen, which
flags a BioProject when it has at least `--min-support` (10) carbapenemase-
positive isolates and fewer than `--coherence-pct` (20 percent) of them are
reported resistant. That screen is **outcome-dependent** and circular, because
it flags a submission by the very quantity under test. It is superseded by the
audit above and is reported for reference only.

## 7. Rule application

- Rules are fitted on the derivation (BV-BRC) set in
  `scripts/17_rule_derivation.py` and written to
  `results/tables/rules_decision_list.tsv`.
- `scripts/19_external_validation.py` applies them unchanged. The external set
  is genotyped with AMRFinderPlus rather than Kleborate, so apparent
  performance is an independent check, and heterogeneity can bias apparent agreement in either direction.
- Outputs, all under `results/tables/`: `external_performance.tsv` (T0/T1/T2
  for each of the four subsets in section 6, plus leave-one-BioProject-out),
  `external_isolates.tsv` (per isolate, carrying `call_as_deposited`,
  `phenotype_source` and `phenotype_corrected`), `external_strata.tsv`,
  `external_stratum_comparison.tsv`, `external_biosample_audit.tsv`,
  `external_biosample_concordance.tsv`.

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
python scripts/19_external_validation.py --drug meropenem --biosample-audit /dev/null
python scripts/39_biosample_antibiogram_audit.py        # audit vs BioSample
python scripts/19_external_validation.py --drug meropenem
python scripts/23_external_summary.py
```

Three passes, because the audit and the validation each need the other's
output. The first settles which isolates are eligible, the audit retrieves and
compares their deposited antibiograms, and the third repeats the validation with
each phenotype taken from the primary record. Pass one's performance tables are
provisional and pass three overwrites them. `run_all.sh` section 15 does this.
