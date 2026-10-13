#!/usr/bin/env bash
# run_all.sh — clean rebuild of every analysis output.
#
# Empties results/, regenerates everything in dependency order from the
# record-level pipeline, and stops at the first failure. Nothing downstream
# ever reads a stale or legacy file.
#
# PREREQUISITES (not run here — they are slow and rarely change):
#   scripts/01_download_bvbrc.sh      phenotype records + assemblies
#   scripts/00_qc_fasta.sh            structural FASTA check (so typing runs)
#   scripts/00c_run_checkm2.sh        completeness + contamination
#   scripts/00b_qc_assembly.py        assembly QUALITY control (QualiBact)
#
#   No-download alternatives to 00b+00c, when the assemblies came from a public
#   archive that already publishes the metrics:
#   scripts/00e_fetch_bvbrc_qc.py     BV-BRC (preferred: complete, and it
#                                     describes the assembly that was typed)
#   scripts/00d_fetch_ncbi_qc.py      NCBI Datasets (partial coverage)
#   scripts/04_run_kleborate.sh       genotyping
#   scripts/05_merge_kleborate.py     merge shards
#   scripts/05b_apply_qc_filter.py    drop QC-failed genomes from the genotypes
#
# The two QC steps are different and both are required. 00_qc_fasta.sh checks
# that a file is a valid FASTA, because Kleborate aborts a whole batch on a
# malformed record. 00b_qc_assembly.py checks whether the assembly is good
# enough to genotype, against species-specific QualiBact thresholds. Only the
# second one bears on whether ompK36 truncation and carbapenemase absence are
# real.
#
# 05b is not optional. Quarantining an assembly moves a file; it does not
# remove that genome from the merged genotype table, which is what the analysis
# actually reads. Skip 05b and the QC step changes nothing downstream.
#
# EXTERNAL VALIDATION (also prerequisites, network-bound):
#   scripts/19b_fetch_ncbi_pathogen.sh         validation metadata
#   scripts/19a_fetch_derivation_accessions.py accessions for overlap exclusion
#
# Requires verified breakpoints in config/breakpoints.csv.
#
# Usage:
#   bash run_all.sh              # full rebuild
#   KEEP=1 bash run_all.sh       # keep existing results (no wipe)

set -uo pipefail          # note: not -e; failures are handled per-step
: "${KEEP:=0}"

FAILED=()
STEP_STATUS=()

trap 'echo "\n!! run_all.sh aborted unexpectedly (line $LINENO)"; write_summary; exit 1' ERR

say() { printf '\n\033[1m== %s ==\033[0m\n' "$1"; }
need() { [ -f "$1" ] || { echo "MISSING REQUIRED INPUT: $1"; exit 1; }; }

WARN=0
# Primary analyses must succeed. Secondary ones (a drug or stratum with no
# usable records) warn and continue, so one empty subset cannot abort the
# whole rebuild — but the omission is reported at the end and the results
# manifest lists what was not generated.
must() {
  if "$@"; then
    STEP_STATUS+=("OK   $*")
  else
    STEP_STATUS+=("FAIL $*")
    FAILED+=("$*")
    echo "   !! REQUIRED STEP FAILED: $*"
    write_summary
    exit 1
  fi
}

opt() {
  if "$@"; then
    STEP_STATUS+=("OK   $*")
  else
    STEP_STATUS+=("SKIP $*")
    echo "   !! SKIPPED (no usable data or analysis failed): $*"
    WARN=$((WARN+1))
    SKIPPED_STEPS+=("$*")
  fi
}

# A run that skipped a step is not a full run. The summary used to record
# SKIP in a list most readers do not reach, so a reported "complete rebuild"
# could silently be missing the external validation. It now says so at the end
# and exits non-zero unless the caller sets ALLOW_SKIPPED=1.
SKIPPED_STEPS=()

# A short, shareable summary: environment, step outcomes, headline numbers and
# the generated files. This is the file to send when asking about a run.
write_summary() {
  {
    echo "kp-multiscale-amr — run summary"
    echo "timestamp   : $STAMP"
    echo "host        : $(hostname 2>/dev/null || echo unknown)"
    echo "python      : $(python -V 2>&1)"
    echo "conda env   : ${CONDA_DEFAULT_ENV:-<none>}"
    echo "full log    : $LOG"
    echo
    echo "--- verified breakpoints ---"
    if [ -f config/breakpoints.csv ]; then
      awk -F, 'NR>1 && $7=="YES" {printf "  %-12s %-7s S%s%-5s R%s%-5s  %s\n", $1,$2,$3,$4,$5,$6,$8}' \
        config/breakpoints.csv 2>/dev/null | head -12
    else
      echo "  config/breakpoints.csv NOT FOUND"
    fi
    echo
    echo "--- step outcomes ---"
    for s in "${STEP_STATUS[@]:-}"; do echo "  $s"; done
    echo
    if [ "${#FAILED[@]}" -gt 0 ]; then
      echo "--- FAILURES ---"
      for f in "${FAILED[@]}"; do echo "  $f"; done
      echo
      echo "--- last 25 log lines ---"
      tail -25 "$LOG"
      echo
    fi
    echo "--- headline numbers (from manuscript_results.json) ---"
    if [ -f results/manuscript_results.json ]; then
      python - <<'PYEOF' 2>/dev/null || echo "  (could not parse manifest)"
import json
d = json.load(open("results/manuscript_results.json"))
for a in d.get("analyses", []):
    p = a["carbapenemase_positive_susceptible_pct_ci"]
    print(f"  {a['analysis']:<38} binary {a['binary_set']:>6,}  "
          f"gene+ susceptible {p[0]}% [{p[1]}-{p[2]}]")
if d.get("missing"):
    print("  NOT GENERATED: " + ", ".join(d["missing"]))
PYEOF
    else
      echo "  manuscript_results.json NOT GENERATED"
    fi
    echo
    echo "--- generated files ---"
    echo "  tables       : $(ls results/tables/*.tsv 2>/dev/null | wc -l)"
    echo "  figures (svg): $(ls results/figures/*.svg 2>/dev/null | wc -l)"
    echo "  supplementary: $(ls results/supplementary/* 2>/dev/null | wc -l)"
    ls results/figures/*.svg 2>/dev/null | sed 's|.*/|    |'
  } > "$SUMMARY"
  echo
  echo "summary written to $SUMMARY"
}

# The exit status of a rebuild. A required stage that did not run — because it
# failed, or because its input was absent — makes the run partial, and a
# partial run must not exit 0, or it will be reported as a complete rebuild.
# ALLOW_SKIPPED=1 is the only way to call a partial run deliberate.
finish_gate() {
  if [ "${#SKIPPED_STEPS[@]}" -gt 0 ]; then
    echo
    echo "INCOMPLETE RUN: ${#SKIPPED_STEPS[@]} step(s) did not run:"
    for st in "${SKIPPED_STEPS[@]}"; do echo "  - $st"; done
    echo "Results from this run are NOT a full rebuild."
    echo "Do not describe this run as complete, and do not cite its outputs as"
    echo "a regenerated result for any stage listed above."
    [ "${ALLOW_SKIPPED:-0}" = "1" ] || exit 3
    echo "ALLOW_SKIPPED=1 was set: exiting 0 for a deliberate partial run."
  fi
  return 0
}

# tests/test_run_all_gate.sh sources this file to exercise finish_gate without
# running the pipeline. Everything above this line only defines things.
if [ "${RUN_ALL_LIB_ONLY:-0}" = "1" ]; then
  return 0 2>/dev/null || exit 0
fi

# ------------------------------------------------- SIDE EFFECTS START HERE
# Logging was set up at the top of the script until the gate test began
# sourcing this file. Sourcing then ran this block too, stamping a second,
# empty log and repointing latest.log at it in the middle of the very run that
# was writing the real one. Nothing below the guard runs when the file is
# sourced, so the stamp, the log and the pointers belong here.
STAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p results/logs
LOG="results/logs/run_${STAMP}.log"
SUMMARY="results/logs/summary_${STAMP}.txt"

# everything printed from here on goes to the terminal AND the log
exec > >(tee -a "$LOG") 2>&1

# ---------------------------------------------------------------- checks
need config/breakpoints.csv
need data/processed/kleborate_klebsiella_pneumo_complex_output.tsv
need data/raw/phenotypes/kp_amr_anchor.tsv

# If assembly QC has been run, the genotype table must have been filtered
# through 05b. Otherwise the analysis silently includes genomes that failed QC.
if [ -f data/interim/assembly_qc_metrics.tsv ] \
   && [ ! -f data/processed/kleborate_klebsiella_pneumo_complex_output.prefilter.tsv ]; then
  echo "ERROR: assembly QC metrics exist but the genotype table was never filtered."
  echo "       Run: python scripts/05b_apply_qc_filter.py"
  echo "       (quarantining assemblies does not remove them from the genotype table)"
  exit 1
fi

echo "kp-multiscale-amr — clean rebuild"
echo "started $(date)  |  log: $LOG"
echo "python $(python -V 2>&1)  |  env ${CONDA_DEFAULT_ENV:-<none>}"

say "0. regression tests (must pass before anything is generated)"
must python tests/test_classification.py
must bash tests/test_run_all_gate.sh

if [ "$KEEP" != "1" ]; then
  say "1. clearing previous outputs"
  rm -rf results/tables results/figures results/supplementary results/manuscript_results.json
  mkdir -p results/tables results/figures results/supplementary
  echo "   results/ emptied — no legacy file can be read downstream"
else
  echo "   KEEP=1: existing results retained"
  mkdir -p results/tables results/figures results/supplementary
fi

# ---------------------------------------------------------------- analyses
say "2. record-level analyses — CLSI (primary)"
must python scripts/08b_discordance_recordlevel.py --drug meropenem \
    --standard CLSI --method "Broth dilution" --breakpoints CLSI --tag _clsi_bmd
opt  python scripts/08b_discordance_recordlevel.py --drug imipenem \
    --standard CLSI --method "Broth dilution" --breakpoints CLSI --tag _clsi_bmd

say "3. record-level analyses — EUCAST reclassification"
opt python scripts/08b_discordance_recordlevel.py --drug meropenem \
    --standard CLSI --method "Broth dilution" --breakpoints EUCAST --tag _clsi_bmd_eucast
opt python scripts/08b_discordance_recordlevel.py --drug imipenem \
    --standard CLSI --method "Broth dilution" --breakpoints EUCAST --tag _clsi_bmd_eucast

say "4. record-level deduplication (clonality sensitivity)"
opt python scripts/08b_discordance_recordlevel.py --drug meropenem \
    --standard CLSI --method "Broth dilution" --breakpoints CLSI \
    --dedup --tag _clsi_bmd_dedup
opt python scripts/08b_discordance_recordlevel.py --drug imipenem \
    --standard CLSI --method "Broth dilution" --breakpoints CLSI \
    --dedup --tag _clsi_bmd_dedup

say "5. unstratified analyses (supplementary)"
opt python scripts/08b_discordance_recordlevel.py --drug meropenem \
    --breakpoints CLSI --tag _unstratified
opt python scripts/08b_discordance_recordlevel.py --drug imipenem \
    --breakpoints CLSI --tag _unstratified

say "6. sensitivity analyses"
must python scripts/09_discordance_sensitivity.py --drug meropenem --tag _clsi_bmd
opt  python scripts/09_discordance_sensitivity.py --drug imipenem  --tag _clsi_bmd

say "7. extension analyses (determinants, predictive value, projection)"
opt python scripts/11_extensions.py --drug meropenem --tag _clsi_bmd

say "8. figures — PNG, PDF and SVG; no captions drawn inside"
# Figure 1 — study flow
must python scripts/12_flow_diagram.py
# Figure 1 again, typeset from the layout template: same numbers, but a
# single text node per label so the figure can be opened and corrected
must python scripts/12b_flow_diagram_typeset.py
# Figure 2 — primary discordance analysis
must python scripts/10_figures.py --drug meropenem --tag _clsi_bmd
# Supplementary figures S1-S3
opt  python scripts/10_figures.py --drug imipenem  --tag _clsi_bmd
opt  python scripts/10_figures.py --drug meropenem --tag _clsi_bmd_eucast
opt  python scripts/10_figures.py --drug meropenem --tag _unstratified
opt  python scripts/10_figures.py --drug meropenem --tag _clsi_bmd_dedup

say "9. results manifest — the single source of truth"
must python scripts/14_results_manifest.py

say "10. Figures 3 and 4 (read from the manifest and record-level tables)"
# Figure 3 — robustness across analysis variants
must python scripts/15_figure3_robustness.py
# Figure 4 — resistance without an acquired carbapenemase
must python scripts/16_figure4_false_resistant.py

say "11. supplementary tables"
must python scripts/13_supplementary.py

# ------------------------------------------------- predictive rules
# These answer the editorial request for usable rules rather than associations.
# 17 derives them with clone-aware cross-validation; 18 interrogates the
# discordant isolates themselves; 19 applies the rules, unchanged, to an
# independent NCBI collection. 19 is optional only because it needs the
# downloaded validation metadata, not because it is secondary.
say "12. predictive rule derivation"
must python scripts/17_rule_derivation.py --drug meropenem
opt  python scripts/17_rule_derivation.py --drug imipenem --tag _imipenem

say "13. discordance investigation"
# Pass --qc data/interim/assembly_qc_metrics.tsv once 00b_qc_assembly.py has
# run, so that assembly quality is tested as a competing explanation.
if [ -f data/interim/assembly_qc_metrics.tsv ]; then
  must python scripts/18_discordance_investigation.py --drug meropenem \
      --qc data/interim/assembly_qc_metrics.tsv
else
  echo "   NOTE: no assembly QC metrics; the technical explanation for"
  echo "   discordance will be reported as untested."
  must python scripts/18_discordance_investigation.py --drug meropenem
fi

say "14. external validation"
# A missing external table used to increment WARN only. WARN is not what the
# exit gate tests, so a run that never validated against NCBI Pathogen
# Detection could still finish with status 0 and be described as a complete
# rebuild. A required stage that did not run is recorded as a skipped step,
# whether it failed or its input was absent.
if [ -f data/external/ncbi_kp_ast.tsv ]; then
  opt python scripts/19_external_validation.py --drug meropenem
else
  echo "   SKIPPED: external input data/external/ncbi_kp_ast.tsv is absent."
  echo "   Run scripts/19b_fetch_ncbi_pathogen.sh first, or set"
  echo "   ALLOW_SKIPPED=1 to declare this a deliberate downstream-only run."
  WARN=$((WARN+1))
  SKIPPED_STEPS+=("external validation (data/external/ncbi_kp_ast.tsv absent)")
fi

say "15. provenance, integrity and release reporting"
# These four were written after the download repair and reported results of
# their own, but they sat outside run_all.sh, so no single logged run covered
# every number the manuscript cites. They run here.
# Effective download integrity for the files that were actually genotyped.
if [ -f results/qc_concordance/fasta_integrity_report.tsv ]; then
  must python scripts/00g_effective_integrity_report.py
else
  echo "   SKIPPED: results/qc_concordance/fasta_integrity_report.tsv is absent."
  echo "   Run scripts/00f_check_download_integrity.py against the assembly"
  echo "   directories first."
  WARN=$((WARN+1))
  SKIPPED_STEPS+=("effective download integrity (integrity report absent)")
fi
# Where every genome left the study, at genome-ID level.
must python scripts/25_genome_ledger.py
# Does the porin rule's gain survive the clonality of the derivation set?
# The grouped CV stops a rule being scored on its own clone but puts no
# interval on the gain, which is the claim the paper makes.
must python scripts/27_derivation_gain_bootstrap.py --drug meropenem
# Species and integrity sensitivity (Supplementary Table S18).
must python scripts/24_species_and_integrity_sensitivity.py
# Plain-text summary of the external validation.
if [ -f results/tables/external_performance.tsv ]; then
  must python scripts/23_external_summary.py
else
  echo "   SKIPPED: results/tables/external_performance.tsv is absent."
  echo "   It is written by scripts/19_external_validation.py."
  WARN=$((WARN+1))
  SKIPPED_STEPS+=("external validation summary (external results absent)")
fi

# The audit trail for this run: environment, script versions, input and output
# checksums, and the steps that did not run. Written last so it records the
# summary this run produced.
say "DONE"
write_summary
# Plain copies, not symlinks: a symlink does not survive a ZIP extracted on
# Windows, and a pointer is only worth having once there is a finished run to
# point at.
cp -f "$LOG" results/logs/latest.log 2>/dev/null || true
cp -f "$SUMMARY" results/logs/latest_summary.txt 2>/dev/null || true
must python scripts/26_release_manifest.py
if [ "$WARN" -gt 0 ]; then
  echo "   $WARN optional step(s) were skipped — see the log above and the"
  echo "   'NOT YET GENERATED' list from the results manifest."
fi
cat <<'MSG'
Reconcile the manuscript against results/manuscript_results.json ONLY.
Any number in the abstract, Results, Discussion or tables that does not
appear there has not been regenerated and must not be cited.

Logs:
  results/logs/latest.log            full transcript of this run
  results/logs/latest_summary.txt    short summary — send THIS when asking
                                     about a run

Generated:
  results/manuscript_results.json   headline values (SINGLE SOURCE OF TRUTH)
  results/tables/recordlevel_*      analysis tables
  results/supplementary/            supplementary tables + xlsx

  results/figures/  (PNG, PDF and SVG; no captions drawn inside)
    Figure1_study_flow                       Figure 1
    Figure2_main_meropenem_CLSI              Figure 2  (primary analysis)
    Figure3_robustness_meropenem             Figure 3  (robustness)
    Figure4_carbapenemase_negative_meropenem Figure 4  (no-carbapenemase direction)
    FigureS1_imipenem_CLSI                   Figure S1
    FigureS2_meropenem_unstratified          Figure S2
    FigureS3_meropenem_deduplicated          Figure S3
    FigureS4_meropenem_EUCAST                Figure S4

  File names now state their role in the manuscript. Previously every
  six-panel figure was prefixed "figure2_" because that names the template
  rather than the manuscript figure, which made them impossible to tell apart.
MSG

# --------------------------------------------------------------- skip gate
finish_gate
