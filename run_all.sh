#!/usr/bin/env bash
# run_all.sh — clean rebuild of every analysis output.
#
# Empties results/, regenerates everything in dependency order from the
# record-level pipeline, and stops at the first failure. Nothing downstream
# ever reads a stale or legacy file.
#
# PREREQUISITES (not run here — they are slow and rarely change):
#   scripts/01_download_bvbrc.sh      phenotype records + assemblies
#   scripts/00_qc_fasta.sh            assembly QC
#   scripts/04_run_kleborate.sh       genotyping
#   scripts/05_merge_kleborate.py     merge shards
#
# Requires verified breakpoints in config/breakpoints.csv.
#
# Usage:
#   bash run_all.sh              # full rebuild
#   KEEP=1 bash run_all.sh       # keep existing results (no wipe)

set -uo pipefail          # note: not -e; failures are handled per-step
: "${KEEP:=0}"

# ---------------------------------------------------------------- logging
STAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p results/logs
LOG="results/logs/run_${STAMP}.log"
SUMMARY="results/logs/summary_${STAMP}.txt"
ln -sf "run_${STAMP}.log" results/logs/latest.log 2>/dev/null || true
ln -sf "summary_${STAMP}.txt" results/logs/latest_summary.txt 2>/dev/null || true

# everything printed from here on goes to the terminal AND the log
exec > >(tee -a "$LOG") 2>&1

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
  fi
}

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

# ---------------------------------------------------------------- checks
need config/breakpoints.csv
need data/processed/kleborate_klebsiella_pneumo_complex_output.tsv
need data/raw/phenotypes/kp_amr_anchor.tsv

echo "kp-multiscale-amr — clean rebuild"
echo "started $(date)  |  log: $LOG"
echo "python $(python -V 2>&1)  |  env ${CONDA_DEFAULT_ENV:-<none>}"

say "0. regression tests (must pass before anything is generated)"
python tests/test_classification.py

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

say "DONE"
write_summary
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
