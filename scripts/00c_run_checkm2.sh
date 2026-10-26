#!/usr/bin/env bash
# 00c_run_checkm2.sh — completeness and contamination for every assembly.
#
# CheckM2 supplies the two QualiBact metrics that cannot be read off a FASTA.
# Contamination is the one that matters most here: a contaminated assembly can
# contribute a carbapenemase gene belonging to another organism, producing a
# genotype-positive, phenotype-susceptible isolate that is an artefact of the
# assembly rather than a biological discordance. That is the exact cell of the
# contingency table this study reports on, so it cannot be left unmeasured.
#
# CheckM2 runs in its own environment; its dependencies conflict with Kleborate's.
#
#   conda create -n checkm2 -c bioconda -c conda-forge checkm2=1.1.0
#   conda activate checkm2
#   checkm2 database --download            # ~2.9 GB, once
#   bash scripts/00c_run_checkm2.sh
#
# Usage:
#   bash scripts/00c_run_checkm2.sh                 # all assemblies
#   THREADS=32 bash scripts/00c_run_checkm2.sh      # set thread count
#   CHUNK=500  bash scripts/00c_run_checkm2.sh      # assemblies per invocation
#
# Output:
#   data/interim/checkm2/quality_report.tsv   merged, for 00b_qc_assembly.py
set -uo pipefail

GEN="${GEN:-data/raw/genomes}"
OUT="${OUT:-data/interim/checkm2}"
: "${THREADS:=$(nproc 2>/dev/null || echo 4)}"
: "${CHUNK:=500}"
: "${EXT:=fna}"

command -v checkm2 >/dev/null || {
  echo "ERROR: checkm2 not found. Run: conda activate checkm2"; exit 1; }

mkdir -p "$OUT/shards"
MERGED="$OUT/quality_report.tsv"

mapfile -t ALL < <(find "$GEN" -maxdepth 1 -name "*.${EXT}" | sort)
total=${#ALL[@]}
[ "$total" -gt 0 ] || { echo "No *.${EXT} assemblies under $GEN"; exit 1; }
echo "assemblies : $total"
echo "threads    : $THREADS"
echo "chunk size : $CHUNK"

shard=0
for ((i = 0; i < total; i += CHUNK)); do
  shard=$((shard + 1))
  sdir="$OUT/shards/shard_$(printf '%04d' "$shard")"

  # resumable: a shard with a report is left alone
  if [ -s "$sdir/quality_report.tsv" ]; then
    echo "shard $shard — already done, skipping"
    continue
  fi

  batch=("${ALL[@]:i:CHUNK}")
  echo "shard $shard — ${#batch[@]} assemblies"
  rm -rf "$sdir"
  if ! checkm2 predict --threads "$THREADS" --force \
        --input "${batch[@]}" --output-directory "$sdir" \
        >"$OUT/shards/shard_$(printf '%04d' "$shard").log" 2>&1; then
    echo "   !! shard $shard FAILED — see its log" >&2
  fi
done

# ---------------------------------------------------------------- merge
echo
echo "merging shard reports"
first=1
: > "$MERGED"
for f in "$OUT"/shards/shard_*/quality_report.tsv; do
  [ -s "$f" ] || continue
  if [ "$first" = 1 ]; then
    cat "$f" >> "$MERGED"; first=0
  else
    tail -n +2 "$f" >> "$MERGED"
  fi
done

rows=$(( $(wc -l < "$MERGED") - 1 ))
echo "merged     : $rows genomes -> $MERGED"
if [ "$rows" -lt "$total" ]; then
  echo
  echo "WARNING: $((total - rows)) assemblies have no CheckM2 result."
  echo "Those genomes will be recorded as not_evaluated for completeness and"
  echo "contamination, and the QC verdict will be reported as PARTIAL."
fi
echo
echo "Next: python scripts/00b_qc_assembly.py --checkm2 $MERGED --apply"
