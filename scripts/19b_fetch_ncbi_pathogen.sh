#!/usr/bin/env bash
# 19b_fetch_ncbi_pathogen.sh — retrieve the external validation metadata.
#
# NCBI Pathogen Detection publishes a per-isolate metadata table for each
# organism group. For Klebsiella it carries submitted AST results, AMRFinderPlus
# genotypes (including ompK35/ompK36 protein changes) and precomputed assembly
# statistics, which is everything 19_external_validation.py needs.
#
# The full table is ~190 MB, so it is filtered during download rather than
# stored: only K. pneumoniae rows that carry any AST result are kept.
#
# The snapshot identifier is recorded in the output directory because NCBI
# replaces "latest" continuously. Pin it in the manuscript; without it the
# validation set cannot be reconstructed.
#
# Usage:
#   bash scripts/19b_fetch_ncbi_pathogen.sh
#   SNAPSHOT=PDG000000012.2539 bash scripts/19b_fetch_ncbi_pathogen.sh
set -uo pipefail

BASE="https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella"
OUT="${OUT:-data/external}"
mkdir -p "$OUT"

# ------------------------------------------------- resolve the snapshot
# Pinned to the snapshot used in the manuscript so the published one-command
# recipe reproduces the exact 986-isolate cohort. NCBI Pathogen Detection
# accessions are not archival and the directory listing changes, so set
# SNAPSHOT=latest explicitly to fetch a newer release instead.
SNAPSHOT="${SNAPSHOT:-PDG000000012.2548}"
if [ "$SNAPSHOT" = "latest" ]; then
  echo "resolving latest Klebsiella snapshot"
  SNAPSHOT=$(curl -s "$BASE/" \
    | grep -oE 'PDG[0-9]+\.[0-9]+/' | tr -d '/' \
    | sort -t. -k2 -n | tail -1)
fi
[ -n "$SNAPSHOT" ] || { echo "Could not determine a snapshot id"; exit 1; }
echo "snapshot : $SNAPSHOT"

URL="$BASE/$SNAPSHOT/Metadata/$SNAPSHOT.metadata.tsv"
TSV="$OUT/ncbi_kp_ast.tsv"

# ------------------------------------------------- download and filter
# Field positions are resolved from the header rather than hard-coded, because
# NCBI has added columns to this table before now and a fixed index would
# silently filter on the wrong field.
echo "downloading and filtering (this streams ~190 MB)"
curl -s "$URL" | awk -F'\t' -v OFS='\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "scientific_name") sn = i
      if ($i == "AST_phenotypes")  ast = i
    }
    if (!sn || !ast) {
      print "ERROR: scientific_name or AST_phenotypes column not found" > "/dev/stderr"
      exit 1
    }
    print
    next
  }
  $sn ~ /Klebsiella pneumoniae/ && $ast != "" && $ast != "null" { print }
' > "$TSV" || { echo "download or filter failed"; exit 1; }

rows=$(( $(wc -l < "$TSV") - 1 ))
echo "$SNAPSHOT" > "$OUT/SNAPSHOT.txt"
date -u +"%Y-%m-%dT%H:%M:%SZ" >> "$OUT/SNAPSHOT.txt"

echo
echo "isolates with any AST : $rows"
echo "table                 : $TSV"
echo "snapshot recorded in  : $OUT/SNAPSHOT.txt"
echo
echo "Next:"
echo "  python scripts/19a_fetch_derivation_accessions.py   # overlap exclusion"
echo "  python scripts/19_external_validation.py --drug meropenem"
