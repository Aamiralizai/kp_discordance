#!/usr/bin/env bash
# 02_download_ncbi.sh — pull K. pneumoniae assemblies from NCBI (diversity/background).
# NOTE: NCBI genome packages do NOT include MIC/AST phenotypes — BV-BRC remains the
# label source. Use this for genome breadth only.
#
# Requires the NCBI `datasets` CLI (in WSL):
#   conda install -n kp -c conda-forge ncbi-datasets-cli -y
#
#   bash scripts/02_download_ncbi.sh
set -euo pipefail

GEN="data/raw/genomes"; TMP="data/interim/ncbi"
mkdir -p "$GEN" "$TMP"

# --- Filters (edit to taste) — keep the download bounded and high quality ---
TAXON="Klebsiella pneumoniae"
SINCE="01/01/2018"                    # date window
LEVELS="complete,chromosome"          # skip fragmented contig/scaffold assemblies

echo "[preview] how many genomes match these filters? (no download yet)"
datasets summary genome taxon "$TAXON" \
  --assembly-level "$LEVELS" --released-since "$SINCE" \
  --exclude-atypical --as-json-lines 2>/dev/null | wc -l || true
read -r -p "Proceed with download? [y/N] " ans; [ "$ans" = "y" ] || { echo "aborted"; exit 0; }

echo "[1/3] Downloading (dehydrated) ..."
datasets download genome taxon "$TAXON" \
  --assembly-level "$LEVELS" --released-since "$SINCE" \
  --exclude-atypical --exclude-multi-isolate \
  --include genome --dehydrated --filename "$TMP/kp_ncbi.zip"

echo "[2/3] Rehydrating (fetches the actual sequences) ..."
unzip -oq "$TMP/kp_ncbi.zip" -d "$TMP/kp_ncbi"
datasets rehydrate --directory "$TMP/kp_ncbi"

echo "[3/3] Flattening .fna into $GEN/ (named by accession) ..."
find "$TMP/kp_ncbi/ncbi_dataset/data" -name "*.fna" | while read -r f; do
  acc=$(basename "$(dirname "$f")")
  [ -e "$GEN/${acc}.fna" ] || cp "$f" "$GEN/${acc}.fna"
done
echo "Done. NCBI assemblies: $(ls "$GEN"/GC*.fna 2>/dev/null | wc -l) files in $GEN/"
