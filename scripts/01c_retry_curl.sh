#!/usr/bin/env bash
# 01c_retry_curl.sh — download missing BV-BRC assemblies using curl (not wget).
#
# Why curl: wget can crash (SIGABRT) on repeated FTPS handshakes, and with -O it
# leaves 0-byte files behind on failure. curl -f fails cleanly without creating
# the file, and --retry adds proper backoff.
#
# Fully resumable. Run it as many times as you like.
#
# Usage (project root, inside WSL):
#   bash scripts/01c_retry_curl.sh
#   JOBS=4 bash scripts/01c_retry_curl.sh          # gentler on the server
#   MAX=500 bash scripts/01c_retry_curl.sh         # only fetch 500 this run
set -uo pipefail

PHEN="data/raw/phenotypes"
GEN="data/raw/genomes"
INT="data/interim"
LIST="$PHEN/genome_list.txt"
BASE="ftp://ftp.bv-brc.org/genomes"

: "${JOBS:=4}"          # parallel transfers (keep modest; FTPS is stateful)
: "${TRIES:=4}"         # curl retries per file
: "${MIN_BYTES:=1000}"  # smaller than this = failed download
: "${MAX:=0}"           # 0 = no cap

command -v curl >/dev/null || { echo "ERROR: curl not installed. sudo apt install curl"; exit 1; }
mkdir -p "$GEN" "$INT"
[ -s "$LIST" ] || { echo "ERROR: $LIST missing. Run 01_download_bvbrc.sh first."; exit 1; }

echo "== step 1: removing failed/empty files (< ${MIN_BYTES} bytes) =="
b=$(find "$GEN" -name '*.fna' | wc -l)
find "$GEN" -name '*.fna' -size -"${MIN_BYTES}"c -delete
a=$(find "$GEN" -name '*.fna' | wc -l)
echo "   removed $((b-a)); $a valid assemblies on disk"

echo "== step 2: diffing against genome list =="
tr -d '\r' < "$LIST" | sed '/^[[:space:]]*$/d' | sort -u > "$INT/_want.txt"
find "$GEN" -name '*.fna' -printf '%f\n' 2>/dev/null | sed 's/\.fna$//' | sort -u > "$INT/_have.txt"
comm -23 "$INT/_want.txt" "$INT/_have.txt" > "$INT/_missing.txt"
miss=$(wc -l < "$INT/_missing.txt")
echo "   wanted: $(wc -l < "$INT/_want.txt") | have: $(wc -l < "$INT/_have.txt") | missing: $miss"
[ "$miss" -eq 0 ] && { echo "All assemblies present."; rm -f "$INT"/_*.txt; exit 0; }

if [ "$MAX" -gt 0 ]; then
  head -n "$MAX" "$INT/_missing.txt" > "$INT/_batch.txt"
  echo "   capping this run at $MAX"
else
  cp "$INT/_missing.txt" "$INT/_batch.txt"
fi

echo "== step 3: downloading with curl (${JOBS} parallel, ${TRIES} retries) =="
export GEN BASE TRIES MIN_BYTES
n=$(wc -l < "$INT/_batch.txt")
cat "$INT/_batch.txt" | xargs -P "$JOBS" -I{} bash -c '
  id="{}"; out="$GEN/$id.fna"
  curl --ssl-reqd --user anonymous:guest \
       -f -sS --retry "$TRIES" --retry-delay 2 --retry-connrefused \
       --connect-timeout 30 --max-time 300 \
       -o "$out" "$BASE/$id/$id.fna" 2>/dev/null
  if [ ! -s "$out" ] || [ "$(stat -c%s "$out" 2>/dev/null || echo 0)" -lt "$MIN_BYTES" ]; then
    rm -f "$out"
  fi
' 
# progress is quiet by design; check the tally below

echo "== step 4: tally =="
find "$GEN" -name '*.fna' -size -"${MIN_BYTES}"c -delete
find "$GEN" -name '*.fna' -printf '%f\n' 2>/dev/null | sed 's/\.fna$//' | sort -u > "$INT/_have2.txt"
comm -23 "$INT/_want.txt" "$INT/_have2.txt" > "$INT/failed_ids.txt"
echo "   assemblies on disk : $(wc -l < "$INT/_have2.txt")"
echo "   still missing      : $(wc -l < "$INT/failed_ids.txt")  -> $INT/failed_ids.txt"
rm -f "$INT"/_want.txt "$INT"/_have.txt "$INT"/_have2.txt "$INT"/_missing.txt "$INT"/_batch.txt
echo
echo "Re-run to retry. IDs failing every pass likely have no .fna at that path."
