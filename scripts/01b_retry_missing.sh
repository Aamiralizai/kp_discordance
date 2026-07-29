#!/usr/bin/env bash
# 01b_retry_missing.sh — download ONLY the assemblies that are missing or empty.
#
# Safe to run repeatedly. Each run:
#   1. removes zero-byte / truncated .fna files (failed downloads leave these)
#   2. diffs genome_list.txt against what's actually on disk
#   3. re-downloads only the gap, in parallel, with retries
#   4. writes any still-failing IDs to data/interim/failed_ids.txt
#
# Usage (from project root, inside WSL):
#   bash scripts/01b_retry_missing.sh              # default 6 parallel, 3 tries
#   JOBS=10 TRIES=5 bash scripts/01b_retry_missing.sh
#   MIN_BYTES=10000 bash scripts/01b_retry_missing.sh   # stricter size floor
set -uo pipefail

PHEN="data/raw/phenotypes"
GEN="data/raw/genomes"
INT="data/interim"
LIST="$PHEN/genome_list.txt"
FTP="ftps://ftp.bv-brc.org"

: "${JOBS:=6}"        # parallel downloads
: "${TRIES:=3}"       # wget attempts per file
: "${MIN_BYTES:=1000}" # anything smaller is treated as a failed download

mkdir -p "$GEN" "$INT"
[ -s "$LIST" ] || { echo "ERROR: $LIST not found. Run 01_download_bvbrc.sh first."; exit 1; }

echo "== step 1: clearing failed/empty downloads (< ${MIN_BYTES} bytes) =="
before=$(find "$GEN" -name '*.fna' | wc -l)
find "$GEN" -name '*.fna' -size -"${MIN_BYTES}"c -delete
after=$(find "$GEN" -name '*.fna' | wc -l)
echo "   removed $((before - after)) bad files; $after valid assemblies remain"

echo "== step 2: computing what's still missing =="
# strip CR (Windows line endings) and blanks from the ID list
tr -d '\r' < "$LIST" | sed '/^[[:space:]]*$/d' | sort -u > "$INT/_want.txt"
find "$GEN" -name '*.fna' -printf '%f\n' 2>/dev/null | sed 's/\.fna$//' | sort -u > "$INT/_have.txt"
comm -23 "$INT/_want.txt" "$INT/_have.txt" > "$INT/_missing.txt"
want=$(wc -l < "$INT/_want.txt"); have=$(wc -l < "$INT/_have.txt"); miss=$(wc -l < "$INT/_missing.txt")
echo "   wanted: $want | have: $have | missing: $miss"
[ "$miss" -eq 0 ] && { echo "Nothing to do — all assemblies present."; exit 0; }

echo "== step 3: downloading $miss missing assemblies (${JOBS} parallel, ${TRIES} tries) =="
export GEN FTP TRIES MIN_BYTES
cat "$INT/_missing.txt" | xargs -P "$JOBS" -I{} bash -c '
  id="{}"; out="$GEN/$id.fna"
  wget --ftp-user=anonymous --ftp-password=guest --secure-protocol=auto \
       -q -t "$TRIES" --waitretry=2 --timeout=30 -O "$out" \
       "$FTP/genomes/$id/$id.fna" 2>/dev/null
  # discard empty/short results so the next run retries them
  if [ ! -s "$out" ] || [ "$(stat -c%s "$out" 2>/dev/null || echo 0)" -lt "$MIN_BYTES" ]; then
    rm -f "$out"; echo "FAIL $id"
  fi
'

echo "== step 4: final tally =="
find "$GEN" -name '*.fna' -size -"${MIN_BYTES}"c -delete
find "$GEN" -name '*.fna' -printf '%f\n' 2>/dev/null | sed 's/\.fna$//' | sort -u > "$INT/_have.txt"
comm -23 "$INT/_want.txt" "$INT/_have.txt" > "$INT/failed_ids.txt"
ok=$(wc -l < "$INT/_have.txt"); bad=$(wc -l < "$INT/failed_ids.txt")
echo "   downloaded OK : $ok"
echo "   still failing : $bad  (listed in $INT/failed_ids.txt)"
rm -f "$INT/_want.txt" "$INT/_have.txt" "$INT/_missing.txt"

if [ "$bad" -gt 0 ]; then
  echo
  echo "Re-run this script to retry. IDs that fail every time usually have no"
  echo ".fna deposited at that BV-BRC path (not a network problem) — those are"
  echo "safe to drop; just exclude them from the phenotype table before modelling."
fi
