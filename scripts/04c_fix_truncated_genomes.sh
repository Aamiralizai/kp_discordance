#!/usr/bin/env bash
# 04c_fix_truncated_genomes.sh — repair the truncated genome downloads and re-type them.
#
# BACKGROUND
#   About 47% of data/raw/genomes/*.fna (3,453 of 7,360) are truncated FTPS
#   downloads: the transfer stopped mid-sequence and wget -q still reported success.
#   Kleborate typed those truncated files. See results/qc_concordance/
#   fasta_integrity_report.tsv for the per-genome evidence.
#
# WHAT THIS DOES (nothing existing is modified or deleted)
#   1. Verified re-download of the truncated genomes from the BV-BRC HTTPS API into
#      data/raw/genomes_redownload/. A genome is kept only if its length and contig
#      count equal the BV-BRC record and the file is complete.
#   2. Kleborate on the re-downloaded genomes, with the SAME call as the original
#      run (kleborate -a ... -p kpsc -r, chunks of 25, one shard per core), into
#      data/interim/kleborate_redownload/. A chunk that fails (e.g. segfault) is
#      retried genome by genome so one bad genome cannot drop 24 others.
#   3. Merges the main table to data/processed/kleborate_redownload_output.tsv.
#   4. Checks that Kleborate's total_size for every re-typed genome equals the
#      BV-BRC genome length, which proves Kleborate saw the complete genome.
#
# RUN (from the project root, in the same conda env used for the original typing)
#   conda activate kleborate
#   bash scripts/04c_fix_truncated_genomes.sh
#
# Options (environment variables):
#   SHARDS=30      parallel Kleborate workers (default: cores - 1)
#   CHUNK=25       genomes per kleborate call (default 25, as originally)
#   DL_WORKERS=2   parallel download connections (use 1 if downloads stall)
#   DL_MODE=single single | bulk  (bulk = 8 genomes per request; faster on a good link)
#   SKIP_DOWNLOAD=1  only run Kleborate on what is already downloaded
#
# Re-running is safe: verified downloads and typed genomes are skipped.
set -uo pipefail

IDS="results/qc_concordance/truncated_ids.txt"
NEW="data/raw/genomes_redownload"
OUT="data/interim/kleborate_redownload"
DONE="$OUT/done.txt"
MERGED="data/processed/kleborate_redownload_output.tsv"

: "${SHARDS:=$(( $(nproc) > 1 ? $(nproc) - 1 : 1 ))}"
: "${CHUNK:=25}"
: "${DL_WORKERS:=2}"
: "${DL_MODE:=single}"
: "${SKIP_DOWNLOAD:=0}"

[ -f "$IDS" ] || { echo "ERROR: $IDS not found. Run from the project root."; exit 1; }
command -v kleborate >/dev/null || { echo "ERROR: kleborate not found. Run: conda activate kleborate"; exit 1; }
PY=$(command -v python3 || command -v python)
mkdir -p "$NEW" "$OUT" "$(dirname "$MERGED")"
touch "$DONE"
# recover progress from an interrupted run (Ctrl+C before the shards finished)
cat "$DONE".part_* 2>/dev/null >> "$DONE"; rm -f "$DONE".part_*
LOG="$OUT/fix_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
echo "== 04c fix run started $(date) | kleborate: $(kleborate --version 2>&1 | head -1) | shards=$SHARDS chunk=$CHUNK"

# ---------------------------------------------------------------- 1. download
if [ "$SKIP_DOWNLOAD" != "1" ]; then
  echo "== step 1: verified re-download ($(wc -l < "$IDS") genomes requested) =="
  "$PY" scripts/01d_redownload_verified.py "$IDS" "$NEW" \
        --workers "$DL_WORKERS" --mode "$DL_MODE" --rounds 10
  echo "   still missing after download: $(grep -c . "$NEW/failed_ids.txt" 2>/dev/null || echo 0)"
fi

# ---------------------------------------------------------------- 2. kleborate
echo "== step 2: Kleborate on verified genomes =="
find "$NEW" -maxdepth 1 -name '*.fna' -printf '%f\n' | sed 's/\.fna$//' | sort -u > "$OUT/_all.txt"
sort -u "$DONE" -o "$DONE"
comm -23 "$OUT/_all.txt" "$DONE" > "$OUT/_todo.txt"
todo=$(wc -l < "$OUT/_todo.txt")
echo "   verified on disk: $(wc -l < "$OUT/_all.txt") | already typed: $(wc -l < "$DONE") | to type: $todo"

if [ "$todo" -gt 0 ]; then
  echo "   started $(date)"
  rm -f "$OUT"/_shard_*
  split -n r/"$SHARDS" -d "$OUT/_todo.txt" "$OUT/_shard_"
  pids=()
  for sf in "$OUT"/_shard_*; do
    n=$(basename "$sf" | sed 's/_shard_//')
    (
      sdir="$OUT/shard_$n"; mkdir -p "$sdir"
      split -l "$CHUNK" -d "$sf" "$sdir/_chunk_"
      for cf in "$sdir"/_chunk_*; do
        files=(); while read -r id; do [ -n "$id" ] && files+=("$NEW/$id.fna"); done < "$cf"
        [ ${#files[@]} -eq 0 ] && { rm -f "$cf"; continue; }
        if kleborate -a "${files[@]}" -o "$sdir" -p kpsc -r >>"$sdir/kleborate.log" 2>&1; then
          cat "$cf" >> "$DONE.part_$n"
        else
          # isolate the failure: retry each genome on its own
          while read -r id; do
            [ -z "$id" ] && continue
            if kleborate -a "$NEW/$id.fna" -o "$sdir" -p kpsc -r >>"$sdir/kleborate.log" 2>&1; then
              echo "$id" >> "$DONE.part_$n"
            else
              echo "$id" >> "$sdir/failed_genomes.txt"
            fi
          done < "$cf"
        fi
        rm -f "$cf"
      done
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p"; done
  cat "$DONE".part_* 2>/dev/null >> "$DONE"; rm -f "$DONE".part_*
  sort -u "$DONE" -o "$DONE"
  echo "   finished $(date)"
fi

# ---------------------------------------------------------------- 3. merge
echo "== step 3: merging main Kleborate table =="
"$PY" - "$OUT" "$MERGED" <<'PY'
import sys, glob, csv
out_dir, merged = sys.argv[1], sys.argv[2]
rows, header, seen = [], None, set()
for f in sorted(glob.glob(f"{out_dir}/shard_*/klebsiella_pneumo_complex_output.txt")):
    with open(f) as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r, None)
        if h is None: continue
        if header is None: header = h
        if h != header:
            sys.exit(f"ERROR: column mismatch in {f}")
        for row in r:
            if row and row[0] not in seen:
                seen.add(row[0]); rows.append(row)
if header is None:
    sys.exit("ERROR: no klebsiella_pneumo_complex_output.txt found in shards")
with open(merged, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(header); w.writerows(rows)
print(f"   {len(rows)} genomes -> {merged}")
PY

# ---------------------------------------------------------------- 4. verify
echo "== step 4: did Kleborate see complete genomes? =="
"$PY" - "$NEW/manifest.tsv" "$MERGED" <<'PY'
import sys, csv
man, merged = sys.argv[1], sys.argv[2]
exp = {}
for r in csv.DictReader(open(man), delimiter="\t"):
    if r["status"] == "ok":
        exp[r["genome_id"]] = int(r["bvbrc_length"])
ok = bad = 0; badlist = []
for r in csv.DictReader(open(merged), delimiter="\t"):
    g = r["strain"]; ts = r.get("general__contig_stats__total_size", "")
    if g in exp and ts and int(float(ts)) == exp[g]: ok += 1
    else: bad += 1; badlist.append(g)
print(f"   total_size == BV-BRC genome_length: {ok}   mismatched or unverified: {bad}")
if badlist: print("   first mismatches:", ", ".join(badlist[:10]))
PY

echo "== summary =="
echo "   requested       : $(wc -l < "$IDS")"
echo "   verified on disk: $(find "$NEW" -maxdepth 1 -name '*.fna' | wc -l)"
echo "   typed           : $(wc -l < "$DONE")"
echo "   failed genomes  : $(cat "$OUT"/shard_*/failed_genomes.txt 2>/dev/null | wc -l)"
echo "   log             : $LOG"
echo "Send back: $MERGED, $NEW/manifest.tsv, $NEW/failed_ids.txt and $LOG"
rm -f "$OUT"/_all.txt "$OUT"/_todo.txt "$OUT"/_shard_*
