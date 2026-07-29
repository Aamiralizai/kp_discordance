#!/usr/bin/env bash
# 04_run_kleborate.sh — type all assemblies with Kleborate 3.x (kpsc preset).
#
# Parallel + resumable:
#   * splits the genome list into SHARDS, one worker per shard, each with its
#     own outdir (avoids write races that -r alone would cause)
#   * records completed genomes in data/interim/kleborate/done.txt
#   * re-running skips anything already typed
#
# Run inside the `kleborate` conda env, from the project root:
#   conda activate kleborate
#   bash scripts/04_run_kleborate.sh
#
# Options:
#   SHARDS=8   bash scripts/04_run_kleborate.sh    # parallel workers (default: nproc-1)
#   LIMIT=200  bash scripts/04_run_kleborate.sh    # only type 200 (timing test)
#   CHUNK=25   bash scripts/04_run_kleborate.sh    # genomes per kleborate call
set -uo pipefail

GEN="data/raw/genomes"
OUT="data/interim/kleborate"
DONE="$OUT/done.txt"
MERGED="data/processed/kleborate_results.tsv"

: "${SHARDS:=$(( $(nproc) > 1 ? $(nproc) - 1 : 1 ))}"
: "${CHUNK:=25}"
: "${LIMIT:=0}"

command -v kleborate >/dev/null || {
  echo "ERROR: kleborate not found. Run: conda activate kleborate"; exit 1; }

mkdir -p "$OUT" "$(dirname "$MERGED")"
touch "$DONE"

echo "== step 0: reconciling done.txt with actual shard outputs =="
# guard against a stale/oversized done.txt causing genomes to be skipped
if ls "$OUT"/shard_*/klebsiella_pneumo_complex_output.txt >/dev/null 2>&1; then
  cat "$OUT"/shard_*/klebsiella_pneumo_complex_output.txt 2>/dev/null \
    | grep -v '^strain' | cut -f1 | sort -u > "$OUT/_actual.txt"
  if [ -s "$OUT/_actual.txt" ]; then
    cp "$OUT/_actual.txt" "$DONE"
    echo "   done.txt set to $(wc -l < "$DONE") genomes actually present in shards"
  fi
fi

echo "== step 1: work list =="
find "$GEN" -name '*.fna' -size +1000c -printf '%f\n' | sed 's/\.fna$//' | sort -u > "$OUT/_all.txt"
sort -u "$DONE" -o "$DONE"
comm -23 "$OUT/_all.txt" "$DONE" > "$OUT/_todo.txt"
[ "$LIMIT" -gt 0 ] && { head -n "$LIMIT" "$OUT/_todo.txt" > "$OUT/_t"; mv "$OUT/_t" "$OUT/_todo.txt"; }
todo=$(wc -l < "$OUT/_todo.txt")
echo "   total: $(wc -l < "$OUT/_all.txt") | already typed: $(wc -l < "$DONE") | to do: $todo"
[ "$todo" -eq 0 ] && { echo "Nothing to type."; }

if [ "$todo" -gt 0 ]; then
  echo "== step 2: typing with $SHARDS parallel workers (chunks of $CHUNK) =="
  echo "   started $(date)"
  # NEVER delete shard_* — they hold results from previous runs.
  # Only clear the transient shard WORK LISTS.
  rm -f "$OUT"/_shard_*
  split -n r/"$SHARDS" -d "$OUT/_todo.txt" "$OUT/_shard_"

  pids=()
  for sf in "$OUT"/_shard_*; do
    n=$(basename "$sf" | sed 's/_shard_//')
    (
      sdir="$OUT/shard_$n"; mkdir -p "$sdir"
      # process this shard in chunks so progress survives interruption
      split -l "$CHUNK" -d "$sf" "$sdir/_chunk_"
      for cf in "$sdir"/_chunk_*; do
        files=(); while read -r id; do
          [ -n "$id" ] && files+=("$GEN/$id.fna")
        done < "$cf"
        [ ${#files[@]} -eq 0 ] && continue
        if kleborate -a "${files[@]}" -o "$sdir" -p kpsc -r >>"$sdir/kleborate.log" 2>&1; then
          # only mark done if the call succeeded
          cat "$cf" >> "$DONE.part_$n"
        else
          echo "chunk failed: $cf" >> "$sdir/failed_chunks.txt"
        fi
        rm -f "$cf"
      done
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p"; done

  cat "$DONE".part_* 2>/dev/null >> "$DONE"
  rm -f "$DONE".part_*
  sort -u "$DONE" -o "$DONE"
  echo "   finished $(date)"
fi

echo "== step 3: merging shard outputs =="
# find the main kpsc result table in each shard (name varies by version)
first=1; : > "$MERGED"
for f in $(find "$OUT" -path '*shard_*' -name '*output.txt' -o -path '*shard_*' -name '*.txt' \
           | grep -v -e kleborate.log -e failed_chunks | sort); do
  case "$(basename "$f")" in
    *output*.txt|*results*.txt) ;;
    *) continue ;;
  esac
  if [ "$first" -eq 1 ]; then cat "$f" >> "$MERGED"; first=0
  else tail -n +2 "$f" >> "$MERGED"; fi
done
rows=$(( $(wc -l < "$MERGED") - 1 ))
echo "   merged rows: $rows -> $MERGED"

echo "== summary =="
echo "   typed genomes : $(wc -l < "$DONE")"
fc=$(cat "$OUT"/shard_*/failed_chunks.txt 2>/dev/null | wc -l)
echo "   failed chunks : $fc"
[ "$fc" -gt 0 ] && echo "   (re-run this script to retry them)"
rm -f "$OUT"/_all.txt "$OUT"/_todo.txt "$OUT"/_shard_*
