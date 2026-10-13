#!/usr/bin/env bash
# 00_qc_fasta.sh — find and quarantine malformed assemblies.
#
# Kleborate aborts the ENTIRE invocation on an invalid FASTA (e.g. a header
# with no sequence beneath it), so one bad file kills a whole chunk of genomes.
# This moves offenders out of the way so typing can proceed.
#
# Detects:
#   * zero-length sequences (header immediately followed by another header/EOF)
#   * files with no '>' header at all
#   * empty / tiny files
#
# Usage:
#   bash scripts/00_qc_fasta.sh            # report only (dry run)
#   APPLY=1 bash scripts/00_qc_fasta.sh    # actually quarantine
set -uo pipefail

GEN="data/raw/genomes"
QUAR="data/interim/bad_fasta"
REPORT="data/interim/fasta_qc_report.tsv"
: "${APPLY:=0}"
: "${MIN_BYTES:=1000}"

mkdir -p "$(dirname "$REPORT")"
[ "$APPLY" = "1" ] && mkdir -p "$QUAR"

printf 'file\tproblem\n' > "$REPORT"
total=0; bad=0

for f in "$GEN"/*.fna; do
  [ -e "$f" ] || continue
  total=$((total+1))
  problem=""

  sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if [ "$sz" -lt "$MIN_BYTES" ]; then
    problem="too_small(${sz}b)"
  elif ! grep -q '^>' "$f"; then
    problem="no_fasta_header"
  else
    # zero-length sequence: a '>' line whose next non-empty line is also '>' or EOF
    if awk '
      /^>/ {
        if (prev_header && !seen_seq) { bad=1; exit }
        prev_header=1; seen_seq=0; next
      }
      { if (length($0) > 0) seen_seq=1 }
      END { if (prev_header && !seen_seq) bad=1; exit bad ? 0 : 1 }
    ' "$f"; then
      problem="zero_length_sequence"
    fi
  fi

  if [ -n "$problem" ]; then
    bad=$((bad+1))
    printf '%s\t%s\n' "$(basename "$f")" "$problem" >> "$REPORT"
    if [ "$APPLY" = "1" ]; then
      mv "$f" "$QUAR/" 2>/dev/null || true
    fi
  fi
done

echo "scanned : $total assemblies"
echo "bad     : $bad"
echo "report  : $REPORT"
if [ "$bad" -gt 0 ]; then
  echo
  echo "problem breakdown:"
  tail -n +2 "$REPORT" | cut -f2 | sort | uniq -c | sort -rn
fi
if [ "$APPLY" = "1" ]; then
  echo
  echo "Quarantined to $QUAR/  (assemblies removed from $GEN/)"
  echo "Now re-run typing:  bash scripts/04_run_kleborate.sh"
else
  echo
  echo "DRY RUN — nothing moved. To quarantine:"
  echo "  APPLY=1 bash scripts/00_qc_fasta.sh"
fi
