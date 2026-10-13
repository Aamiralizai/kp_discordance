#!/usr/bin/env bash
# 01_download_bvbrc.sh
# Fetch BV-BRC *Klebsiella pneumoniae* LAB AMR phenotypes for the anchor drugs,
# then download the matching genome assemblies.
#
# RUN ON YOUR MACHINE (needs outbound FTPS to ftp.bv-brc.org; institutional
# firewalls must allow explicit FTPS). This is NOT run inside the sandbox.
#
#   bash scripts/01_download_bvbrc.sh
#
# Outputs:
#   data/raw/phenotypes/PATRIC_genomes_AMR.txt   full lab AMR table (big)
#   data/raw/phenotypes/kp_amr_anchor.tsv        filtered slice (species+drugs, lab only)
#   data/raw/phenotypes/genome_list.txt          unique genome IDs
#   data/raw/genomes/<id>.fna                     assemblies
set -euo pipefail

FTP="ftps://ftp.bv-brc.org"
WGET="wget --ftp-user=anonymous --ftp-password=guest --secure-protocol=auto -q"
PHEN="data/raw/phenotypes"; GEN="data/raw/genomes"
mkdir -p "$PHEN" "$GEN"

SPECIES="Klebsiella pneumoniae"
# Drugs to retrieve. Override without editing this file:
#   DRUGS="meropenem|imipenem|ciprofloxacin|gentamicin" bash scripts/01_download_bvbrc.sh
# Widening this is cheap; the cost is that EVERY added drug needs verified
# breakpoints in config/breakpoints.csv before it can be analysed.
: "${DRUGS:=ceftazidime/avibactam|ceftazidime-avibactam|meropenem|imipenem}"

echo "[1/4] Downloading master lab-AMR table ..."
# filename has appeared as both spellings in BV-BRC docs; try both
$WGET -O "$PHEN/PATRIC_genomes_AMR.txt" "$FTP/RELEASE_NOTES/PATRIC_genomes_AMR.txt" \
  || $WGET -O "$PHEN/PATRIC_genomes_AMR.txt" "$FTP/RELEASE_NOTES/PATRIC_genome_AMR.txt"

echo "[2/4] Filtering: $SPECIES + anchor drugs, LAB phenotypes only ..."
awk -v FS='\t' -v OFS='\t' -v species="$SPECIES" -v drugs="$DRUGS" '
NR==1 {
  for (i=1;i<=NF;i++){ key=$i; gsub(/\r/,"",key); col[key]=i }
  gn=col["genome_name"]; ab=col["antibiotic"];
  rp=col["resistant_phenotype"]; mv=col["measurement_value"];
  lm=col["laboratory_typing_method"];
  if(gn==0||ab==0){ print "ERROR: expected columns not found" > "/dev/stderr"; exit 1 }
  print $0; next
}
{
  name=$(gn); drug=tolower($(ab)); meth=(lm?$(lm):"");
  if (index(name, species)==0) next;             # species
  if (drug !~ drugs) next;                        # anchor drugs
  if (meth ~ /[Cc]omputational/) next;            # drop predicted, keep lab-measured
  if ($(rp)=="" && (mv?$(mv):"")=="") next;        # need R/S or an MIC
  print $0
}' "$PHEN/PATRIC_genomes_AMR.txt" > "$PHEN/kp_amr_anchor.tsv"
echo "    kept $(($(wc -l < "$PHEN/kp_amr_anchor.tsv") - 1)) records -> $PHEN/kp_amr_anchor.tsv"

echo "[3/4] Extracting unique genome IDs ..."
awk -v FS='\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="genome_id")g=i; next}{print $g}' \
  "$PHEN/kp_amr_anchor.tsv" | sort -u > "$PHEN/genome_list.txt"
echo "    $(wc -l < "$PHEN/genome_list.txt") unique genomes"

echo "[4/4] Downloading assemblies (.fna) ..."
# MAX_GENOMES caps the download (default: all). Skips files already present;
# removes empty files left by failed fetches so re-runs resume cleanly.
: "${MAX_GENOMES:=0}"
LIST="$PHEN/genome_list.txt"
total=$(wc -l < "$LIST")
[ "$MAX_GENOMES" -gt 0 ] && total=$MAX_GENOMES
i=0
{ [ "$MAX_GENOMES" -gt 0 ] && head -n "$MAX_GENOMES" "$LIST" || cat "$LIST"; } | \
while read -r id; do
  i=$((i+1)); [ -z "$id" ] && continue
  out="$GEN/$id.fna"
  if [ -s "$out" ]; then echo "[$i/$total] skip $id"; continue; fi
  echo "[$i/$total] $id"
  $WGET -O "$out" "$FTP/genomes/$id/$id.fna" || { echo "  ! failed $id"; rm -f "$out"; }
done
echo "Done. $(ls "$GEN"/*.fna 2>/dev/null | wc -l) assemblies in $GEN/"
echo "  (tip: cap the run with  MAX_GENOMES=400 bash scripts/01_download_bvbrc.sh )"
