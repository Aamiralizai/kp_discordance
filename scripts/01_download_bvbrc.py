#!/usr/bin/env python
"""01_download_bvbrc.py  —  Windows-friendly BV-BRC downloader.

Cross-platform (Windows / macOS / Linux). Does the same job as
scripts/01_download_bvbrc.sh but in Python, using `curl` for the FTPS
transfers. curl.exe ships with Windows 10 (1803+) and Windows 11, so no
extra install is needed on a modern Windows box.

Run from the project root:
    python scripts/01_download_bvbrc.py
    python scripts/01_download_bvbrc.py --max-genomes 200   # cap for a first slice

Outputs:
    data/raw/phenotypes/PATRIC_genomes_AMR.txt   full lab AMR table (large)
    data/raw/phenotypes/kp_amr_anchor.tsv        filtered slice
    data/raw/phenotypes/genome_list.txt          unique genome IDs
    data/raw/genomes/<id>.fna                     assemblies
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

FTP_HOST = "ftp.bv-brc.org"
# BV-BRC uses explicit FTPS. curl forces it with --ssl-reqd on an ftp:// URL.
# -k skips certificate verification (BV-BRC docs note verify can fail); drop it
# if your environment has the CA chain and you prefer strict verification.
CURL = ["curl", "--ssl-reqd", "-k", "--user", "anonymous:guest", "-sS", "-o"]

SPECIES = "Klebsiella pneumoniae"
# Override with: --drugs meropenem imipenem ciprofloxacin ...
ANCHOR_DRUGS = {"ceftazidime/avibactam", "ceftazidime-avibactam",
                "meropenem", "imipenem"}

PHEN = Path("data/raw/phenotypes")
GEN = Path("data/raw/genomes")


def curl_download(remote_path: str, dest: Path) -> bool:
    url = f"ftp://{FTP_HOST}/{remote_path}"
    try:
        subprocess.run(CURL + [str(dest), url], check=True)
        return dest.exists() and dest.stat().st_size > 0
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  ! curl failed for {remote_path}: {e}", file=sys.stderr)
        return False


def filter_amr(src: Path, out: Path) -> pd.DataFrame:
    """Species + anchor drugs, LAB phenotypes only, must carry R/S or an MIC."""
    df = pd.read_csv(src, sep="\t", dtype=str, low_memory=False).fillna("")
    cols = {c.lower().strip(): c for c in df.columns}
    gn = cols.get("genome_name"); ab = cols.get("antibiotic")
    rp = cols.get("resistant_phenotype"); mv = cols.get("measurement_value")
    lm = cols.get("laboratory_typing_method")
    if not gn or not ab:
        sys.exit("ERROR: expected columns (genome_name, antibiotic) not found.")

    m = df[gn].str.contains(SPECIES, case=False, na=False)
    m &= df[ab].str.lower().str.strip().isin(ANCHOR_DRUGS)
    if lm:                                   # drop computational predictions
        m &= ~df[lm].str.contains("computational", case=False, na=False)
    has_pheno = df[rp].str.strip().ne("") if rp else False
    has_mic = df[mv].str.strip().ne("") if mv else False
    m &= (has_pheno | has_mic)

    kept = df[m].copy()
    kept.to_csv(out, sep="\t", index=False)
    return kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-genomes", type=int, default=None,
                    help="cap number of assemblies to download (first slice)")
    ap.add_argument("--drugs", nargs="*", default=None,
                    help="antibiotics to retrieve (default: the carbapenem set)")
    args = ap.parse_args()
    if args.drugs:
        global ANCHOR_DRUGS
        ANCHOR_DRUGS = {d.strip().lower() for d in args.drugs}
        print(f"retrieving phenotypes for: {sorted(ANCHOR_DRUGS)}")
    PHEN.mkdir(parents=True, exist_ok=True)
    GEN.mkdir(parents=True, exist_ok=True)

    print("[1/4] Downloading master lab-AMR table ...")
    amr = PHEN / "PATRIC_genomes_AMR.txt"
    ok = curl_download("RELEASE_NOTES/PATRIC_genomes_AMR.txt", amr) or \
        curl_download("RELEASE_NOTES/PATRIC_genome_AMR.txt", amr)
    if not ok:
        sys.exit("Could not download AMR table. Check outbound FTPS / firewall, "
                 "or pull it manually from the portal (Route A).")

    print("[2/4] Filtering to species + anchor drugs (lab only) ...")
    kept = filter_amr(amr, PHEN / "kp_amr_anchor.tsv")
    print(f"    kept {len(kept)} records -> {PHEN/'kp_amr_anchor.tsv'}")

    print("[3/4] Extracting unique genome IDs ...")
    gid_col = next((c for c in kept.columns if c.lower().strip() == "genome_id"), None)
    ids = sorted(kept[gid_col].astype(str).unique()) if gid_col else []
    (PHEN / "genome_list.txt").write_text("\n".join(ids))
    print(f"    {len(ids)} unique genomes")

    if args.max_genomes:
        ids = ids[: args.max_genomes]
        print(f"    (capped to first {len(ids)} for this run)")

    print("[4/4] Downloading assemblies (.fna) ...")
    for i, gid in enumerate(ids, 1):
        dest = GEN / f"{gid}.fna"
        if dest.exists() and dest.stat().st_size > 0:
            continue
        print(f"    ({i}/{len(ids)}) {gid}")
        curl_download(f"genomes/{gid}/{gid}.fna", dest)
    print(f"Done. Assemblies -> {GEN}/ ; phenotypes -> {PHEN}/")


if __name__ == "__main__":
    main()
