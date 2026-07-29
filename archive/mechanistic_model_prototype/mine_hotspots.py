"""Stage 1 — mine hotspot mutations from public genomes (population-structure aware).

INPUT : data/raw/genomes/ (assemblies), data/raw/phenotypes/ (MIC / S-I-R)
OUTPUT: data/interim/hotspots.csv  (locus, variant, effect_size, homoplasy_count)

Steps (implement):
  1. Kleborate-type all assemblies -> ST, AMR determinants, virulence loci.
  2. Build unitig / SNP matrix.
  3. pyseer with a phylogeny/kinship matrix (LMM) -> association, CORRECTING
     for population structure. treeWAS as a cross-check.
  4. Homoplasy / independent-recurrence count across the tree to separate true
     convergent hotspots from lineage markers.
  5. Keep mechanistically interpretable loci: ompK35/ompK36, blaKPC/OXA-48/NDM
     alleles + copy number, ramR/acrR (efflux).

DO NOT skip step 3's structure correction — it is the difference between a real
hotspot and a clonal artifact (and between acceptance and desk-reject).
"""
raise NotImplementedError("Fill in once data/raw is populated. See data/DATA_SOURCES.md")
