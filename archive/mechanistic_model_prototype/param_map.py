"""Stage 2 — map each hotspot class to a mechanistic kinetic parameter.

INPUT : data/interim/hotspots.csv, data/raw/enzyme_kinetics/kinetic_constants_curated.csv
OUTPUT: data/processed/genotype_parameters.csv

Mapping (fix values from INDEPENDENT sources; do not tune to time-kill curves):
  porin loss (ompK36 GD/TD, ompK35 stop)  -> permeability P  (down)
  carbapenemase allele + copy number      -> enzyme level E, kcat, Km
  efflux up (ramR / acrR loss of function)-> efflux term
  CZA-evading KPC variants                -> avibactam acylation / off-rate
"""
raise NotImplementedError("Fill in once enzyme constants are curated. See data/DATA_SOURCES.md")
