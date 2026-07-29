"""Stage 4 — the headline test: mechanistic vs correlational ML on UNSEEN combinations.

INPUT : data/processed/genotype_parameters.csv, phenotype labels
OUTPUT: results/tables/ood_metrics.csv, results/figures/ood_generalization.png

Protocol:
  - Leave-combination-out split: hold out genotypes whose mutation COMBINATION
    is unseen in training (e.g. carbapenemase + porin loss never co-occur in train),
    mirroring lineage co-occurrence structure.
  - Give the mechanistic model and the ML baselines the SAME training labels.
  - Compare log2-MIC MAE on the held-out convergent genotypes.
  - Baselines: Kleborate-features + gradient boosting; unitig + logistic regression.
  - Propagate kinetic-constant and MIC-label uncertainty (bootstrap) — report CIs.

The synthetic version of exactly this test is in experiments/toy_model.py.
"""
raise NotImplementedError("Fill in once real parameters + labels exist.")
