# Deprecated — do not use

`08_discordance_analysis_INVALID.py` applies testing-method and interpretive-standard
restrictions AFTER aggregating multiple MIC records per genome, and reconciles
multiple records by taking a median with the logical disjunction of censoring
flags. Both are invalidating errors:

* the resulting "CLSI broth microdilution" stratum is not a record-level stratum;
* the median-with-OR rule can construct an interval supported by no original record.

Superseded by `scripts/08b_discordance_recordlevel.py`. Retained only for
provenance. Its outputs (`results/tables/discordance_*.tsv`) must not be used.
