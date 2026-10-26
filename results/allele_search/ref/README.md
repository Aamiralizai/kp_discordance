# Reference set for the uncatalogued-allele search

`AMRProt.fa` is the AMRFinderPlus protein reference set, release `2026-08-07.1`
as recorded in `version.txt`, fetched from the NCBI pathogen detection
reference gene catalogue FTP area. It is pinned here rather than downloaded at
run time so that the search can be repeated against the same alleles.

```
sha256  41d5ebf4f807c9590f27de7dd23c9f037afb4d8efd613790c54847fb1fc2896b
records 10078
```

Each header is pipe-delimited. Field 4 is the gene symbol, field 5 the family,
field 8 the drug subclass and field 9 the drug class, so the carbapenemase
references are those whose subclass is `CARBAPENEM`. The companion
`ReferenceGeneCatalog.txt` is not kept here. The copy retrieved during this
work was a truncated transfer, and every field the search needs is already in
the FASTA headers, so it was not required.
