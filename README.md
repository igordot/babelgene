# babelgene

[![R build status](https://github.com/igordot/babelgene/workflows/R-CMD-check/badge.svg)](https://github.com/igordot/babelgene/actions)
[![codecov](https://codecov.io/gh/igordot/babelgene/branch/master/graph/badge.svg)](https://codecov.io/gh/igordot/babelgene)

Genomic analysis of model organisms often requires the use of databases based on human data or making comparisons to patient-derived resources.
This requires converting genes between human and non-human equivalents.
The `babelgene` R package provides provides gene orthologs/homologs:

* in an R-friendly tidy/long format with one gene pair per row
* for multiple frequently studied model organisms, such as mouse, rat, fly, and zebrafish
* sourced from multiple databases
* as gene symbols, NCBI Entrez, or Ensembl IDs
* can be queried without an active internet connection

Check the [documentation website](https://igordot.github.io/babelgene) for more information.
