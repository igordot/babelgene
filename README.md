# babelgene

[![CRAN](https://www.r-pkg.org/badges/version/babelgene)](https://cran.r-project.org/package=babelgene)
[![R build status](https://github.com/igordot/babelgene/workflows/R-CMD-check/badge.svg)](https://github.com/igordot/babelgene/actions)
[![codecov](https://codecov.io/gh/igordot/babelgene/branch/main/graph/badge.svg?token=j2n6FRGaZ7)](https://codecov.io/gh/igordot/babelgene)

Genomic analysis of model organisms often requires the use of databases based on human data or making comparisons to patient-derived resources.
This requires converting genes between human and non-human equivalents.
The `babelgene` R package helps to simplify the process.
It provides provides gene orthologs/homologs:

* in an R-friendly tidy/long format with one gene pair per row
* for multiple frequently studied model organisms, such as mouse, rat, fly, and zebrafish
* sourced from multiple databases
* as gene symbols, NCBI Entrez, and Ensembl IDs
* without a need for an active internet connection

Check the [documentation website](https://igordot.github.io/babelgene/) for more information.
