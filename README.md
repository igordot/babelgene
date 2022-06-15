# babelgene

<!-- badges: start -->
[![CRAN](https://www.r-pkg.org/badges/version/babelgene)](https://cran.r-project.org/package=babelgene)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/last-month/babelgene)](https://cran.r-project.org/package=babelgene)
[![R-CMD-check](https://github.com/igordot/babelgene/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/igordot/babelgene/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/igordot/babelgene/branch/main/graph/badge.svg?token=j2n6FRGaZ7)](https://app.codecov.io/gh/igordot/babelgene)
<!-- badges: end -->

Genomic analysis of model organisms frequently requires the use of databases based on human data or making comparisons to patient-derived resources.
This requires harmonization of gene names into the same gene space.
The `babelgene` R package helps to simplify the conversion process.
It provides gene orthologs/homologs:

* for multiple frequently studied model organisms, such as mouse, rat, fly, and zebrafish
* sourced from multiple databases
* as gene symbols, NCBI Entrez, and Ensembl IDs
* without accessing external resources and requiring an active internet connection
* in an R-friendly "[tidy](https://r4ds.had.co.nz/tidy-data.html)" format with one gene pair per row

Check the [documentation website](https://igordot.github.io/babelgene/) for more information.
