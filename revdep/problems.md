# orthogene

<details>

* Version: 1.2.0
* GitHub: https://github.com/neurogenomics/orthogene
* Source code: https://github.com/cran/orthogene
* Date/Publication: 2022-04-26
* Number of recursive dependencies: 210

Run `revdep_details(, "orthogene")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      [ FAIL 4 | WARN 0 | SKIP 0 | PASS 141 ]
      
      ══ Failed tests ════════════════════════════════════════════════════════════════
      ── Failure (test-map_orthologs_babelgene.R:12:5): map_orthologs_babelgene works ──
      nrow(gene_map_b1) is not more than 13100. Difference: -34
      ── Failure (test-map_orthologs_babelgene.R:29:5): map_orthologs_babelgene works ──
      nrow(gene_map_b3) is not more than 15900. Difference: -27
      ── Failure (test-map_orthologs_babelgene.R:42:5): map_orthologs_babelgene works ──
      nrow(gene_map1) is not more than 29700. Difference: -49
      ── Failure (test-map_orthologs_babelgene.R:60:5): map_orthologs_babelgene works ──
      nrow(gene_map2) is not more than 29700. Difference: -49
      
      [ FAIL 4 | WARN 0 | SKIP 0 | PASS 141 ]
      Error: Test failures
      Execution halted
    ```

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘DelayedMatrixStats’
      All declared Imports should be used.
    There are ::: calls to the package's namespace in its code. A package
      almost never needs to use ::: for its own objects:
      ‘aggregate_rows’
    ```

