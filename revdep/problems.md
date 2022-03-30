# orthogene

<details>

* Version: 1.0.0
* GitHub: https://github.com/neurogenomics/orthogene
* Source code: https://github.com/cran/orthogene
* Date/Publication: 2021-10-26
* Number of recursive dependencies: 167

Run `revdep_details(, "orthogene")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
         3,849 / 4,492 (86%)
      Total genes remaining after convert_orthologs :
         643 / 4,492 (14%)
      Finisheddmelanogasterin0.232minutes.
      Saving benchmarking results ==> /var/folders/lj/80zz6n91631_lb6h3l2m3hq4f503zw/T//Rtmpmg1Vml/file37766fe58c9d.csv
      WARNING: Species ' monkeytypo ' not found in taxa dict.
      [ FAIL 1 | WARN 0 | SKIP 0 | PASS 100 ]
      
      ══ Failed tests ════════════════════════════════════════════════════════════════
      ── Failure (test-map_orthologs_babelgene.R:18:5): map_orthologs_babelgene works ──
      nrow(gene_map2) is not more than 13000. Difference: -1.3e+04
      
      [ FAIL 1 | WARN 0 | SKIP 0 | PASS 100 ]
      Error: Test failures
      Execution halted
    ```

