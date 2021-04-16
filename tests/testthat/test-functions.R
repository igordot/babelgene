
test_that("orthologs()", {
  orthologs_df <- orthologs()
  expect_s3_class(orthologs_df, "data.frame")
  expect_lt(nrow(orthologs_df), 100)
})
