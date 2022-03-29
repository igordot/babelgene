
test_that("orthologs one gene", {
  o1 <- orthologs(genes = "PTPRC", species = "mouse")
  expect_s3_class(o1, "data.frame")
  expect_equal(ncol(o1), 9)
  expect_equal(nrow(o1), 1)
  expect_identical(colnames(o1[, 1:3]), c("human_symbol", "human_entrez", "human_ensembl"))
  expect_identical(colnames(o1[, 5:7]), c("symbol", "entrez", "ensembl"))
})

test_that("orthologs multiple genes", {
  o5 <- orthologs(genes = c("PTPRC", "EPCAM", "CD34", "CD3D", "CD19"), species = "mouse")
  expect_s3_class(o5, "data.frame")
  expect_equal(nrow(o5), 5)
  expect_equal(o5$symbol, c("Cd19", "Cd34", "Cd3d", "Epcam", "Ptprc"))
})

test_that("orthologs multiple matches", {
  bglap <- orthologs(genes = "BGLAP", species = "mouse")
  expect_s3_class(bglap, "data.frame")
  expect_equal(nrow(bglap), 3)
})

test_that("orthologs all mouse genes", {
  mgi_human_genes <- babelgene:::mgi_orthologs_df$human_symbol
  mgi_human <- orthologs(genes = mgi_human_genes, species = "mouse", human = TRUE)
  expect_s3_class(mgi_human, "data.frame")
  expect_gt(nrow(mgi_human), 18000)
  mgi_mouse_genes <- babelgene:::mgi_orthologs_df$mouse_symbol
  mgi_mouse <- orthologs(genes = mgi_mouse_genes, species = "mouse", human = FALSE)
  expect_s3_class(mgi_mouse, "data.frame")
  expect_gt(nrow(mgi_mouse), 18000)
})

test_that("orthologs identifier types", {
  ptprc_sym <- orthologs(genes = "PTPRC", species = "mouse")
  ptprc_ent <- orthologs(genes = 5788, species = "mouse")
  ptprc_ens <- orthologs(genes = "ENSG00000081237", species = "mouse")
  ptprc_ens_m <- orthologs(genes = "ENSMUSG00000026395", species = "mouse", human = FALSE)
  expect_identical(ptprc_sym, ptprc_ent)
  expect_identical(ptprc_sym, ptprc_ens)
  expect_identical(ptprc_sym, ptprc_ens_m)
})

test_that("orthologs unusual symbols", {
  # gene symbol that can also be Entrez ID
  int_312 <- orthologs(genes = 312, species = "fruit fly", human = FALSE)
  str_312 <- orthologs(genes = "312", species = "fruit fly", human = FALSE)
  id_312 <- orthologs(genes = 38161, species = "fruit fly", human = FALSE)
  expect_equal(nrow(int_312), 1)
  expect_identical(int_312, str_312)
  expect_identical(int_312, id_312)
  # gene symbol that can also be Ensembl ID (not present in HCOP yet)
  ens_mix <- orthologs(genes = c("ENSG00000083622", "ENSG00000145063", "ENSA"), species = "mouse")
  expect_equal(nrow(ens_mix), 1)
  # multiple Ensembl IDs corresponding to the same symbol
  ndst_sym_top <- orthologs(genes = "NDST2", species = "mouse")
  ndst_ens_all <- orthologs(genes = "ENSG00000166507", species = "mouse", top = FALSE)
  ndst_sym_all <- orthologs(genes = "NDST2", species = "mouse", top = FALSE)
  expect_identical(ndst_sym_top, ndst_ens_all)
  expect_lt(nrow(ndst_sym_top), nrow(ndst_sym_all))
  expect_lt(nrow(ndst_ens_all), nrow(ndst_sym_all))
})

test_that("orthologs species", {
  o_m <- orthologs(genes = "PTPRC", species = "mouse")
  o_r <- orthologs(genes = "PTPRC", species = "rat")
  o_z <- orthologs(genes = "PTPRC", species = "zebrafish")
  expect_equal(nrow(o_m), 1)
  expect_equal(nrow(o_r), 1)
  expect_equal(nrow(o_z), 1)
  expect_identical(o_m[, 1:3], o_r[, 1:3])
  expect_identical(o_m[, 1:3], o_z[, 1:3])
})

test_that("orthologs non-human input", {
  ptprc_m <- orthologs(genes = "Ptprc", species = "mouse", human = FALSE)
  ptprc_r <- orthologs(genes = "Ptprc", species = "rat", human = FALSE)
  ptprc_z <- orthologs(genes = "ptprc", species = "zebrafish", human = FALSE)
  expect_equal(nrow(ptprc_m), 1)
  expect_equal(nrow(ptprc_r), 1)
  expect_equal(nrow(ptprc_z), 1)
  expect_identical(ptprc_m[, 1:3], ptprc_r[, 1:3])
  expect_identical(ptprc_m[, 1:3], ptprc_z[, 1:3])
})

test_that("orthologs top", {
  gapdh_top <- orthologs(genes = "GAPDH", species = "mouse")
  gapdh_top_t <- orthologs(genes = "GAPDH", species = "mouse", top = TRUE)
  gapdh_top_f <- orthologs(genes = "GAPDH", species = "mouse", top = FALSE)
  expect_equal(nrow(gapdh_top), nrow(gapdh_top_t))
  expect_equal(nrow(gapdh_top_t), 1)
  expect_lt(nrow(gapdh_top_t), nrow(gapdh_top_f))
})

test_that("orthologs min_support", {
  lilrb_default <- orthologs(genes = "LILRB3", species = "mouse", top = FALSE)
  lilrb_min1 <- orthologs(genes = "LILRB3", species = "mouse", min_support = 1, top = FALSE)
  lilrb_min3 <- orthologs(genes = "LILRB3", species = "mouse", min_support = 3, top = FALSE)
  lilrb_min5 <- orthologs(genes = "LILRB3", species = "mouse", min_support = 5, top = FALSE)
  lilrb_min9 <- orthologs(genes = "LILRB3", species = "mouse", min_support = 9, top = FALSE)
  expect_equal(nrow(lilrb_min5), 1)
  expect_equal(nrow(lilrb_min9), 0)
  expect_equal(nrow(lilrb_min3), nrow(lilrb_default))
  expect_lt(nrow(lilrb_default), nrow(lilrb_min1))
})

test_that("orthologs wrong input", {
  expect_error(orthologs(genes = ".", species = "mouse"))
  expect_error(orthologs(genes = "?", species = "mouse"))
  expect_error(orthologs(genes = "xxxxx", species = "mouse"))
  expect_error(orthologs(genes = "xxxxx", species = "mouse", human = FALSE))
  expect_error(orthologs(genes = "PTPRC", species = "xxxxx"))
  expect_error(orthologs(genes = "PTPRC", species = 10090))
  expect_error(orthologs(genes = "PTPRC", species = "mouse", human = "xxxxx"))
  expect_error(orthologs(genes = "PTPRC", species = "mouse", top = "xxxxx"))
  expect_error(orthologs(genes = c("PTPRC", "ENSG00000081237"), species = "mouse"))
  expect_error(orthologs(genes = c("Ptprc", "ENSMUSG00000026395"), species = "mouse", human = FALSE))
  expect_error(orthologs(genes = data.frame(gene = "PTPRC"), species = "mouse", top = "xxxxx"))
})

test_that("species basics", {
  sp <- species()
  expect_s3_class(sp, "data.frame")
  expect_equal(nrow(sp), 19)
  spm <- species("mouse")
  expect_s3_class(spm, "data.frame")
  expect_equal(nrow(spm), 1)
  spr <- species("rat")
  expect_s3_class(spr, "data.frame")
  expect_equal(nrow(spr), 1)
  expect_error(species(10090))
})
