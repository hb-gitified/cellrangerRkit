library(cellrangerRkit)
context('GeneBCMatrix')

pipestance_path <- system.file("extdata", "test_pipestance", package="cellrangerRkit")

filt_gbm <- load_cellranger_matrix(pipestance_path, barcode_filtered=T)
test_that("filtered matrix has correct dimensions", {
  expect_equal(nrow(filt_gbm), 35635)
  expect_equal(ncol(filt_gbm), 114)
})

test_that("bc filter flag is set correctly", {
  expect_equal(filt_gbm@barcode_filtered, TRUE)
})

# Try subselecting barcodes
filt_gbm <- filt_gbm[,1:50]
test_that("bc-subselected matrix has correct dimensions", {
  expect_equal(nrow(filt_gbm), 35635)
  expect_equal(ncol(filt_gbm), 50)
})

# Try subselecting genes
filt_gbm <- filt_gbm[1:500,]
test_that("gene-subselected matrix has correct dimensions", {
  expect_equal(nrow(filt_gbm), 500)
  expect_equal(ncol(filt_gbm), 50)
  expect_equal(length(row.names(fData(filt_gbm))), 500)
})

test_that("can load secondary analysis results", {
  load_cellranger_analysis_results(pipestance_path)
})

