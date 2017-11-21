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

# Try loading H5
pipestance_path_h5 <- system.file("extdata", "test_pipestance_h5", package="cellrangerRkit")

filt_gbm_h5 <- load_cellranger_matrix_h5(pipestance_path_h5, barcode_filtered=T)
test_that("filtered matrix has correct dimensions", {
  expect_equal(nrow(filt_gbm_h5), 343)
  expect_equal(ncol(filt_gbm_h5), 34)
})

# Test on a Cell Ranger 1.3.0 pipestance
pipestance_path_cr130 <- system.file("extdata", "test_pipestance_cr130", package="cellrangerRkit")
filt_gbm_cr130 <- load_cellranger_matrix(pipestance_path_cr130, barcode_filtered=T)
test_that("filtered matrix has correct dimensions", {
  expect_equal(nrow(filt_gbm_cr130), 33694)
  expect_equal(ncol(filt_gbm_cr130), 4122)
})

an <- load_cellranger_analysis_results(pipestance_path_cr130)
test_that("can load secondary analysis results (cr 1.3.0)", {
  expect_true(!is.null(an$clustering$graphclust))
  expect_true(!is.null(an$clustering$kmeans_2_clusters))
})

# Test matrix H5 loading
filt_gbm_h5_cr130 <- load_cellranger_matrix_h5(pipestance_path_cr130, barcode_filtered=T)
test_that("filtered matrix H5 has correct dimensions", {
  expect_equal(nrow(filt_gbm_h5_cr130), 33694)
  expect_equal(ncol(filt_gbm_h5_cr130), 4122)
})
