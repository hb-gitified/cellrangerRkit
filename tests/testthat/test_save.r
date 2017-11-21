library(cellrangerRkit)
context('Save matrix H5')

pipestance_path <- system.file("extdata", "test_pipestance_cr130", package="cellrangerRkit")
filt_gbm <- load_cellranger_matrix(pipestance_path, barcode_filtered=T)

test_that("filtered matrix has correct dimensions", {
  expect_equal(nrow(filt_gbm), 33694)
  expect_equal(ncol(filt_gbm), 4122)
})

filt_gbm <- filt_gbm[,1:(ncol(filt_gbm)-100)]
if (file.exists('test_save.h5')) {
  unlink('test_save.h5')
}
save_cellranger_matrix_h5(filt_gbm, 'test_save.h5', 'GRCh38')
reloaded_gbm <- get_matrix_from_h5('test_save.h5')
test_that("saved matrix has correct dimensions", {
  expect_equal(nrow(reloaded_gbm), 33694)
  expect_equal(ncol(reloaded_gbm), 4122-100)
})

test_that("sums are equal", {
  expect_equal(sum(filt_gbm), sum(reloaded_gbm))
})


