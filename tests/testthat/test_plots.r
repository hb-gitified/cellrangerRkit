library(cellrangerRkit)
context('Plotting')

pipestance_path <- system.file("extdata", "test_pipestance", package="cellrangerRkit")

gbm <- load_cellranger_matrix(pipestance_path)

test_that("can plot bc counts", {
  plot_barcode_counts(gbm)
})
