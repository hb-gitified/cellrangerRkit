library(cellrangerRkit)
context('Subsampling')

pipestance_path <- system.file("extdata", "test_subsample",
                               package="cellrangerRkit")

set.seed(0)

gbm <- load_cellranger_matrix(pipestance_path)
gbm <- load_molecule_info(gbm)

gbm_ss <- subsample_gene_bc_matrix(gbm, 1000000)
test_that("subsampling returns orig matrix if target depth higher", {
  expect_equal(sum(gbm_ss), sum(gbm))
})

rpc <- get_mean_raw_reads_per_cell(gbm)

gbm_ss <- subsample_gene_bc_matrix(gbm, floor(rpc))
test_that("subsampling returns orig matrix at same depth", {
  expect_equal(sum(gbm_ss), sum(gbm))
})

gbm_ss <- subsample_gene_bc_matrix(gbm, floor(rpc/2))
test_that("subsampling to 1/2", {
  expect_true(abs(round(sum(gbm_ss)/sum(gbm),2) - 0.5) < 0.01)
})

gbm_ss <- subsample_gene_bc_matrix(gbm, floor(rpc/10))
test_that("subsampling to 1/10", {
  expect_true(abs(round(sum(gbm_ss)/sum(gbm),2) - 0.1) < 0.01)
})


