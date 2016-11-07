library(cellrangerRkit)
context('Concatenation')

pipestance_path <- system.file("extdata", "test_pipestance", package="cellrangerRkit")

filt_gbm <- load_cellranger_matrix(pipestance_path, barcode_filtered=T)

print("loaded")
gbm1 <- filt_gbm[,1:30]
gbm2 <- filt_gbm[,31:60]
gbm3 <- filt_gbm[,61:114]

big_gbm <- concatenate_gene_bc_matrices(list(gbm1, gbm2, gbm3))

test_that("concatenated gbm has correct content", {
  expect_true(nrow(big_gbm) == nrow(filt_gbm))
  expect_true(all(fData(big_gbm) == fData(filt_gbm)))
  expect_true(ncol(big_gbm) == ncol(gbm1) + ncol(gbm2) + ncol(gbm3))
  expect_true(sum(big_gbm) == sum(gbm1) + sum(gbm2) + sum(gbm3))
})

test_that("can't concat different gene sets", {
  expect_error(concatenate_gene_bc_matrices(list(gbm1, gbm1[1:100,])))
})

gbm4 <- filt_gbm[,1:50]
test_that("can concatenate intersecting barcode sets", {
  expect_error(big_gbm <- concatenate_gene_bc_matrices(list(gbm1, gbm4), NA))

  # TODO For whatever reason, these subsetted matrices have diff sample names now
  # Not sure why, but the sample names are same length, but differ between exprs
  # and phenoData during this call between gbm1 and gbm4
  # I think the reason may have to do with fact that barcodes should actually be
  # the columns in an expression dataset typically, but are rep as rows here
  expect_true(ncol(concatenate_gene_bc_matrices(list(gbm1, gbm4))) == ncol(gbm1) + ncol(gbm4))
})
