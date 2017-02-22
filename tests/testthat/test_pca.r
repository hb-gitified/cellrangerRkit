library(cellrangerRkit)
context('PCA')

test_matrix = Matrix(matrix(rbinom(200*100, 10, 0.1), 200, 100))

grp1 <- 1:40
grp2 <- 41:100
test_matrix[1:5,grp1] <- 5
test_matrix[6:10,grp2] <- 6

row.names(test_matrix) = sprintf("gene%d", 1:200)
colnames(test_matrix) = sprintf("sample%d", 1:100)

test_genes = data.frame(id=sprintf("gene%d", 1:200))
row.names(test_genes) = test_genes$id

test_samples = data.frame(barcode=sprintf("sample%d", 1:100))
row.names(test_samples) = test_samples$barcode

gbm <- newGeneBCMatrix(test_matrix, test_genes, test_samples)

pca <- run_pca(gbm, 2)

test_that("pca results make sense", {
  expect_true(nrow(pca$x) == ncol(gbm))
  expect_equal(ncol(pca$x), 2)
  expect_true(nrow(pca$rotation) == length(pca$use_genes))
  expect_equal(ncol(pca$rotation), 2)
  expect_true(all(pca$x[grp1,1] < 0) && all(pca$x[grp2,1] > 0) ||
              all(pca$x[grp1,1] > 0) && all(pca$x[grp2,1] < 0))
} )
