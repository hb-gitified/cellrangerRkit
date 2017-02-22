library(cellrangerRkit)
context('kmeans')

set.seed(0)
test_matrix = Matrix(matrix(rbinom(200*100, 10, 0.1), 200, 100))

grp1 <- 1:40
grp2 <- 41:100
test_matrix[1:5,grp1] <- 5
test_matrix[6:10,grp2] <- 6

row.names(test_matrix) = sprintf("gene%d", 1:200)
colnames(test_matrix) = sprintf("cell%d", 1:100)

genes = data.frame(id=sprintf("gene%d", 1:200), symbol=sprintf("symbol%d", 1:200))
row.names(genes) = genes$id

cells = data.frame(barcode=sprintf("cell%d", 1:100))
row.names(cells) = cells$barcode

gbm <- newGeneBCMatrix(test_matrix, fd=genes, pd=cells)

pca <- run_pca(gbm, 2)
set.seed(0)
km <- run_kmeans_clustering(pca, 2, nstart=10)

test_that("kmeans works", {
  expect_true(all(km$cluster[grp1] == 1) && all(km$cluster[grp2] == 2) ||
                all(km$cluster[grp1] == 2) && all(km$cluster[grp2] == 1))
} )
