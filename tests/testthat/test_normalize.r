library(cellrangerRkit)
context('Normalization')

gene_names = c('a', 'b')
sample_names = c('sample1', 'sample2', 'sample3')

# Generate a test dataset
test_matrix = Matrix(matrix(c(1,2,10,20,100,200), ncol=3))
row.names(test_matrix) = gene_names
colnames(test_matrix) = sample_names

test_genes = data.frame(id=gene_names, symbol=gene_names)
row.names(test_genes) = gene_names

test_sample_data = data.frame(barcode=sample_names)
row.names(test_sample_data) = sample_names

# Create GBM from it
test_gbm <- newGeneBCMatrix(test_matrix, test_genes, test_sample_data)

# Test median normalization
test_gbm_bcnorm <- normalize_barcode_sums_to_median(test_gbm)

test_that("bc normalization gives correct result", {
  expect_true(all(colSums(test_gbm_bcnorm) == 30))
})

# Test log normalization
test_gbm <- newGeneBCMatrix(test_matrix, test_genes, test_sample_data)
test_gbm_log2 <- log_gene_bc_matrix(test_gbm, base=2)
test_that("bc normalization gives correct result", {
  expect_true(all(test_gbm_log2 == log2(1+test_gbm)))
})
