library(cellrangerRkit)
context('Util')


gene_names = c('gene_a', 'gene_b', 'gene_c')
sample_names = c('sample_a', 'sample_b', 'sample_c')

test_matrix = Matrix(matrix(c(1,0,3,10,0,30,100,0,300), ncol=3))
row.names(test_matrix) = gene_names
colnames(test_matrix) = sample_names

test_gene_matrix = data.frame(id=gene_names, symbol=gene_names)
row.names(test_gene_matrix) = test_gene_matrix$id

test_samples = data.frame(barcode=sample_names)
row.names(test_samples) = test_samples$barcode

test_gbm <- newGeneBCMatrix(test_matrix, test_gene_matrix, test_samples)
test_genes <- get_nonzero_genes(test_gbm)


test_that("get_nonzero_genes gives correct result", {
  expect_equal(length(test_genes), 2)
  expect_equal(test_genes[1], 1)
  expect_equal(test_genes[2], 3)
})

pipestance_path <- system.file("extdata", "test_pipestance", package="cellrangerRkit")
filt_gbm <- load_cellranger_matrix(pipestance_path, barcode_filtered=T)

test_that("can get mean reads per cell", {
  get_mean_mapped_reads_per_cell(filt_gbm)
})
