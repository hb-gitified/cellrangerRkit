library(cellrangerRkit)
context('differential expression')

set.seed(0)

# Simulate a gene expression matrix
n_cells <- 1e3
n_transcripts <- 10000
n_genes <- 10000
plaw_alpha <- 0.7
n_de_genes <- 50
l2fc <- 2

txome0 <- (1:n_genes)^-plaw_alpha
txome0 <- txome0/sum(txome0)
de_idx <- sample(length(txome0), n_de_genes)
txome1 <- txome0
txome1[de_idx] <- l2fc * txome1[de_idx]
txome1 <- txome1/sum(txome1)

size_factors <- exp(rnorm(n_cells, 0, 0.1))

test_matrix <- do.call(rbind, lapply(1:(n_cells/2), function(i) rpois(n_genes, size_factors[i]*n_transcripts*txome0)))
test_matrix <- rbind(test_matrix, do.call(rbind, lapply((n_cells/2+1):n_cells, function(i) rpois(n_genes, size_factors[i]*n_transcripts*txome1))))
test_matrix <- t(test_matrix)

row.names(test_matrix) = sprintf("gene%d", 1:n_genes)
colnames(test_matrix) = sprintf("cell%d", 1:n_cells)

genes = data.frame(id=sprintf("gene%d", 1:n_genes), symbol=sprintf("symbol%d", 1:n_genes))
row.names(genes) = genes$id

cells = data.frame(barcode=sprintf("cell%d", 1:n_cells))
row.names(cells) = cells$barcode

gbm <- newGeneBCMatrix(test_matrix, fd=genes, pd=cells)

de_result <- run_differential_expression(gbm, 1:(n_cells/2), (1+n_cells/2):n_cells)

act_pos <- 1:n_genes %in% de_idx
called_pos <- de_result$p_adj < 0.05
tp <- act_pos & called_pos

sens <- sum(tp) / sum(act_pos)
ppv <- sum(tp) / sum(called_pos)

test_that("DE results make sense", {
  expect_equal(sens, 1)
  expect_equal(ppv, 1)
} )


