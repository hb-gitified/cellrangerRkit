#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

#' Normalization functions
#'

#' Normalize barcodes by sum to the median barcode by sum
#'
#' @param gbm a GeneBCMatrix to normalize
#' @return A GeneBCMatrix with barcode-normalized values
#' @export
#' @examples
#' gbm_norm_by_bc <- normalize_barcode_sums_to_median(gbm1)
#'
normalize_barcode_sums_to_median <- function(gbm) {
  bc_sums <- colSums(gbm)
  median_sum <- median(bc_sums)
  new_matrix <- sweep(exprs(gbm), 2, median_sum/bc_sums, '*')

  return(newGeneBCMatrix(new_matrix, fData(gbm), pData(gbm), template=gbm))
}

#' Take the element-wise logarithm of a gene-barcode matrix while maintaining sparsity
#'
#' @param gbm a GeneBCMatrix to log
#' @param base base of logarithm to use
#' @return A GeneBCMatrix with logged counts
#' @export
#' @examples
#' gbm_log <- log_gene_bc_matrix(gbm1)
#'
log_gene_bc_matrix <- function(gbm, base=2) {
  x <- uniqTsparse(as(exprs(gbm), 'dgTMatrix'))
  slot(x, 'x') <- log(1 + slot(x, 'x'), base=base)
  return(newGeneBCMatrix(x, fData(gbm), pData(gbm), template=gbm))
}
