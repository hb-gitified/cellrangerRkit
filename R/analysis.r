#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

#' Analysis functions
#'

#' Run principal components analysis to reduce the dimensionality of a gene-barcode matrix
#'
#' @param gbm a GeneBCMatrix to run PCA on
#' @param n_pcs The number of principal components to compute
#' @return A list containing \itemize{
#' \item "x" - The rotated data matrix where rows are barcodes and columns are PCs
#' \item "rotation" - The loadings (eigenvectors) where each column is a PC
#' \item "use_genes" - The list of gene indices that were used
#' \item "normalized_mat" - The matrix that was used as input to the PCA
#' }
#' @export
#' @examples
#' pca <- run_pca(gbm1, n_pcs=2)
#' plot(PC2 ~ PC1, pca$x)
#'
run_pca <- function(gbm, n_pcs=10) {
  use_genes <- get_nonzero_genes(gbm)
  gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
  gbm_log <- log_gene_bc_matrix(gbm_bcnorm)

  pca <- sparse_pca(t(exprs(gbm_log)), n_pcs)

  pc_names <- sprintf("PC%d", 1:ncol(pca$x))
  rownames(pca$x) <- colnames(gbm)
  colnames(pca$x) <- pc_names

  rownames(pca$rotation) <- rownames(gbm)[use_genes]
  colnames(pca$rotation) <- pc_names
  cat('Variance explained by PCs',pca$var_explained,'\n')
  return(list(x=pca$x,
              rotation=pca$rotation,
              use_genes=use_genes,
              normalized_mat=exprs(gbm_log),
              var_explained=pca$var_explained))
}

#' Project dimensionally reduced gene-barcode matrix into a visualizable space via t-SNE
#'
#' @param pca A PCA result from run_pca()
#' @param ... parameters to pass to Rtsne
#' @return See ?Rtsne in the Rtsne package. The item of interest is the 'Y' element which contains the coordinates.
#' @importFrom Rtsne Rtsne
#' @export
#' @examples
#' \dontrun{
#' pca <- run_pca(gbm1, n_pcs=10)
#' tsne_res <- run_tsne(pca)
#' plot(tsne_res$Y[,1], tsne_res$Y[,2])
#'}
run_tsne <- function(pca, ...) {
  return(Rtsne(pca$x, ...))
}

#' Run k-means clustering
#'
#' @param pca A PCA result from run_pca()
#' @param k number of clusters
#' @param ... parameters to pass to kmeans
#' @return See ?kmeans. The item of interest is the 'cluster' element
#' @export
#' @seealso kmeans
#' @examples
#' pca <- run_pca(gbm1, n_pcs=10)
#' km_res <- run_kmeans_clustering(pca, k=5)
#'
run_kmeans_clustering <- function(pca, k, ...) {
  return(kmeans(pca$x, k, ...))
}

#' Run a fast PCA algorithm while maintaining sparsity
#'
#' @param x a sparse matrix to run PCA on
#' @param n_pcs number of principal components to compute
#' @param mu column means
#' @param s column standard deviations
#' @param center_scale perform centering and scaling
#' @return a list containing the rotated matrix ('x') and the principal component vectors ('rotation')
#' @importFrom irlba irlba
#' @examples
#' \dontrun{
#' pca <- sparse_pca(x, n_pcs=10)
#' plot(pca$x[,1], pca$x[,2], xlab='PC1', ylab='PC2')
#' }
sparse_pca <- function(x, n_pcs, mu=NULL, s=NULL, center_scale=TRUE) {
  if (is.null(mu) && center_scale) mu <- colMeans(x)
  if (is.null(s) && center_scale) s <- apply(x, 2, sd, na.rm=TRUE)

  if (center_scale) {
    s[s == 0] <- min(s[s > 0])
    svd_res <- irlba::irlba(x, n_pcs, center=mu, scale=s)
  } else {
    svd_res <- irlba::irlba(x, n_pcs)
  }
  # compute explained variance
  variance_sum <- sum(apply(x,2,var,na.rm=TRUE)) # sample variance sum
  pca_sum <- sum(svd_res$d^2)/(dim(x)[1]-1)
  var_explained <- pca_sum / variance_sum

  return(list(x=svd_res$u %*% diag(svd_res$d), rotation=svd_res$v, var_explained=var_explained))
}

#' Get ordered gene list by mean difference between clusters
#'
#' @param gbm A GeneBCMatrix
#' @param cluster_labels An array of integer labels from 1 to k (length = number of cells)
#' @param test_ID A label between 1 to k to be tested
#' @return an array of genes sorted by the mean differential expression between testID and the rest
#'
get_ordered_gene_list_mean_diff <- function(gbm,cluster_labels,test_ID) {
  # return gene list sorted by mean difference between a cluster verse the rest
  cat('Comparing Cluster',test_ID, 'against the other clusters...\n')
  mean_0 <- rowMeans(gbm[,cluster_labels!=test_ID])
  mean_1 <- rowMeans(gbm[,cluster_labels==test_ID])
  mean_diff <- mean_1 - mean_0
  sort_index <- sort(mean_diff,index.return = TRUE, decreasing=TRUE)$ix
  return(sort_index)
}

#' Get ordered gene list by gene uniqueness
#'
#' @param gbm A GeneBCMatrix
#' @param cluster_labels An array of integer labels from 1 to k (length = number of cells)
#' @param test_ID A label between 1 to k to be tested
#' @return an array of genes sorted by the mean differential expression between testID and the rest
#'
get_ordered_gene_list_unique <- function(gbm,cluster_labels,test_ID) {
  # return gene list sorted by mean difference between a cluster verse the rest
  cat('Comparing Cluster',test_ID, 'against the other clusters...\n')
  mean_0 <- rowMeans(exprs(gbm[,cluster_labels!=test_ID]) > 0 )
  mean_1 <- rowMeans(exprs(gbm[,cluster_labels==test_ID]) > 0 )
  mean_diff <- mean_1 - mean_0
  sort_index <- sort(mean_diff,index.return = TRUE, decreasing=TRUE)$ix
  return(sort_index)
}

#' Order cells based on clusters
#'
#' @param gbm A GeneBCMatrix
#' @param clu An array of integer labels from 1 to k (length = number of cells)
#' @return A list of (n_top_genes * number of clusters) genes (may include duplicates)
#' @import Biobase
#' @export
#' @examples
#' \dontrun{
#'  order_cell_by_clusters(gbm_log,clu)
#' }
order_cell_by_clusters <- function(gbm,clu) {
  uni_labels <- sort(unique(clu))
  sorted_clu_id <- sort(clu,index.return=TRUE)
  index_list <- lapply(uni_labels, function(x) sorted_clu_id$ix[sorted_clu_id$x==x])
  cell_info <- as.character(pData(gbm)$barcode)
  ordered_cells <- lapply(index_list, function(x) {
                          list(barcode=cell_info[x],ix=x)} )
  names(ordered_cells) <- uni_labels
  return(ordered_cells)
}

#' Prioritize genes that are most significant markers for each cluster
#'
#' @param gbm A GeneBCMatrix (recommended to be normalized by UMI counts and log-transformed)
#' @param clu An array of integer labels from 1 to k (length = number of cells)
#' @param n_top_genes Number of top genes for each cluster group to consider
#' @param method Method to priortize genes: currently only supports mean difference
#' @import Biobase
#' @export
#' @examples
#' \dontrun{
#'  prioritize_top_genes(gbm_log,clu,3,'mean_diff')
#' }
prioritize_genes<-function(gbm,clu,method) {
  uni_labels <- sort(unique(clu))
  if (method == 'mean_diff') {
    index_list <- lapply(uni_labels, function(x) get_ordered_gene_list_mean_diff(gbm,clu,x))
  } else if (method == 'unique') {
    index_list <- lapply(uni_labels, function(x) get_ordered_gene_list_unique(gbm,clu,x))
  } else {
    stop('Method: ',method,' not recognized\n')
  }
  genes <- fData(gbm)
  ordered_genes <- lapply(index_list, function(x) {
                          list(id=genes$id[x],symbol=genes$symbol[x],ix=x)} )
  names(ordered_genes) <- uni_labels
  return(ordered_genes)
}

#' Tabulate the proportion of each cell type in respective classes
#'
#' @param clu An array of integer labels from 1 to k (length = number of cells)
#' @param anno Annotation for the k classes (can be an array of strings)
#' @param n_digits Number of digits for to display
#' @return A table of the break down of the numbers of each label (with optional annotation)
#' @export
#' @examples
#' \dontrun{
#'  cell_composition(clu)
#' }
cell_composition <- function(clu,anno=NULL,n_digits=5) {
  cat("Cell composition: \n")
  num_cells <- table(factor(clu,levels=min(clu):max(clu)))
  num_cells <- round(num_cells, 0)
  proportion <- num_cells / sum(num_cells)
  proportion <- round(proportion, n_digits)
  out <- rbind(num_cells,proportion)
  if (!is.null(labels)) {
    # colnames(out) <- anno
    annotation <- anno
    out <- rbind(annotation,out)
  }
  return(out)
}
