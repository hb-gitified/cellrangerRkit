#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

#' Check for a broken version of irlba
check_irlba <- function() {
  irlba_ver <- try(packageVersion('irlba'), silent=T)
  if (inherits(irlba_ver, 'try-error')) {
    warning("Could not find 'irlba' package.")
  } else if (irlba_ver == '2.1.1') {
    warning("Your installed version of the 'irlba' package (2.1.1) has a critical bug in it that causes an error when you try to run PCA. Please see https://github.com/bwlewis/irlba/issues/6 for the issue thread.")
  }
}

check_irlba()

#' Analysis functions
#'

#' Run principal components analysis to reduce the dimensionality of a gene-barcode matrix
#'
#' @param gbm a GeneBCMatrix to run PCA on (must set logscale=TRUE if it is already log-transformed)
#' @param n_pcs The number of principal components to compute (default: 10)
#' @param logscale Logical if the input is log-scale or not (default: FALSE)
#' @return A list containing \itemize{
#' \item "x" - The rotated data matrix where rows are barcodes and columns are PCs
#' \item "sdev" - the standard deviations of the principal components (i.e., sqrt of eigvals of the covariance matrix)
#' \item "rotation" - The loadings (eigenvectors) where each column is a PC
#' \item "tot_var" - The total variation in the scaled and centered matrix (this is also the effective rank of the matrix)
#' \item "var_pcs" - The proportion of variance explained by each principle comoponent
#' \item "use_genes" - The list of gene indices that were used
#' \item "normalized_mat" - The matrix that was used as input to the PCA
#' }
#' @export
#' @examples
#' pca <- run_pca(gbm1, n_pcs=2)
#' plot(PC2 ~ PC1, pca$x)
#'
run_pca <- function(gbm, n_pcs=10, logscale=FALSE, ...) {
  if (logscale) {
    gbm_log <- gbm
  } else {
    use_genes <- get_nonzero_genes(gbm)
    gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
    gbm_log <- log_gene_bc_matrix(gbm_bcnorm, ...)
  }

  pca <- sparse_pca(t(exprs(gbm_log)), n_pcs)

  pc_names <- sprintf("PC%d", 1:ncol(pca$x))
  rownames(pca$x) <- colnames(gbm)
  colnames(pca$x) <- pc_names

  rownames(pca$rotation) <- rownames(gbm)[use_genes]
  colnames(pca$rotation) <- pc_names
  var_explained <- sum(pca$var_pcs)
  cat('Variance explained by PCs:',var_explained,'\n')
  return(list(x=pca$x,
              rotation=pca$rotation,
              sdev=pca$sdev,
              tot_var=pca$tot_var,
              var_pcs=pca$var_pcs,
              use_genes=use_genes,
              normalized_mat=exprs(gbm_log),
              var_explained=var_explained))
}


#' Project dimensionally reduced gene-barcode matrix into a visualizable space via t-SNE
#'
#' @param pca_res A PCA result from run_pca()
#' @param n_pcs The number of top pca components to use (default is to use all pcs)
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
run_tsne <- function(pca_res, n_pcs=NULL, ...) {
  if (!is.null(n_pcs)) {
    n_p <- min(dim(pca_res$x)[2], n_pcs)
  } else {
    n_p <- dim(pca_res$x)[2]
  }
  pca_loadings <- pca_res$x[,1:n_p]
  cat('Using number of of pcs:',n_p,'\n')
  return(Rtsne(pca_loadings, pca=FALSE, ...))
}

#' Run k-means clustering
#'
#' @param pca A PCA result from run_pca()
#' @param n_pcs The number of top pca components to use (default is to use all pcs)
#' @param k number of clusters
#' @param ... parameters to pass to kmeans
#' @return See ?kmeans. The item of interest is the 'cluster' element
#' @export
#' @seealso kmeans
#' @examples
#' pca <- run_pca(gbm1, n_pcs=10)
#' km_res <- run_kmeans_clustering(pca, k=5)
#'
run_kmeans_clustering <- function(pca, k, n_pcs=NULL, ...) {
  if (!is.null(n_pcs)) {
    n_p <- min(dim(pca$x)[2], n_pcs)
  } else {
    n_p <- dim(pca$x)[2]
  }
  pca_loadings <- pca$x[,1:n_p]
  cat('Using number of of pcs:', n_p, '\n')
  return(kmeans(pca_loadings, k, ...))
}


#' Run a fast PCA algorithm while maintaining sparsity
#'
#' @param x a sparse matrix to run PCA on
#' @param n_pcs number of principal components to compute
#' @param mu column means
#' @param s column standard deviations
#' @param center_scale perform centering and scaling
#' @return A list containing \itemize{
#' \item "x" - The rotated data matrix where rows are barcodes and columns are PCs
#' \item "sdev" - the standard deviations of the principal components (i.e., sqrt of eigvals of the covariance matrix)
#' \item "rotation" - The loadings (eigenvectors) where each column is a PC
#' \item "tot_var" - The total variation in the scaled and centered matrix (this is also the effective rank of the matrix)
#' \item "var_pcs" - The proportion of variance explained by each principle comoponent
#' }
#' @importFrom irlba irlba
#' @examples
#' \dontrun{
#' pca <- sparse_pca(x, n_pcs=10)
#' plot(pca$x[,1], pca$x[,2], xlab='PC1', ylab='PC2')
#' }
sparse_pca <- function(x, n_pcs, mu=NULL, s=NULL, center_scale=TRUE) {
  if (is.null(mu) && center_scale) mu <- colMeans(x)
  if (is.null(s) && center_scale) s <- apply(x, 2, sd, na.rm=TRUE)

  # irlba version 2.1.+ is borked, but 'fastpath=FALSE' is a workaround.
  # Older versions of irlba don't need the workaround, but don't have the arg.
  irlba_wrapper <- function(...) {
    # Fails if irlba version is pre-fastpath
    res <- try(irlba::irlba(fastpath=FALSE, ...), silent=T)
    if (inherits(res, 'try-error')) {
      # Pre-fastpath, so don't bother with option
      res <- irlba::irlba(...)
    }
    res
  }

  if (center_scale) {
    s[s == 0] <- min(s[s > 0])
    svd_res <- irlba_wrapper(x, n_pcs, center=mu, scale=s)
  } else {
    svd_res <- irlba_wrapper(x, n_pcs)
  }

  # compute explained variance
  n <- dim(x)[1]
  variance_sum <- sum(apply(x,2,var,na.rm=TRUE)/(s^2)) # sample variance sum
  var_pcs <- svd_res$d^2/(n-1)/variance_sum

  return(list(x=svd_res$u %*% diag(svd_res$d), rotation=svd_res$v, sdev=svd_res$d/sqrt(n-1),
              tot_var=variance_sum, var_pcs=var_pcs))
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

#' Get ordered gene list by the sSeq differential expression method
#'
#' @param gbm A GeneBCMatrix
#' @param cluster_labels An array of integer labels from 1 to k (length = number of cells)
#' @param test_ID A label between 1 to k to be tested
#' @param sseq_params List returned by compute_sseq_params
#' @param t_mat BC x gene matrix
#' @return an array of genes sorted by the mean differential expression between testID and the rest
#'
get_ordered_gene_list_sseq <- function(gbm, cluster_labels, test_ID, sseq_params,
                                       t_mat) {
  # return gene list sorted by L2FC in a cluster vs the rest
  cat('Comparing Cluster', test_ID, 'against the other clusters...\n')
  group0 <- cluster_labels == test_ID
  group1 <- cluster_labels != test_ID

  de_result <- sseq_differential_expression(t_mat, group0, group1,
                                 sseq_params,
                                 gene_ids=fData(gbm)$id,
                                 gene_symbols=fData(gbm)$symbol)
  return(de_result)
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
order_cell_by_clusters <- function(gbm, clu) {
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
#' @param gbm A GeneBCMatrix (recommended to be normalized by UMI counts and log-transformed for the  "mean_diff" method). For the "sseq" method, this should be a raw (unnormalized) GeneBCMatrix.
#' @param clu An array of integer labels from 1 to k (length = number of cells)
#' @param logscale Logical if the input is log-scale or not (default: TRUE)
#' @param n_top_genes Number of top genes for each cluster group to consider
#' @param method Method to prioritize genes: currently only supports mean difference ('mean_diff', 'unique', or 'sseq')
#' @param p_cutoff If method is "sseq," only consider genes w/ adjusted p-value below this value
#' @param min_mean If method is "sseq," only consider genes w/ mean normalized UMI counts/cell exceeding this value
#' @param min_log2fc If method is "sseq," only consider genes with at least this log2 fold-change
#' @param order_by If method is "sseq," sort genes by 'pvalue' or by 'log2fc'
#' @import Biobase
#' @export
#' @examples
#' \dontrun{
#'  prioritize_top_genes(gbm, cluster_labels, 'sseq')
#' }
prioritize_top_genes <- function(gbm, clu, method, logscale=TRUE,
                                 p_cutoff=0.05, min_mean=1, min_log2fc=0,
                                 order_by='pvalue') {
  uni_labels <- sort(unique(clu))

  if (!(order_by %in% c('pvalue', 'log2fc'))) {
    stop('Order-by: ', order_by, ' not recognized\n')
  }

  if (method %in% c('mean_diff, unique')) {
    if (logscale) {
      gbm_log <- gbm
    } else {
      use_genes <- get_nonzero_genes(gbm)
      gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
      gbm_log <- log_gene_bc_matrix(gbm_bcnorm)
    }
  }

  if (method == 'mean_diff') {
    index_list <- lapply(uni_labels, function(x) get_ordered_gene_list_mean_diff(gbm_log, clu, x))
    genes <- fData(gbm_log)
    ordered_genes <- lapply(index_list, function(x) {
      list(id=genes$id[x], symbol=genes$symbol[x], ix=x)} )

  } else if (method == 'unique') {
    index_list <- lapply(uni_labels, function(x) get_ordered_gene_list_unique(gbm_log, clu, x))
    genes <- fData(gbm_log)
    ordered_genes <- lapply(index_list, function(x) {
      list(id=genes$id[x], symbol=genes$symbol[x], ix=x)} )

  } else if (method == 'sseq') {
    t_mat <- t(exprs(gbm))
    cat('Computing differential expression parameters...\n')
    sseq_params <- compute_sseq_params(t_mat)
    ordered_genes <- lapply(uni_labels, function(cluster_label) {
      de_result <- get_ordered_gene_list_sseq(gbm, clu, cluster_label, sseq_params, t_mat)
      de_result$ix <- 1:nrow(de_result)
      de_result$id <- de_result$gene_id
      de_result$symbol <- de_result$gene_name
      de_result$significant <- with(de_result, common_mean >= min_mean &
                                      p_adj < p_cutoff &
                                      log2fc >= min_log2fc)
      if (order_by == 'pvalue') {
        ord <- order(with(de_result, ifelse(significant, p, Inf)))
      } else if (order_by == 'log2fc') {
        ord <- order(-with(de_result, ifelse(significant, log2fc, -Inf)))
      }
      de_result[ord,]
    })
  } else {
    stop('Method: ', method, ' not recognized\n')
  }

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

#' Train a multinomial distribution from a vector
#'
#' @param x An input vector of values
#' @return A scaled vector of frequencies
#' @export
#' @examples
#' train_multinomial(1:10)
train_multinomial <- function(x) {
  (1+colSums(x)) / sum(1+colSums(x))
}

#' Merge clusters among cluster assignments
#'
#' @param cluster_result A vector of integer class labels for each cell
#' @param cluster_select A vector of cluster (integer) ids to select
#' @param colour (Optional) colours assignments for original clusters
#' @return A list containing: \itemize{
#' \item ids - new (consecutive) id labels for each cell after merging integer labels
#' \item col - new colours (mapping to the ids) corresponding to the integer labels if colour is specified, otherwise equals NULL
#' }
#' @export
#' @examples
#' merge_clusters_test(rep(1:5,4),c(1,3,4))
merge_clusters <- function(cluster_result,cluster_select,colour=NULL) {
  old_labs <- sort(unique(cluster_result))
  if (length(cluster_select) < 2 )
    stop('cluster_select should have more than two cluster IDs for mergeing')
  if (!is.integer(old_labs))
    stop('The cluster ids need to be integer types, e.g., outputs from clustering algorithms.')
  if (!(all(cluster_select %in% old_labs)))
    stop('The labels in {', cluster_select,'} include one or more invalid cluster IDs!\n')
  if (!is.null(colour)) {
    if (length(colour) != length(old_labs))
      stop('The length of colour is not the same as the number of unique labels!\n')
  }
  tmp_labs_map <- old_labs
  tmp_labs_map[old_labs %in% cluster_select] <- cluster_select[length(cluster_select)]
  uniq_tmp_labs <- sort(unique(tmp_labs_map))
  new_labs_map <- unlist(lapply(tmp_labs_map,function(x) which(uniq_tmp_labs %in% x)))
  names(new_labs_map) <- old_labs # map from old_labs to new_labs_map
  new_clusters <- unlist(lapply(cluster_result, function(x) new_labs_map[x] )) # new cluster assignments
  names(new_clusters) <- names(cluster_result)
  if (!is.null(colour)) { # output color assignment
    new_colour <- unlist(lapply(uniq_tmp_labs, function(x) colour[x] ))
    names(new_colour) <- 1:length(uniq_tmp_labs)
  } else {
    new_colour <- NULL
  }
  return(list(ids=new_clusters, col=new_colour))
}

#' Select cells according to cluster IDs
#'
#' @param cluster_result A vector of class labels for each cell
#' @param cluster_select A vector (or scalar) of cluster ids to select
#' @return Logical whether or not a cell is selected
#' @export
#' @examples
#' select_cells_by_clusterID(rep(1:5,4),c(1,3,4))
select_cells_by_clusterID <- function(cluster_result,cluster_select) {
  if (!(all(cluster_select %in% cluster_result))) {
    stop('The labels in {', cluster_select,'} include one or more invalid cluster IDs!\n')
  }
  cell_select_indices <- cluster_result %in% cluster_select
  return(cell_select_indices)
}

#' Select cells by clusters according to the fraction of cells that express at least one of the input gene markers
#'
#' @param gbm A GeneBCMatrix (recommended to be normalized by UMI counts and log-transformed)
#' @param cluster A vector of class labels for each cell
#' @param gene_symbols A list of gene symbols
#' @param frac Minimum fraction of cells in a cluster required for the cluster to be selected
#' @return Logical whether or not a cell is selected
#' @export
#' @examples
#' select_cells_by_clusterID(rep(1:5,4),c(1,3,4))
select_clusters_by_gene_markers <- function(gbm, gene_symbols, cluster, frac) {
  gbm_trunc <- trunc_gbm_by_genes(gbm, gene_symbols)
  uniq_clu <- sort(unique(cluster))
  # select the clusters where over "frac" cells express any of the genes listed above
  total_umi_trunc <- colSums(exprs(gbm_trunc))
  n_cells_with_genes <- unlist(lapply(uniq_clu, function(x) sum(total_umi_trunc[cluster==x] > 0)))
  n_cells_per_cluster <- unlist(lapply(uniq_clu, function(x) length(total_umi_trunc[cluster==x])))
  cell_frac <- n_cells_with_genes / n_cells_per_cluster
  cID_sel <- uniq_clu[which(cell_frac > frac)]
  cat("Selected populations with >",frac,"that express",gene_symbols,":",cID_sel," \n")
  return(select_cells_by_clusterID(cluster,cID_sel))
}

#' Compute the correlation score with a bulk population profile
#'
#' @param gbm_in A GeneBCMatrix (recommended to be normalized by UMI counts)
#' @param compare_exprs A matrix that includes expressions with one typical expression per row
#' @param compare_gene_ids The gene ids corresponding to the gene expression
#' @return A matrix of correlation scores (each row is a cell in gbm_in) based on the intersection of genes between gbm_in and compare_exprs
#' @export
#' @examples
#' \dontrun{
#'  get_correletion_scores(gbm, bulk_profile$avg_expr, bulk_profile$genes$id)
#' }
get_correletion_scores <- function(gbm_in, compare_exprs, compare_gene_ids, method='spearman') {
  use_gene_ids <- fData(gbm_in)$id
  sig_genes <- intersect(use_gene_ids, compare_gene_ids)
  in_X <- as.matrix(t(exprs(gbm_in)[which(use_gene_ids %in% sig_genes), ])) # predictors
  bulk_filt <- compare_exprs[ ,which(compare_gene_ids %in% sig_genes)]

  z <- lapply(1:nrow(bulk_filt), function(j) sapply(1:nrow(in_X),
                                 function(i) cor(in_X[i,], bulk_filt[j,], method=method)))
  z <- do.call(cbind, z)
  colnames(z) <- rownames(bulk_filt)
  z
}


#' Compute the normalized dispersion for each gene
#'
#' @param gbm_norm A GeneBCMatrix normalized by UMI counts (or output from normalize_barcode_sums_to_median())
#' @param num_bins The number of bins to calculate the dispersion norm
#' @return A data frame containing \itemize{
#' \item "mean" - The mean expression each gene
#' \item "var" - The variance of each gene
#' \item "dispersion" - The computed dispersion (var/mean) of each gene
#' \item "dispersion_norm" - The dispersion noramlized by binning
#' \item "id" - gene ids
#' \item "symbol" - gene symbols
#' }
#' @import Biobase
#' @export
#' @examples
#' \dontrun{
#'  get_gene_dispersion(gbm1)
#' }
get_gene_dispersion <- function(gbm_in, num_bins=20) {
  m <- t(exprs(gbm_in))
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,length.out=num_bins)),Inf)))

  var_by_bin <- with(df, aggregate(dispersion ~ mean_bin, FUN =  function(x) c( bin_median = median(x), bin_mad=mad(x) )))
  var_by_bin$bin_median <- var_by_bin$dispersion[,'bin_median']
  var_by_bin$bin_mad <- var_by_bin$dispersion[,'bin_mad']

  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df$id <- fData(gbm_in)$id
  df$symbol <- fData(gbm_in)$symbol
  df
}

#' Compute sSeq parameters
#'
#' @param x (cell x gene) count matrix
#' @param zeta_quantile Take this quantile of the method-of-moments dispersion estimates as the shrinkage target
compute_sseq_params <- function(x, zeta_quantile=0.995) {
  # Condition i, replicate j, gene g
  params <- list()
  N <- nrow(x)
  G <- ncol(x)

  # Estimate size factors
  grand_median <- median(rowSums(x))
  s_ij <- rowSums(x) / grand_median
  params$s_ij <- s_ij

  # Method-of-moment esimators for per-gene mean and variance (negative binomial)
  x_size_norm <- x / s_ij

  mu_g <- colMeans(x_size_norm)
  v_g <- apply(x_size_norm, 2, var)

  # Only use genes with non-zero variance
  use_g <- v_g > 0
  params$mu_g <- mu_g
  params$v_g <- v_g
  params$use_g <- use_g

  # Method-of-moment estimates of per-gene dispersions
  phi_mm_g <- pmax(0, (N*v_g - mu_g*sum(1/s_ij)) / (mu_g^2*sum(1/s_ij)))
  params$phi_mm_g <- phi_mm_g

  # Compute zeta_hat, the optimal target dispersion
  zeta_hat <- quantile(params$phi_mm_g, zeta_quantile, na.rm=T)
  params$zeta_hat <- zeta_hat

  # Compute delta, the shrinkage
  mean_phi_mm_g <- mean(phi_mm_g[use_g])
  delta <- (sum((phi_mm_g[use_g] - mean_phi_mm_g)^2)/(G-1)) / (sum((phi_mm_g[use_g] - zeta_hat)^2)/(G-2))
  params$delta <- delta

  # Compute the shrunken dispersion estimate
  phi_g <- rep(NA, G)
  phi_g[use_g] <- (1-delta)*phi_mm_g[use_g] + delta*zeta_hat
  params$phi_g <- phi_g

  params
}

#' Negative binomial exact test as in sSeq
#'
#' @param x_a Total count for a single gene in group A
#' @param x_b Total count for the same gene in group B
#' @param s_a Size factor for group A
#' @param s_b Size factor for group B
#' @param mu Common mean for gene
#' @param phi Common dispersion for gene
nb_exact_test <- function(x_a, x_b, s_a, s_b, mu, phi) {
  # Compute p-value for pairwise exact test with negative binomial distribution
  # Note: runtime is O(n) in the max(count)

  all_x_a <- seq(0, x_a+x_b, 1)
  all_x_b <- seq(x_a+x_b, 0, -1)

  .prob <- function(x, s) dnbinom(x, mu=s*mu, size=1/(phi/s))
  p_obs <- .prob(x_a, s_a) * .prob(x_b, s_b)
  p_all <- .prob(all_x_a, s_a) * .prob(all_x_b, s_b)

  # Probability that random value under null hypothesis is more extreme than observed
  sum(p_all[p_all <= p_obs]) / sum(p_all)
}

#' Compute differential expression using the sSeq method
#'
#' @param x (cell x gene) Count matrix
#' @param cond0 (integer) Indices of cells in group A
#' @param cond1 (integer) Indices of cells in group B
#' @param sseq_params Params precomputed from whole sample
sseq_differential_expression <- function(x, cond0, cond1, sseq_params, gene_ids, gene_symbols) {
  x_a <- x[cond0,,drop=F]
  x_b <- x[cond1,,drop=F]

  G <- ncol(x)

  # size factors
  s_a <- sum(sseq_params$s_ij[cond0])
  s_b <- sum(sseq_params$s_ij[cond1])

  # gene sums
  x_ga <- colSums(x_a)
  x_gb <- colSums(x_b)

  # compute p-values
  p_res <- rep(NA, G)
  p_res[sseq_params$use_g] <- sapply(which(sseq_params$use_g),
                                     function(g) nb_exact_test(x_ga[g], x_gb[g],
                                                               s_a, s_b,
                                                               sseq_params$mu_g[g], sseq_params$phi_g[g]))

  # multiple testing correction
  p_adj <- rep(1.0, G)
  p_adj[sseq_params$use_g] <- p.adjust(p_res[sseq_params$use_g], method='BH')

  # summarize results
  data.frame(gene_id=gene_ids, gene_name=gene_symbols,
             tested=sseq_params$use_g,
             sum_a=x_ga, sum_b=x_gb,
             common_mean=sseq_params$mu_g,
             dispersion=sseq_params$phi_g,
             mean_a_sizenorm=x_ga/s_a,
             mean_b_sizenorm=x_gb/s_b,
             # Introduce pseudocount to de-emphasize lower counts
             log2fc=log2((1+x_ga)/(1+s_a)) - log2((1+x_gb)/(1+s_b)),
             p=p_res,
             p_adj=p_adj)
}

#' Run differential expression analysis between two sets of cells
#'
#' @param gbm GeneBCMatrix
#' @param cell_indices0 (integer) Indices of cells in group A
#' @param cell_indices1 (integer) Indices of cells in group B
#' @return A data frame containing differential expression results for each gene.
#' @export
#' @examples
#' \dontrun{
#' de_result <- run_differential_expression(gbm, 1:500, 501:1000)
#' top_results <- de_result %>% filter(p_adj < 0.05) %>% arrange(desc(log2fc))
#' head(top_results)
#' }
run_differential_expression <- function(gbm, cell_indices0, cell_indices1) {
  # Cell x gene matrix
  mat <- t(exprs(gbm))

  # Compute sSeq params
  cat('Comparing parameters...\n')
  sseq_params <- compute_sseq_params(mat)

  # Run DE tests
  cat('Testing genes for differential expression...\n')
  sseq_differential_expression(mat, cell_indices0, cell_indices1,
                               sseq_params,
                               gene_ids=fData(gbm)$id,
                               gene_symbols=fData(gbm)$symbol)
}


