#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

#' Plotting functions
#'

#' Plot barcode distribution
#'
#' Plot the barcode counts in descending order on a log-scale
#'
#' @param gene_bc_matrix A GeneBCMatrix object
#' @return none. plot is generated.
#' @export
#' @examples
#' \dontrun{
#' plot_barcode_counts(gbm1)
#' }
#'
plot_barcode_counts <- function(gene_bc_matrix) {
  bc_counts <- sort(colSums(exprs(gene_bc_matrix)), decreasing=TRUE)
  nonzero_bc_counts <- bc_counts[bc_counts > 0]
  plot(nonzero_bc_counts ~ I(1:length(nonzero_bc_counts)), t='l', log='xy',
       xlab='Barcode rank', ylab='Barcode UMI counts')
}

#' Visualize total umi counts under a 2D projection
#'
#' Generate a ggplot object that highlights UMI counts
#'
#' @param gbm A GeneBCMatrix object (NOT log-transformed)
#' @param projection A two column matrix projection of each cell
#' @param limits (min,max) saturates values on the color bar
#' @return A ggplot object with facets corresponding to each gene symbol
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' visualize_umi_counts(gbm,tsne_proj)
#' }
#'
visualize_umi_counts<- function(gbm,projection,limits=c(0,10),marker_size=0.1) {
  gene_values <- log(colSums(exprs(gbm)),base=10)
  gene_values[gene_values<limits[1]] <- limits[1]
  gene_values[gene_values>limits[2]] <- limits[2]
  projection_names <-  colnames(projection)
  colnames(projection) <- c('Component.1', 'Component.2')
  proj_gene <- data.frame(cbind(projection,gene_values))
  p <- ggplot(proj_gene, aes(Component.1, Component.2)) +
    geom_point(aes(colour=gene_values),size=marker_size) +
    scale_colour_gradient(low="blue",high="red",name = "log10") +
    labs(x=projection_names[1],y=projection_names[2]) + ggtitle('Total number of UMIs') + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(p)
}

#' Visualize gene markers under a 2D projection
#'
#' Generate a ggplot object that highlights gene markers
#'
#' @param gbm A GeneBCMatrix object (recommended to be normalized by UMI counts and log-transformed)
#' @param gene_probes A vector of gene symbols to highlight
#' @param projection A two column matrix projection of each cell
#' @param limits (min,max) saturates values on the color bar
#' @return A ggplot object with facets corresponding to each gene symbol
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' visualize_gene_markers(gbm,c('gene1','gene2'),tsne_proj)
#' }
#'
visualize_gene_markers<- function(gbm,gene_probes,projection,limits=c(0,10),marker_size=0.1,title=NULL) {
  gbm_trunc <- trunc_gbm_by_genes(gbm,gene_probes)
  gene_values <- t(as.matrix(exprs(gbm_trunc)))
  gene_values[gene_values<limits[1]] <- limits[1]
  gene_values[gene_values>limits[2]] <- limits[2]
  colnames(gene_values) <- gene_probes
  projection_names <-  colnames(projection)
  colnames(projection) <- c('Component.1', 'Component.2')
  proj_gene <- data.frame(cbind(projection,gene_values))
  proj_gene_melt <- melt(proj_gene, id.vars=c("Component.1", "Component.2"))
  p<- ggplot(proj_gene_melt, aes(Component.1, Component.2)) +
    geom_point(aes(colour=value),size=marker_size) + facet_wrap(~variable) +
    scale_colour_gradient(low="grey",high="red",name = "val") +
    labs(x=projection_names[1],y=projection_names[2])
  if (!is.null(title)) {  p <- p + ggtitle(title) }
  p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(p)
}

#' Visualize class labels under a 2D projection
#'
#' Generate a ggplot object that highlights specified class labels
#'
#' @param cluster_result A vector or matrix of class labels for each cell (row)
#' @param projection A two column matrix projection of each cell
#' @param colour Customized colors for each class label
#' @return A ggplot object with facets corresponding to each gene symbol
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' visualize_clusters(cluster_result,tsne_proj)
#' }
#'
visualize_clusters <-function(cluster_result,projection,colour=NULL,alpha=1,marker_size=0.1,title=NULL) {
  projection_names <-  colnames(projection)
  colnames(projection) <- c('Component.1', 'Component.2')
  proj_clu <- data.frame(cbind(projection,cluster_result))
  proj_clu_melt <- melt(proj_clu, id.vars=c('Component.1','Component.2'))
  p <- ggplot(proj_clu_melt, aes(Component.1, Component.2))+
    geom_point(aes(colour = factor(cluster_result)),size=marker_size,alpha=alpha)+
    guides(col = guide_legend(title="ID",override.aes = list(size=3)))+facet_wrap(~variable)+
    labs(x=projection_names[1],y=projection_names[2]) + theme_bw()
  if (!is.null(title)) {  p <- p + ggtitle(title) }
  if ((!is.null(colour))) {
    names(colour) <- 1:length(colour)
    p <- p + scale_color_manual(values=colour)
  }

  if (is.vector(cluster_result)) {
    p <- p + theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  } else {
    p <- p + theme(plot.title = element_text(hjust = 0.5),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }

  return(p)
}

#' Plot heatmap with specified gene and cell ordering
#'
#' Generate a ggplot object that highlights specified class labels
#'
#' @param gbm A GeneBCMatrix object (recommended to be normalized by UMI counts and log-transformed)
#' @param genes_to_plot A list of gene symbols (or output from function: prioritize_genes)
#' @param cells_to_plot A dataframe with cell names and cell indices (or output from function: order_cell_by_clusters)
#' @param n_genes Number of genes to include in the heatmap
#' @param colour Customized colors for each class label on the colour axes
#' @param limits (min,max) saturates values on the color bar
#' @return none. a pheatmap is generated.
#' @import pheatmap
#' @export
#' @examples
#' \dontrun{
#' gbm_pheatmap(gbm,genes_to_plot,cells_to_plot)
#' }
#'
gbm_pheatmap <- function(gbm,genes_to_plot,cells_to_plot,n_genes=5,colour=NULL,limits=c(-3,3)) {
  if (!is.list(genes_to_plot)) { # custmoized gene list
    cat('Plotting one gene set instead of multiple cluster-specific gene sets\n')
    gene_indices <- sapply(genes_to_plot, function(x) get_gene_index(gbm,x))
    gene_annotation <- NULL
  } else {
    gene_indices <- unlist(lapply(genes_to_plot,function(x) x$ix[1:n_genes]))
    gene_grouping <- rep(names(genes_to_plot),each=n_genes)
    gene_annotation <- data.frame(ClusterID = as.factor(gene_grouping))
  }
  cell_indices <- unlist(lapply(cells_to_plot,function(x) x$ix))
  # scaling by genes and saturating the gene expression values and saturate at limits
  value <- t(scale(t(as.matrix(exprs(gbm))[gene_indices,cell_indices])))
  value[value<limits[1]] <- limits[1]
  value[value>limits[2]] <- limits[2]
  rownames(value) <- make.unique(fData(gbm)$symbol[gene_indices])
  cell_grouping <- unlist(lapply(1:length(cells_to_plot), function(x) {
                          rep(names(cells_to_plot)[x], length(cells_to_plot[[x]]$barcode))}))
  # create annotation of x and y axes
  cell_annotation <- data.frame(ClusterID = as.factor(cell_grouping))
  rownames(cell_annotation) <- colnames(value)
  if (!is.null(gene_annotation)) { rownames(gene_annotation) <- rownames(value) }
  if (is.null(colour)) {
    anno_colors <- NULL
  } else {
    names(colour) <- names(cells_to_plot)
    anno_colors <- list(ClusterID = colour)
  }
  pheatmap(value, cluster_rows = FALSE,cluster_cols = FALSE, show_colnames = F,
           annotation_row=gene_annotation,annotation_col=cell_annotation,
           annotation_names_row = FALSE, annotation_names_col = FALSE,
           annotation_colors = anno_colors)
}

