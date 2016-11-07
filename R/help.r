#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

#' Analyze data from 10x Genomics Single Cell pipeline
#'
#' Load, manipulate, analyze, and plot data from the 10x Genomics Single Cell pipeline, Cell Ranger
#'
#' For a full example, run \code{vignette('rkit-vignette')}
#'
#' The following functions are provided by this package to assist users of 10x Genomics Single Cell
#' in secondary analysis of their data.
#'
#' @section Loading data from Cell Ranger's output:
#' \itemize{
#' \item{\code{\link{load_cellranger_matrix}}}
#' \item{\code{\link{load_cellranger_analysis_results}}}
#' }
#'
#' @section Plotting data from Cell Ranger's output:
#' \itemize{
#' \item{\code{\link{plot_barcode_counts}}}
#' }
#'
#' @section Analyzing data from Cell Ranger's output:
#' \itemize{
#' \item{\code{\link{run_pca}}}
#' \item{\code{\link{run_tsne}}}
#' \item{\code{\link{run_kmeans_clustering}}}
#' }
#'
#' @section Getting info from Cell Ranger's output:
#' \itemize{
#' \item{\code{\link{get_mean_mapped_reads_per_cell}}}
#' \item{\code{\link{get_mean_raw_reads_per_cell}}}
#' }
#'
#' @section Normalizing samples to make them comparable:
#' \itemize{
#' \item{\code{\link{load_molecule_info}}}
#' \item{\code{\link{subsample_gene_bc_matrix}}}
#' }
#'
#' @section Combining data from multiple samples:
#' \itemize{
#' \item{\code{\link{concatenate_gene_bc_matrices}}}
#' }
#'
#' @author 10x Genomics <software@@10xgenomics.com>
#' @docType package
#' @name cellrangerRkit
#'
NULL

