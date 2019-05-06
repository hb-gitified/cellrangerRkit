#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

#' Input/Output functions
#'
#' @include constants.r

#' Ugly but apparently recommended way to avoid data.table
#' calls from giving no visible binding to global variable errors
#' @name bind variables to fix data.table errors when packaging
globalVariables(c("reads", "barcode", "gene", "gem_group"))

#' Check if a directory exists
#'
#' @param path_name Directory name
#' @param stop (Optional) stop if the file does not exist
#' @return Logical TRUE if exists and FALSE otherwise
directory.exists <- function(dir_name, stop_message = NULL) {
  if (is.na(file.info(dir_name)$isdir)) {
    exist <- FALSE
  } else {
    exist <- file.info(dir_name)$isdir
  }
  if (!is.null(stop_message)) {
    if (!exist) { stop(stop_message) }
  }
  return(exist)
}

#' Get path to matrices in Cell Ranger output
#'
#' @param pipestance_path Path to the output directory produced by the Cell Ranger pipeline
#' @param barcode_filtered If true, get the path to the barcode filtered matrices; else to the raw filtered matrices
#' @return The path to the matrices directory
get_matrix_dir_path <- function(pipestance_path, barcode_filtered) {
  fname <- file.path(pipestance_path, 'outs', ifelse(barcode_filtered, 'filtered_gene_bc_matrices', 'raw_gene_bc_matrices'))
  if (directory.exists(fname)) {
    return(fname)
  }

  # cellranger 1.2.0 changed to *_mex
  fname <- file.path(pipestance_path, 'outs', ifelse(barcode_filtered, 'filtered_gene_bc_matrices_mex', 'raw_gene_bc_matrices_mex'))
  if (directory.exists(fname)) {
    return(fname)
  }

  stop(sprintf("Could not find matrix folder: %s\n", fname))
}

#' Get summary csv in Cell Ranger output
#'
#' @param pipestance_path Path to the output directory produced by the Cell Ranger pipeline
#' @return Path to the summary CSV
get_summary_csv_path <- function(pipestance_path) {
  file.path(pipestance_path, 'outs', 'metrics_summary.csv')
}

#' Get list of genomes in Cell Ranger output
#'
#' @param matrix_path Path to the output directory produced by the Cell Ranger pipeline
#' @param barcode_filtered Path to filtered matrix if true, otherwise raw matrix
#' @return A character vector of the genomes found
get_genome_in_matrix_path <- function(matrix_path, genome=NULL) {
  cat("Searching for genomes in:",matrix_path,"\n")
  genomes <- dir(matrix_path)
  if (is.null(genome)) {
    if (length(genomes) == 1) {
      genome <- genomes[1]
    } else {
      stop(sprintf("Multiple genomes found; please specify one. \n Genomes present: %s",paste(genomes, collapse=", ")))
    }
  } else if (!(genome %in% genomes)) {
    stop(sprintf("Could not find specified genome: '%s'. Genomes present: %s",
                 genome,paste(genomes, collapse=", ")))
  }
  cat("Using",genome,"in folder:",file.path(matrix_path,genome),"\n")
  return(genome)
}


#' Load data from the 10x Genomics Cell Ranger pipeline by specifying individual file paths
#'
#' @param mat_fn Path to a .mtx file including matirx entries
#' @param gene_fn  Path to a .tsv file including gene information
#' @param barcode_fn Path to a .tsv file including cell barcode information
#' @param summary_fn (Optional) Path to a .csv file including cell barcode information
#' @return A GeneBCMatrix sparse matrix where rows are genes and columns are cell-barcodes
#' @export
#' @import Matrix
#' @examples
#' \dontrun{
#' gbm <- load_cellranger_matrix_from_files("matrix.mtx","genes.tsv","barcodes.tsv")
#' }
load_cellranger_matrix_from_files <- function(mat_fn, gene_fn, barcode_fn, summary_fn=NULL) {

  # Load the matrix file (Mandatory for gbm object)
  if (file.exists(mat_fn)) {
    mat <- readMM(mat_fn)
    cat('Loaded matrix information\n')
  } else {
    stop(sprintf("Could not find matrix file \n\t%s\n", mat_fn))
  }

  # Load the gene file (Mandatory for gbm object)
  if (file.exists(gene_fn)) {
    gene_info <- read.delim(gene_fn, stringsAsFactors=FALSE, sep="\t", header=FALSE)
    if (dim(mat)[1] != length(gene_info[,1])) {
      stop(sprintf("Mismatch dimension between gene file: \n\t %s\n and matrix file: \n\t %s\n", gene_fn,mat_fn))
    } else {
      rownames(mat) <- gene_info[,1]
      gene_symbols <- gene_info
      row.names(gene_symbols) = gene_symbols[, 1]
      colnames(gene_symbols) = c("id", "symbol")
      cat('Loaded gene information\n')
    }
  } else {
    stop(sprintf("Could not find gene file: \n\t %s\n", gene_fn))
  }

  # Load the barcode file (Mandatory for gbm object)
  if (file.exists(barcode_fn)) {
    barcodes <- read.delim(barcode_fn, stringsAsFactors=FALSE, sep="\t", header=FALSE)
    if (dim(mat)[2] != length(barcodes[,1])) {
      stop(sprintf("Mismatch dimension between barcode file: \n\t %s\n and matrix file: \n\t %s\n", barcode_fn,mat_fn))
    } else {
      colnames(mat) <- barcodes[,1]
      pd = data.frame(id=barcodes[,1], row.names=barcodes[,1])
      colnames(pd) = c("barcode")
      cat('Loaded barcode information\n')
    }
  } else {
    stop(sprintf("Could not find barcode file: \n\t %s\n", barcode_fn))
  }

  # Build GeneBCMatrix object
  res <- newGeneBCMatrix(mat=mat, fd=gene_symbols, pd=pd)
  res@subsampled <- FALSE

  # Load the summary csv file (Optional)
  if (is.null(summary_fn)) {
    warning("No summary file provided. Some functions may be disabled without the metrics in summary csv.\n")
  } else {
    if (file.exists(summary_fn)) {
      summary <- read.csv(summary_fn, as.is=TRUE, header=TRUE)
      res@summary <- summary
      cat('Loaded summary information\n')
    } else {
      cat(sprintf("Could not find summary csv: \n\t %s.\nThis file is only necessary if you are performing depth-normalization (calling the equalize_gbms function) in R.\nIf this pipestance was produced by `cellranger aggr` with the default parameters, depth-normalization in R (via equalize_gbms) is not necessary.\n",summary_fn))
    }
  }

  # res@barcode_filtered <- barcode_filtered
  # if (!is.null(pipestance_path)) {
  #   res@pipestance_path <- pipestance_path
  # }
  return(res)
}

#' Load data from the 10x Genomics Cell Ranger pipeline
#'
#' @param pipestance_path Path to the output directory produced by Cell Ranger
#' @param genome The desired genome (e.g., 'hg19' or 'mm10')
#' @param barcode_filtered Load only the cell-containing barcodes
#' @return A GeneBCMatrix sparse matrix where rows are genes and columns are cell-barcodes
#' @export
#' @import Matrix
#' @examples
#' \dontrun{
#' # Load from a Cell Ranger output directory
#' gene_bc_matrix <- load_cellranger_matrix("/home/user/cellranger_output")
#' }
load_cellranger_matrix <- function(pipestance_path=NULL, genome=NULL, barcode_filtered=TRUE) {
  # check for correct directory structure
  if (!directory.exists(pipestance_path))
    stop("Could not find the pipestance path: '", pipestance_path,"'. Please double-check if the directory exists.\n")
  if (!directory.exists(file.path(pipestance_path,'outs')))
    stop("Could not find the pipestance output directory: '", file.path(pipestance_path,'outs'),"'. Please double-check if the directory exists.\n")

  matrix_path <- get_matrix_dir_path(pipestance_path, barcode_filtered)
  genome <- get_genome_in_matrix_path(matrix_path, genome=genome)

  mat_fn <- file.path(matrix_path, genome, "matrix.mtx")
  gene_fn <- file.path(matrix_path, genome, "genes.tsv")
  barcode_fn <- file.path(matrix_path, genome, "barcodes.tsv")
  summary_fn <- get_summary_csv_path(pipestance_path)

  res <- load_cellranger_matrix_from_files(mat_fn,gene_fn,barcode_fn,summary_fn)
  res@barcode_filtered <- barcode_filtered
  res@subsampled <- FALSE
  res@pipestance_path <- pipestance_path

  return(res)
}

#' Load gene-barcode matrix from a specific H5 file
#'
#' @param Path to the gene-barcode matrix H5 file (NOT the pipestance path)
#' @param genome The desired genome (e.g., 'hg19' or 'mm10')
#' @return A GeneBCMatrix sparse matrix where rows are genes and columns are cell-barcodes
#' @export
#' @import Matrix
#' @import rhdf5
#' @examples
#' \dontrun{
#' # Load hg19 from a filtered matrix h5
#' gene_bc_matrix <- get_matrix_from_h5("/home/user/cellranger_output/outs/filtered_gene_bc_matrices_h5.h5", "hg19")
#' }
get_matrix_from_h5 <- function(filename, genome=NULL) {
  if(!file.exists(filename)) {
    stop(sprintf("Could not find matrix H5 file: \n\t %s\n", filename))
  }

  genomes <- attributes(h5dump(filename, load=FALSE))$names

  if (is.null(genome)) {
    if (length(genomes) > 1) {
      stop(sprintf("Multiple genomes found; please specify one. \n Genomes present: %s",paste(genomes, collapse=", ")))
    }
    genome <- genomes[1]
  }

  if (!(genome %in% genomes)) {
    stop(sprintf("Genome %s not found in file. \n Genomes present: %s",genome,paste(genomes, collapse=", ")))
  }

  dset <- h5read(filename, name = genome)
  sparse_mat <- sparseMatrix(i = dset$indices + 1, p = dset$indptr, x = as.numeric(dset$data), dims = dset$shape, giveCsparse = FALSE)
  genes <- data.frame(id = dset$genes, symbol = dset$gene_names, row.names = dset$genes)
  barcodes <- data.frame(barcode = dset$barcodes, row.names = dset$barcodes)
  gbm <- newGeneBCMatrix(sparse_mat, genes, barcodes)
  H5close()
  return (gbm)
}

#' Load HDF5 data from the 10x Genomics Cell Ranger pipeline
#'
#' @param pipestance_path Path to the output directory produced by Cell Ranger
#' @param genome The desired genome (e.g., 'hg19' or 'mm10')
#' @param barcode_filtered Load only the cell-containing barcodes
#' @return A GeneBCMatrix sparse matrix where rows are genes and columns are cell-barcodes
#' @export
#' @import Matrix
#' @import rhdf5
#' @examples
#' \dontrun{
#' # Load from a Cell Ranger output directory
#' gene_bc_matrix <- load_cellranger_matrix_h5("/home/user/cellranger_output")
#' }
load_cellranger_matrix_h5 <- function(pipestance_path=NULL, genome=NULL, barcode_filtered=TRUE) {
  # check for correct directory structure
  if (!directory.exists(pipestance_path))
    stop("Could not find the pipestance path: '", pipestance_path,"'. Please double-check if the directory exists.\n")
  if (!directory.exists(file.path(pipestance_path,'outs')))
    stop("Could not find the pipestance output directory: '", file.path(pipestance_path,'outs'),"'. Please double-check if the directory exists.\n")

  h5_path <- file.path(pipestance_path, 'outs', ifelse(barcode_filtered, 'filtered_gene_bc_matrices_h5.h5', 'raw_gene_bc_matrices_h5.h5'))
  if (!file.exists(h5_path)) {
    stop(sprintf("Could not find matrix H5: %s\n", h5_path))
  }

  res <- get_matrix_from_h5(h5_path, genome)
  res@barcode_filtered <- barcode_filtered
  res@subsampled <- FALSE
  res@pipestance_path <- pipestance_path

  return (res)
}

#' Load Cell Ranger secondary analysis results
#'
#' @param pipestance_path Path to the output directory produced by Cell Ranger
#' @return A list containing: \itemize{
#' \item pca - pca projection
#' \item tsne - tsne projection
#' \item kmeans - kmeans results for various values of K
#' }
#' @export
#' @examples
#' \dontrun{
#' analysis <- load_cellranger_analysis_results("/home/user/cellranger_output")
#' }
load_cellranger_analysis_results <- function(pipestance_path) {
  analysis_path <- file.path(pipestance_path,'outs','analysis')
  if (directory.exists(analysis_path, stop_message=paste("Could not find:", analysis_path))) {
    return(load_analysis_results_from_path(analysis_path))
  }
}

#' Load Cell Ranger PCA results from a specified analysis folder
#'
#' @param analysis_path Path to the analysis output directory produced by Cell Ranger
#' @return Path to PCA results
#' @examples
#' \dontrun{
#' pca_path <- get_pca_path(analysis_path)
#' }
get_pca_path <- function(analysis_path) {
  pca_dir <- file.path(analysis_path, 'pca')

  # Prior to Cell Ranger 1.2.0
  if (file.exists(file.path(pca_dir, 'projection.csv'))) {
    return(pca_dir)
  }

  # Cell Ranger 1.2.0+
  pca_dir <- file.path(pca_dir, dir(pca_dir)[1])
  if (directory.exists(pca_dir)) {
    return(pca_dir)
  }
  return('')
}

#' Load Cell Ranger clustering results from a specified analysis folder
#'
#' @param analysis_path Path to the analysis output directory produced by Cell Ranger
#' @return Path to clustering results
#' @examples
#' \dontrun{
#' pca_path <- get_pca_path(analysis_path)
#' }
get_clustering_path <- function(analysis_path) {
  # Prior to Cell Ranger 1.3.0
  clustering_dir <- file.path(analysis_path, 'kmeans')
  if (directory.exists(clustering_dir)) {
    return(clustering_dir)
  }

  # Cell Ranger 1.3.0+
  return(file.path(analysis_path, 'clustering'))
}

#' Load Cell Ranger t-SNE results from a specified analysis folder
#'
#' @param analysis_path Path to the analysis output directory produced by Cell Ranger
#' @return Path to t-SNE results
#' @examples
#' \dontrun{
#' tsne_path <- get_tsne_path(analysis_path)
#' }
get_tsne_path <- function(analysis_path) {
  tsne_dir <- file.path(analysis_path, 'tsne')

  # Prior to Cell Ranger 1.2.0
  if (file.exists(file.path(tsne_dir, 'projection.csv'))) {
    return(tsne_dir)
  }

  # Cell Ranger 1.2.0+
  tsne_dir <- file.path(tsne_dir, '2_components')
  if (directory.exists(tsne_dir)) {
    return(tsne_dir)
  }

  return('')
}

#' Load Cell Ranger secondary analysis results from a specified analysis folder
#'
#' @param analysis_path Path to the analysis output directory produced by Cell Ranger
#' @return A list containing: \itemize{
#' \item pca - pca projection
#' \item tsne - tsne projection
#' \item clustering - clustering results for various clustering methods
#' }
#' @export
#' @examples
#' \dontrun{
#' analysis <- load_analysis_results_from_path("/home/user/cellranger_output/outs/analysis")
#' }
load_analysis_results_from_path <- function(analysis_path) {
  clustering_dir <- get_clustering_path(analysis_path)
  pca_dir <- get_pca_path(analysis_path)
  tsne_dir <- get_tsne_path(analysis_path)

  if (directory.exists(clustering_dir)) {
    if (basename(clustering_dir) == 'kmeans') {
      # Pre CR-1.3.0 there was just one clustering type called "kmeans."
      # Append the string "kmeans_" to the clusterings to harmonize them with CR-1.3.0+
      clustering_paths <- dir(clustering_dir)
      clustering_names <- paste0('kmeans_', clustering_paths)
    } else {
      clustering_names <- clustering_paths <- dir(clustering_dir)
    }
    clustering <- lapply(clustering_paths, function(x) read.csv(file.path(clustering_dir, x, 'clusters.csv')))
    names(clustering) <- clustering_names
  } else {
    clustering <- NULL
    warning("Clustering results not found in ", clustering_dir)
  }

  if (directory.exists(pca_dir)) {
    pca_results <- read.csv(file.path(pca_dir, 'projection.csv'))
  } else {
    pca_results <- NULL
    warning("PCA results not found in ", pca_dir)
  }

  if (directory.exists(tsne_dir)) {
    tsne_results <- read.csv(file.path(tsne_dir, 'projection.csv'))
  } else {
    tsne_results <- NULL
    warning("tSNE results not found in ", tsne_dir)
  }

  return(list(pca=pca_results, tsne=tsne_results, clustering=clustering))
}

#' Load molecule info h5 file from the Cell Ranger pipeline
#'
#' @param gbm GeneBCMatrix object to load into
#'
#' @export
#' @return A GeneBCMatrix object containing the loaded molecule info
#' @examples
#' \dontrun{
#' gbm1 = load_molecule_info(gbm1)
#' # Note molecule info stored as gbm1@molecule_info
#' }
load_molecule_info <- function(gbm) {
  warning("Loading the molecule info is only necessary if you are subsampling reads. This method of normalizing matrices is deprecated.\nPlease use the `cellranger aggr` pipeline (new in cellranger 1.2.0), which can combine arbitrary gene-barcode matrices\nand produce a combined, depth-normalized matrix.")
  if (identical(gbm@pipestance_path, character(0))) {
    stop("GeneBCMatrix does not have pipestance_path. Use 'load_molecule_info_from_path(gbm, f_path))' instead")
  } else {
    f_path <- file.path(gbm@pipestance_path, "outs", "molecule_info.h5")
    if (!file.exists(f_path)) {
      warning(sprintf("Could not find %s. Checking for 1.2.0 aggr file raw_molecule_info.h5", f_path))
      # Cell Ranger v1.2.0 `aggr` output
      f_path <- file.path(gbm@pipestance_path, "outs", "raw_molecule_info.h5")
    }
  }
  return(load_molecule_info_from_path(gbm, f_path))
}

#' Load molecule info h5 file from a specified .h5 file
#'
#' @param gbm GeneBCMatrix object to load into
#' @param f_path Path to the molecule_info.h5
#'
#' @import rhdf5
#' @import data.table
#' @export
#' @return A GeneBCMatrix object containing the loaded molecule info
#' @examples
#' \dontrun{
#' gbm1 = load_molecule_info(gbm1, f_path)
#' # Note molecule info stored as gbm1@molecule_info
#' }
load_molecule_info_from_path <- function(gbm, f_path) {

  if(! "bit64" %in% installed.packages()) # check package
    stop("Suggested package bit64 is required to use load_molecule_info. Please install and try again.")

  if (!file.exists(f_path)) # check if the h5 file exists
    stop('File does not exist: ', f_path)

  required_cols <- c("barcode", "gem_group", "gene", "umi", "reads")
  cols <- h5ls(f_path)$name
  for (col in required_cols) {
    if (!(col %in% cols)) {
      stop(sprintf("Missing '%s' column in molecule info h5", col))
    }
  }
  col_data <- lapply(required_cols, function(col) {
    h5read(f_path, col, bit64conversion='bit64')
  })
  names(col_data) <- required_cols

  dt <- as.data.table(data.frame(col_data, stringsAsFactors=FALSE))
  dt <- dt[reads > 0]

  # Aggregate reads per-molecule
  setkeyv(dt, c("barcode", "gem_group", "gene", "umi"))
  dt <- dt[, j=list(reads=sum(reads)), by=c('barcode', 'gem_group', 'gene', 'umi')]

  gbm@molecule_info <- dt
  cat('Loaded molecule information\n')
  return(gbm)
}

#' Download file to local and check if it exists
#'
#' @param url_name File URL
#' @return loc_name Local file path where the file is stored
download_file_and_check <- function(url_name, loc_name) {
  cat('Downloading:',url_name,'\n')
  download.file(url_name, destfile = loc_name, quiet = TRUE)
  if (!file.exists(loc_name)) {
    stop(paste('Failed to download to :',loc_name))
  } else {
    cat(paste('Successfully downloaded to :',loc_name),'\n')
  }
}

#' Download a sample and relevant information from software.10xgenomics.com
#'
#' NOTE: This pipestance format is for v1.1 and earlier
#'
#' @param sample_name A name of the sample to be downloaded
#' @param sample_dir Directory where the data is saved
#' @param host The prefix of the url where the data is downloaded from
#' @param lite (Optional, default = TRUE) If TRUE, then molecule_info will not be downloaded
#' @param barcode_unfiltered (Optional, default = FALSE) If TRUE, then raw_gene_bc_matrix will be downloaded
#' @return None The sample_dir can be used to create a gbm object
#' @export
#' @examples
#' \dontrun{
#'  download_sample('pbmc3k','/home/Downloads','https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/')
#'  }
download_sample <- function(sample_name, sample_dir, host, lite=TRUE, barcode_unfiltered=FALSE) {
  # NOTE: This pipestance format is for v1.1 and earlier

  # header to download the data
  header <- paste(host,sample_name,"/",sample_name,"_",sep="")

  # create sample directory
  sample_dir <- path.expand(sample_dir)
  dir.create(sample_dir,showWarnings = FALSE)
  # create outs directory
  outs_dir <- file.path(sample_dir,"outs")
  if (directory.exists(outs_dir)) { unlink(outs_dir, recursive = TRUE) } # remove outs directory if already exists
  dir.create(outs_dir,showWarnings = FALSE)
  directory.exists(outs_dir, stop_message = paste('Failed to create outs directory in ',sample_dir))
  # --------------------------------------
  # download csv file
  # e.g., "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_metrics_summary.csv"
  download_file_and_check(paste(header,"metrics_summary.csv",sep=""), file.path(outs_dir,"metrics_summary.csv"))
  # --------------------------------------
  # download gene barcode matrix
  # e.g., "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
  filt_gbm_tar_path <- file.path(outs_dir,"filtered_gene_bc_matrices.tar.gz")
  download_file_and_check(paste(header,"filtered_gene_bc_matrices.tar.gz",sep=""),filt_gbm_tar_path)
  untar(filt_gbm_tar_path, exdir = outs_dir, verbose=TRUE)
  filt_gbm_dir <- file.path(outs_dir,"filtered_gene_bc_matrices")
  if (directory.exists(file.path(outs_dir,"filtered_matrices_mex"))) {  # handle alias
    file.rename(file.path(outs_dir,"filtered_matrices_mex"), filt_gbm_dir)
  }
  if (directory.exists(filt_gbm_dir, stop_message=paste('Failed to untar:',filt_gbm_tar_path))) {
    cat('Successfully untared to: ',filt_gbm_tar_path,'\n')
  }
  # --------------------------------------
  # download raw gene barcode matrix (optional)
  # e.g., "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_raw_gene_bc_matrices.tar.gz"
  if (barcode_unfiltered) {
    raw_gbm_tar_path <- file.path(outs_dir,"raw_gene_bc_matrices.tar.gz")
    download_file_and_check(paste(header,"raw_gene_bc_matrices.tar.gz",sep=""),raw_gbm_tar_path)
    untar(raw_gbm_tar_path, exdir = outs_dir, verbose=TRUE)
    raw_gbm_dir <- file.path(outs_dir,"raw_gene_bc_matrices")
    if (directory.exists(file.path(outs_dir,"matrices_mex"))) {  # handle alias
      file.rename(file.path(outs_dir,"matrices_mex"), raw_gbm_dir)
    }
    if (directory.exists(raw_gbm_dir, stop_message=paste('Failed to untar:',raw_gbm_tar_path))) {
      cat('Successfully untared to: ',raw_gbm_tar_path,'\n')
    }
  }
  # --------------------------------------
  # download cluster analysis file
  # e.g., "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_analysis.tar.gz"
  analysis_tar_path <- file.path(outs_dir,"analysis.tar.gz")
  download_file_and_check(paste(header,"analysis.tar.gz",sep=""),analysis_tar_path)
  untar(analysis_tar_path, exdir = outs_dir, verbose=TRUE)
  analysis_path <- file.path(outs_dir,"analysis")
  if (directory.exists(file.path(outs_dir,"analysis_csv"))) { # handle alias
    file.rename(file.path(outs_dir,"analysis_csv"), analysis_path)
  }
  if (directory.exists(analysis_path, stop_message=paste('Failed to untar:',analysis_tar_path))) {
    cat('Successfully untared to: ',analysis_tar_path,'\n')
  }
  # --------------------------------------
  # download molecule information (optional if no sample merging is required) *warning* this is slow!
  # e.g., "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_molecule_info.h5"
  if (!lite) {
    download_file_and_check(paste(header,"molecule_info.h5",sep=""), file.path(outs_dir,"molecule_info.h5"))
  } else {
    cat("molecule_info.h5 was not downloaded because 'lite=TRUE'. Note option 'lite=FALSE' is necessary if you want to merge gbms.\n")
  }
  return(cat('[DONE] Successfully downloaded all data to:',sample_dir,'\n'))
}

#' Save a single-genome GeneBCMatrix as an h5 file for input into `cellranger reanalyze`
#'
#' @param gbm GeneBCMatrix to write
#' @param out_filename The HDF5 filename to write to
#' @param genome_name The original genome name (e.g., GRCh38)
#' @param template Path to the original matrix.h5 file.
#' @param overwrite Overwrite the HDF5 file if it exists
#' @export
#' @import Matrix
#' @import rhdf5
#' @examples
#' \dontrun{
#' # Save a gbm as a matrix HDF5 file
#' save_cellranger_matrix_h5(gbm, "new_matrix.h5", "GRCh38")
#' }
save_cellranger_matrix_h5 <- function(gbm, out_filename, genome_name, template=NULL, overwrite=FALSE) {
  H5close()
  mat <- as(exprs(gbm), "dgCMatrix")

  if (file.exists(out_filename) && overwrite) {
    unlink(out_filename)
  } else if (file.exists(out_filename)) {
    cat(sprintf("Output file %s already exists; provide overwrite=TRUE to overwrite it.\n", out_filename))
    return()
  }

  h5createFile(out_filename)
  root <- H5Fopen(out_filename, flags= "H5F_ACC_RDWR")
  if (!is.null(template)) {
    library_ids <- h5readAttributes(template, '/')$library_ids
    ggs <- h5readAttributes(template, '/')$original_gem_groups
    h5writeAttribute(library_ids, root, 'library_ids')
    h5writeAttribute(ggs, root, 'original_gem_groups')
  }
  H5Fclose(root)
  h5createGroup(out_filename, genome_name)

  suppressWarnings({
    h5write(as.integer(mat@x), out_filename, sprintf('%s/data', genome_name))
    h5write(mat@i, out_filename, sprintf('%s/indices', genome_name))
    h5write(mat@p, out_filename, sprintf('%s/indptr', genome_name))
    h5write(c(nrow(mat), ncol(mat)), out_filename, sprintf('%s/shape', genome_name))
    h5write(as.character(fData(gbm)$id), out_filename, sprintf('%s/genes', genome_name))
    h5write(as.character(fData(gbm)$symbol), out_filename, sprintf('%s/gene_names', genome_name))
    h5write(as.character(pData(gbm)$barcode), out_filename, sprintf('%s/barcodes', genome_name))
  })
}
