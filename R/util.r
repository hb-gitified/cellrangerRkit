#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

#' Utility functions
#'
#' @include constants.r

#' Ugly but apparently recommended way to avoid data.table
#' calls from giving no visible binding to global variable errors
#' @name bind variables to fix data.table errors when packaging
globalVariables(c("reads", "barcode", "gene", "gem_group"))

#' Subsample a matrix to a given sequencing depth
#'
#' Subsample a matrix to a given sequencing depth so it is comparable to matrices from other samples
#'
#' @param gbm GeneBCMatrix object to subsample
#' @param target_reads_per_cell Target depth specified as reads per cell
#' @param target_type Type of reads specified in target_reads_per_cell; one of: \itemize{
#' \item \code{"raw_reads"} - Raw reads per cell.
#' \item \code{"mapped_reads"} Confidently mapped, barcoded reads per cell. Raw reads are sufficient in most cases, but if the valid barcode fraction or confidently mapped read fraction vary greatly between samples, this option should be used.
#' }
#' @param verbose Print status messages along the way
#' @export
#' @return A new GeneBCMatrix object that has been subsampled to the specified read depth
#' @examples
#' \dontrun{
#' new_reads_per_cell = 10000 # downsampling read depth per cell to 10000
#' gbm1 = load_molecule_info(gbm1)
#' gbm1_subsampled = subsample_gene_bc_matrix(gbm1, 10000)
#' }
subsample_gene_bc_matrix <- function(gbm, target_reads_per_cell,
                                     target_type="raw_reads", verbose=TRUE) {
  if (gbm@subsampled) {
    stop("Can't re-subsample a matrix - please use the original matrix to subsample from")
  }
  if (!(target_type %in% c("raw_reads", "mapped_reads"))) {
    stop(sprintf("Unsupported target_type: %s", target_type))
  }
  if (length(gbm@molecule_info) < 1) {
    stop("Molecule info not loaded. Please use load_molecule_info() before subsampling")
  }

  n_cells <- get_metric(gbm, CELL_COUNT_METRIC)
  if (verbose) {
    cat(sprintf("Cells: %d\n", n_cells))
  }

  tot_reads <- get_metric(gbm, TOTAL_READS_METRIC)
  if (verbose) {
    cat(sprintf("Reads: %d\n", tot_reads))
  }

  # Compress original barcodes and store as "N-G"
  # Where N is a compressed 2-bit-rep integer and G is the gem group
  original_bc_seqs <- get_barcode_sequence(colnames(exprs(gbm)))
  original_bc_gem_groups <- get_gem_group(colnames(exprs(gbm)))
  compressed_original_bc_seqs <- compress_sequences(original_bc_seqs)
  original_bc_keys <- format_barcode_string(
    as.character(compressed_original_bc_seqs),
    original_bc_gem_groups)

  candidate_reads <- sum(gbm@molecule_info$reads)
  candidate_read_frac <- candidate_reads / tot_reads

  if (target_type == 'raw_reads') {
    tgt_candidate_reads <- as.integer(round(target_reads_per_cell * n_cells * candidate_read_frac))
    candidate_rpc <- tgt_candidate_reads / n_cells
    raw_reads_per_cell <- target_reads_per_cell

  } else if (target_type == 'mapped_reads') {
    tgt_candidate_reads <- as.integer(round(target_reads_per_cell * n_cells))
    candidate_rpc <- tgt_candidate_reads / n_cells
    raw_reads_per_cell <- candidate_rpc / candidate_read_frac
  }

  if (tgt_candidate_reads > candidate_reads) {
    warning(sprintf("Target candidate reads (%d) is greater than available candidate reads (%s). Returning original matrix.",
                    tgt_candidate_reads, as.character(candidate_reads)))
    return(gbm)
  }
  subsample_rate <- tgt_candidate_reads / candidate_reads

  if (verbose) cat(sprintf("Subsampling: %0.2f\n", subsample_rate))
  bc_gene_umi_subsampled <- copy(gbm@molecule_info)
  bc_gene_umi_subsampled[,reads := rbinom(length(reads), as.integer(reads), subsample_rate)]

  if (verbose) cat("Sorting\n")
  setkey(bc_gene_umi_subsampled, barcode, gene)

  if (verbose) cat("Counting genes\n")
  # Create a column that can be compared to the keys for the original barcodes
  bc_gene_umi_subsampled[,bc_key := format_barcode_string(as.character(barcode), as.integer(gem_group))]

  # Select only barcodes in the matrix (might be filtered) and sum molecules per (bc,gene) tuple
  bc_gene_counts <- bc_gene_umi_subsampled[bc_key %in% original_bc_keys,
                                           j=list(count=sum(reads > 0)), by=c('bc_key', 'gene')]

  if (verbose) cat("Building subsampled matrix\n")
  new_matrix <- with(bc_gene_counts, sparseMatrix(i = 1 + gene,
                                                  j = match(bc_key, original_bc_keys),
                                                  x = count,
                                                  dims=dim(exprs(gbm))))
  gbm = newGeneBCMatrix(new_matrix, fData(gbm), pData(gbm), template=gbm)
  gbm@subsampled <- TRUE
  return(gbm)
}

#' Equalize a list of gbms by subsampling
#'
#' @param gbm_list A list of GeneBCMatrix objects to merge
#' @export
#' @return A list of GeneBCMatrix objects after subsampling
#' @examples
#' \dontrun{
#' equalize_gbms(list(gbm1,gbm2))
#' }
equalize_gbms <- function(gbm_list) {
  warning("This method of normalizing matrices is deprecated.\nPlease use the `cellranger aggr` pipeline (new in cellranger 1.2.0), which can combine arbitrary gene-barcode matrices\nand produce a combined, depth-normalized matrix.")
  n_gbms <- length(gbm_list)
  mean_reads <- sapply(gbm_list,get_mean_mapped_reads_per_cell)
  cat('Mean reads per cell before subsampling:\n')
  lapply(1:n_gbms, function(x) cat('gbm',x,':',mean_reads[x],'\n'))
  min_read <-  min(mean_reads)
  gbm_idx <- which(mean_reads == min_read)[1] # only use the first to normalize
  cat('Normalizing gbms to mean: ',min_read,' (',gbm_idx,')\n')
  gbm_others <- gbm_list[-gbm_idx]
  gbm_others <- lapply(gbm_others, function(x) subsample_gene_bc_matrix(x,
                                   target_reads_per_cell=min_read, target_type="mapped_reads"))
  gbm_out_list <- append(gbm_others,gbm_list[gbm_idx],gbm_idx-1)
  lapply(1:n_gbms, function(x) cat('gbm',x,':',median(colSums(gbm_list[[x]])),'\n'))
  cat('Mean UMIs per cell after subsampling:\n')
  lapply(1:n_gbms, function(x) cat('gbm',x,':',median(colSums(gbm_out_list[[x]])),'\n'))
  return(gbm_out_list)
}


#' Get a gene index from a valid gene symbol
#'
#' @param gbm A GeneBCMatrix to get the non-zero genes from
#' @param gene_symbol A gene symbol string
#' @return The indices from a gene index
get_gene_index <- function(gbm, gene_symbol) {
  index <- match(gene_symbol,fData(gbm)$symbol)
  if (is.na(index)) {
    warning(paste('did not find and will ignore gene symbol:', gene_symbol))
  } else if (length(index) > 1) {
    cat('Warning: found',length(index),'instances of',gene_symbol)
    cat('. Selecting the first instance by default.\n')
    index <- index[1]
  }
  return(index)
}

#' Truncate a GeneBCMatrix according to an array of gene symbols
#'
#' @param gbm A GeneBCMatrix
#' @param gene_symbols A list of gene symbols
#' @return A GeneBCMatrix truncated to genes with specific gene symbols
#' @export
#' @examples
#' gbm_trunc <- trunc_gbm_by_genes(gbm1,c('CD79A','NKG7'))
#'
trunc_gbm_by_genes <- function(gbm,gene_symbols) {
  gene_indices <- sapply(gene_symbols, function(x) get_gene_index(gbm,x))
  gene_indices <- gene_indices[!is.na(gene_indices)]
  gbm_trunc <- gbm[gene_indices,]
  return(gbm_trunc)
}

#' Get the non-zero gene indices from a GeneBCMatrix
#'
#' @param gbm a GeneBCMatrix to get the non-zero genes from
#' @return The gene (row) indices for genes with at least 1 count in the matrix
#' @export
#' @examples
#' non_zero_gene_ids <- get_nonzero_genes(gbm1)
#'
get_nonzero_genes <- function(gbm) {
  # Get the unique column indices in triplet representation
  return(sort(unique(1 + slot(as(exprs(gbm), 'dgTMatrix'), 'i'))))
}

# TODO this can have duplicate column/row names... (HOW FIX?)
# TODO this also does not properly handle most slots in the object...
#' Concatenate gene-barcode matrices
#'
#' @param gbms a list of GeneBCMatrix objects to concatenate
#' @return A new GeneBCMatrix containing all the barcodes in the given matrices
#' @export
#' @examples
#' big_gbm <- concatenate_gene_bc_matrices(list(gbm1, gbm2))
#'
concatenate_gene_bc_matrices <- function(gbms) {
  if (class(gbms) != 'list') {
    stop("'gbms' must be a list")
  }
  if (length(gbms) < 1) {
    stop("Fewer than one matrix was supplied for concatentation")
  }
  if (!all(unlist(lapply(gbms, class)) == 'GeneBCMatrix')) {
    stop("'gbms' must contain GeneBCMatrix objects")
  }
  if (!all(unlist(lapply(gbms, nrow)) == nrow(gbms[[1]]))) {
    stop("All matrices to concatenate must contain the same number of genes")
  }

  new_matrix = do.call(cBind, lapply(gbms, function(gbm) exprs(gbm)))
  new_pdata = do.call(rBind, lapply(gbms, function(gbm) pData(gbm)))
  colnames(new_matrix) = row.names(new_pdata)  # when samples overlap, rBind modifies names so reset to match

  result <- newGeneBCMatrix(mat=new_matrix,
                            fd=fData(gbms[[1]]),
                            pd=new_pdata)

  # TODO none of the other slots are handled correctly in original code...
  # Not sure what exactly the right way to combine molecule info, pipestance
  # paths, etc. is or if it really matters that this information is lost
  # on concat

  result@subsampled <- FALSE
  return(result)
}


#' Format a barcode string
#'
#' @param barcode_seq Cell barcode sequence
#' @param gem_group Gem group number
#' @return A formatted barcode string
format_barcode_string <- function(barcode_seq, gem_group) {
  sprintf("%s-%d", barcode_seq, gem_group)
}


#' Read value from metrics CSV
#'
#' @param x String to convert
#' @return Numeric version of string
csv_value_to_numeric <- function(x) {
  x <- gsub(",", "", x)
  if (substr(x, nchar(x), nchar(x)) == "%") {
    x <- sub("%", "", x)
    return(as.numeric(x) / 100)
  } else {
    return(as.numeric(x))
  }
}

#' Get a metric from the metric summary
#'
#' @param gbm A GeneBCMatrix object
#' @param metric_name Metric name
#' @return numeric metric value
get_metric <- function(gbm, metric_name) {
  if (!(metric_name %in% names(gbm@summary))) stop(sprintf("Could not find %s in summary", metric_name))
  return(csv_value_to_numeric(gbm@summary[[metric_name]][1]))
}

#' Get mean raw reads per cell
#'
#' @param gbm A GeneBCMatrix object
#' @return Mean raw reads per cell
#' @export
#' @examples
#' \dontrun{
#' raw_rpc <- get_mean_raw_reads_per_cell(gbm1)
#' }
get_mean_raw_reads_per_cell <- function(gbm) {
  if (length(gbm@summary) < 1) {
    stop("Object does not contain a metrics summary. This is likey because metrics_summary.csv was not specified when creating a gbm object.")
  }

  n_cells <- get_metric(gbm, CELL_COUNT_METRIC)
  tot_reads <- get_metric(gbm, TOTAL_READS_METRIC)

  return(as.numeric(tot_reads)/as.numeric(n_cells))
}

#' Get mean mapped reads per cell
#'
#' @param gbm A GeneBCMatrix object
#' @return Mean confidently mapped, barcoded reads per cell
#' @export
#' @examples
#' \dontrun{
#' mapped_rpc <- get_mean_mapped_reads_per_cell(gbm)
#' }
get_mean_mapped_reads_per_cell <- function(gbm) {
  raw_rpc <- get_mean_raw_reads_per_cell(gbm)

  conf_mapped_frac <- get_metric(gbm, CONF_MAPPED_READS_METRIC)
  valid_bc_frac <- get_metric(gbm, VALID_BARCODES_METRIC)

  return(raw_rpc * conf_mapped_frac * valid_bc_frac)
}

#' Extract the cell barcode sequence from one or more barcode strings
#'
#' @param barcode_str Formatted cell barcode sequence
#' @return A cell barcode sequence string
get_barcode_sequence <- function(barcode_str) {
  do.call(rbind, strsplit(barcode_str, '-'))[,1]
}

#' Extract the gem group from one or more barcode strings
#'
#' @param barcode_str Formatted cell barcode sequence
#' @return A gem group integer
get_gem_group <- function(barcode_str) {
  as.integer(do.call(rbind, strsplit(barcode_str, '-'))[,2])
}

#' Compress a vector of sequences into 2-bit representation
#'
#' Compress a vector of sequences into integer64 objects containing
#' a 2-bit representation. Ns are not allowed
#' @param seqs Vector of nucleotide sequences to compress
#' @return A vector of integer64 objects
compress_sequences <- function(seqs) {
  if (any(grepl('[^ACGT]', seqs))) {
    stop("At least one sequence contains Ns")
  }
  nuc_to_int <- as.integer(0:3)
  names(nuc_to_int) <- c('A', 'C', 'G', 'T')

  chars <- do.call(rbind, strsplit(seqs, ''))
  nuc_ints <- matrix(nuc_to_int[chars], nrow=length(seqs))
  result <- integer64(length(seqs))
  for(i in 1:ncol(nuc_ints)) {
    result <- result * as.integer(4) + nuc_ints[,i]
  }

  result
}

#' Write cluster-specific genes to file
#'
#' @param gene_clusters A list of genes that are prioritized (output from function: prioritize_genes)
#' @param gene_folder_path Path to folder to store the cluster-specific genes
#' @param n_genes Number of genes to include for each cluster
#' @param output_type The format of the output gene list "ensembl" or "symbol
#' @return matrix that includes all the genes written to a folder
#' @export
#' @examples
#' \dontrun{
#'  write_cluster_specific_genes(genes_to_plot,GENE_RES_FOLDER, n_genes=10,output_type='ensembl')
#'  }
write_cluster_specific_genes <- function(gene_clusters, gene_folder_path, n_genes=10,output_type='symbol') {
  cluster_names <- names(gene_clusters)
  output_files <- lapply(cluster_names,function(x) {
    file.path(gene_folder_path,paste('Cluster_',x,'_',output_type,'.tsv',sep=""))
    })

  if (output_type=='ensembl') {
    write_dummy <- lapply(1:length(cluster_names), function(x) {
      write.table(gene_clusters[[x]]$id[1:n_genes], file=output_files[[x]], quote=FALSE,col.names=FALSE,row.names=FALSE)
    })
  } else {
    write_dummy <- lapply(1:length(cluster_names), function(x) {
      write.table(gene_clusters[[x]]$symbol[1:n_genes], file=output_files[[x]], quote=FALSE,col.names=FALSE,row.names=FALSE)
    })
  }
  all_genes <- matrix(unlist(lapply(gene_clusters, function(x) x$symbol[1:n_genes])),ncol=length(cluster_names))
  colnames(all_genes) <-  lapply(cluster_names,function(x) paste('Cluster_',x,sep=""))
  write.table(all_genes,file=file.path(gene_folder_path,'Cluster_ALL.csv'),sep=",",quote=FALSE,row.names=FALSE)
  cat('Written files to folder:',gene_folder_path,'\n')
  return(all_genes)
}

#' Get the read count matrix
#'
#' @param gbm A GeneBCMatrix object
#' @return A dgTMatrix with the same dimension as exprs(gbm)
#' @export
#' @examples
#' \dontrun{
#'  get_read_count_mtx(gbm1)
#'  }
get_read_count_mtx <- function(gbm) {
  if (length(gbm@molecule_info) < 1) {
    stop("Molecule info not loaded. Please use load_molecule_info() before constructing read count matrix")
  }
  original_bc_seqs <- get_barcode_sequence(colnames(exprs(gbm)))
  original_bc_gem_groups <- get_gem_group(colnames(exprs(gbm)))
  compressed_original_bc_seqs <- compress_sequences(original_bc_seqs)
  original_bc_keys <- format_barcode_string(as.character(compressed_original_bc_seqs),original_bc_gem_groups)

  bc_molecule_info <- copy(gbm@molecule_info)
  cat("Sorting...\n")
  setkey(bc_molecule_info, barcode, gene)
  bc_molecule_info[,bc_key := format_barcode_string(as.character(barcode), as.integer(gem_group))]

  cat("Counting genes...\n")
  bc_gene_counts <- bc_molecule_info[bc_key %in% original_bc_keys, j=list(count=sum(as.integer(reads))), by=c('bc_key', 'gene')]

  cat("Building matrix...\n")
  new_matrix <- with(bc_gene_counts, sparseMatrix(i = 1 + gene,
                                                  j = match(bc_key, original_bc_keys), #  <- ensures the new columns match the original columns
                                                  x = count,
                                                  dims=dim(exprs(gbm))))
  return(new_matrix)
}


