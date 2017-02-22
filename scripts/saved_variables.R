library(cellrangerRkit)

cat(sprintf("Loaded cellrangerRkit version %s\n", packageVersion('cellrangerRkit')))
cat("Building cellrangerRkit vignette...\n")

# Set RKIT_VIGNETTE_PIPESTANCE_BASE_PATH
# Set RKIT_VIGNETTE_CACHE_PATH
# 1) Run download_samples.R
# 2) Run saved_variables.R (the vignette also calls this)

## Preload data for the vignette.
# Take as input a cellranger pipestance and load variables used by the rkit vignette.
# Caches these data in an RDS file at a specified location.

#Sys.setenv(RKIT_VIGNETTE_PIPESTANCE_BASE_PATH='/path/to/rkit_vignette_data/pipestances')
#Sys.setenv(RKIT_VIGNETTE_CACHE_PATH='/path/to/rkit_vignette_data')

PIPESTANCE_BASE_PATH <- Sys.getenv('RKIT_VIGNETTE_PIPESTANCE_BASE_PATH')
CACHE_PATH <- Sys.getenv('RKIT_VIGNETTE_CACHE_PATH')

if (PIPESTANCE_BASE_PATH == '') {
  stop('Please set the RKIT_VIGNETTE_PIPESTANCE_BASE_PATH environment variable. (to, e.g., /path/to/rkit_vignette_data/pipestances)')
}
if (CACHE_PATH == '') {
  stop('Please set the VIGNETTE_CACHE_PATH environment variable. (to, e.g., /path/to/rkit_vignette_data). The output of this script goes there.')
}

if (!file.exists(CACHE_PATH)) {
  dir.create(CACHE_PATH, recursive=T)
}

get_pipestance_data <- function(name, pipestance_base_path=PIPESTANCE_BASE_PATH, cache_path=CACHE_PATH, force=FALSE) {
  # See if we need to generate a new rds
  rds_filename <- file.path(cache_path, sprintf('%s.rds', name))
  if (!force && file.exists(rds_filename)) {
    cat(sprintf('Using cached data for %s at %s\n', name, rds_filename))
    return(readRDS(rds_filename))
  }

  # Regenerate the RDS
  cat(sprintf('Did not find cached data for %s at %s\n', name, rds_filename))
  cat('Regenerating cached data...\n')
  pipestance_path <- file.path(pipestance_base_path, name)
  gbm <- load_cellranger_matrix(pipestance_path)
  print(gbm@pipestance_path)
  analysis_results <- load_cellranger_analysis_results(pipestance_path)
  saveRDS(list(gbm=gbm, analysis_results=analysis_results), file=rds_filename)
  cat('done.\n')
  return(readRDS(rds_filename))
}

# Test CR1.2.0 data
#pbmc4k_data <- get_pipestance_data('pbmc4k')

# Load CR1.1.0 vignette data (pbmc3k,pbmc6k)
pbmc3k_data <- get_pipestance_data('pbmc3k')
pbmc6k_data <- get_pipestance_data('pbmc6k')


# The vignette requires:
# gbm
# analysis_results
# merged_tsne_clust

gbm <- pbmc3k_data$gbm
analysis_results <- pbmc3k_data$analysis_results

# Generate merged_tsne_clust

get_merged_analysis <- function(gbm1, gbm2, cache_filename=merged_rds_filename) {
  if (!file.exists(cache_filename)) {
    cat(sprintf('Did not find cached data for pbmc3k_pbmc6k_merged.rds at %s\n', merged_rds_filename))
    cat('Regenerating cached data...\n')
    set.seed(0)
    gbm_list <- list(gbm1, gbm2)
    gbm_list <- lapply(gbm_list,load_molecule_info) # load sample molecule information
    gbm_list_equalized <- equalize_gbms(gbm_list)   # equalize the gene-barcode matrices
    merged_gbm <- concatenate_gene_bc_matrices(gbm_list_equalized)
    merged_ID <- unlist(lapply(1:length(gbm_list), function(x) rep(x,dim(gbm_list[[x]])[2])))

    set.seed(0)
    n_clust <- 5
    pca_result <- run_pca(merged_gbm)
    tsne_result <- run_tsne(pca_result)
    kmeans_result <- run_kmeans_clustering(pca_result, k=n_clust)
    # included re-named cell IDs, t-SNE and k-means result in a merged data frame
    merged_tsne_clust <- data.frame(Barcode=as.factor(1:tsne_result$N),
                                    TSNE.1=tsne_result$Y[,1],TSNE.2=tsne_result$Y[,2],
                                    Cluster=kmeans_result$cluster,Batch=merged_ID)
    saveRDS(merged_tsne_clust, file=cache_filename)
  }
  return(readRDS(merged_rds_filename))
}
gbm1 <- gbm
gbm2 <- pbmc6k_data$gbm
merged_rds_filename <- file.path(CACHE_PATH, 'pbmc3k_pbmc6k_merged.rds')
merged_tsne_clust <- get_merged_analysis(gbm1,gbm2, merged_rds_filename)
