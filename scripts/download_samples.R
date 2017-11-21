library(cellrangerRkit)

# Download data used in the vignette

#Sys.setenv(RKIT_VIGNETTE_PIPESTANCE_BASE_PATH='/path/to/cellranger/pipestances')

PIPESTANCE_BASE_PATH <- Sys.getenv('RKIT_VIGNETTE_PIPESTANCE_BASE_PATH')

if (PIPESTANCE_BASE_PATH == '') {
  stop('Please set the RKIT_VIGNETTE_PIPESTANCE_BASE_PATH environment variable. (to, e.g., /path/to/cellranger/pipestances)')
}
if (!file.exists(PIPESTANCE_BASE_PATH)) {
  dir.create(PIPESTANCE_BASE_PATH, recursive=T)
}

download_sample(sample_name="pbmc3k", sample_dir=file.path(PIPESTANCE_BASE_PATH, 'pbmc3k'), host="http://s3-us-west-2.amazonaws.com/10x.files/samples/cell/", lite=F)
download_sample(sample_name="pbmc6k", sample_dir=file.path(PIPESTANCE_BASE_PATH, 'pbmc6k'), host="http://s3-us-west-2.amazonaws.com/10x.files/samples/cell/", lite=F)
