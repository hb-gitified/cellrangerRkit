library(cellrangerRkit)
context('tsne')

pipestance_path <- system.file("extdata", "test_pipestance", package="cellrangerRkit")

gbm <- load_cellranger_matrix(pipestance_path)
pca <- run_pca(gbm, 2)
tsne_res <- run_tsne(pca, theta=0.75)
