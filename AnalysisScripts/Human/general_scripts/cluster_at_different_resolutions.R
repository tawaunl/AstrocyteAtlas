library(scran.chan)
library(dplyr)
library(scran)
library(SingleCellExperiment)

message("Loading the file")
se_ast_for_pca_post_clean <- readRDS("se_ast_for_pca.post_clean.cells_to_keep.rds")
message("Done loading")
resolutions <- seq(from = 0.1, to = 1, by = 0.1)

for (res in resolutions) {
message(paste("Starting to run resolution", res, sep=" "))

clusters <- clusterSNNGraph.chan(t(reducedDim(se_ast_for_pca_post_clean, "harmony_theta1.5")), num.threads = 7, resolution = res)
message(paste("Saving data for resolution", res, sep=" "))
saveRDS(clusters, paste("Clusters_res_", res, ".rds",sep=""))
} 



