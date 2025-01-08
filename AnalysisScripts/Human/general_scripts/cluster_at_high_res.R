library(SingleCellExperiment)
library(scater)
library(scran)
library(scran.chan)

print("Reading combined file")
out <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.with_logcounts.v3.liddelow_cellranger_rest_cellbender.rds")


print("Reading harmony output file")

se_ast_for_pca_post_clean <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.onlychanHVG_for_PCA.with_harmony.liddelow_cellranger_rest_cellbender.rds")

reducedDims(out) <- reducedDims(se_ast_for_pca_post_clean)

print("CLustering at high res...")
clusters <- clusterSNNGraph.chan(t(reducedDim(se_ast_for_pca_post_clean, "harmony_theta_1.5")), num.threads = 14, resolution = 3)

saveRDS(clusters, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/Clusters_res3.not_cleaned_space.rds")




