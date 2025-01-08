library(SingleCellExperiment)
library(scater)
print("Reading combined file")
out <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.with_logcounts.v3.liddelow_cellranger_rest_cellbender.rds")


print("Reading harmony output file")

se_ast_for_pca_post_clean <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.onlychanHVG_for_PCA.with_harmony.liddelow_cellranger_rest_cellbender.rds")

reducedDims(out) <- reducedDims(se_ast_for_pca_post_clean)

clusters <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/Clusters_res3.not_cleaned_space.rds")

scores <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/Cell_types_scores_for_space_cleaning.rds")

out$clusters_res3 <- clusters$membership


png("/gne/web/dev/apache/htdocs/people/novikovg/Astrocytes_meta/Harmony_AD_MS_PD_integration_2024/UMAP_clusters_resolution_3.png", width = 2400, height = 1500, res = 250)
plotReducedDim(out, dimred = "UMAP_theta1.5", colour_by = "clusters_res3", point_size = 0.1, point_alpha = 0.2)
dev.off()




