library(SingleCellExperiment)
library(scran.chan)
library(scran)
se_ast_for_pca_post_clean <- fixedPCA(se_ast_for_pca_post_clean,BPPARAM=MulticoreParam(workers = 14))
print("done computing PCA, saving file!")

saveRDS(se_ast_for_pca_post_clean, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.onlychanHVG_for_PCA.liddelow_cellranger_rest_cellbender.rds")


### Let's try multiple integrations

resolutions <- c(0.5, 1, 1.5, 2, 2.5, 3)

print("Running through multiple resolutions")

for (res in resolutions) {

se_ast_for_pca_post_clean <- RunHarmony(se_ast_for_pca_post_clean, group.by.vars = "donor",theta=res, reduction.save = paste("harmony_theta_", res, sep=""))
se_ast_for_pca_post_clean <- runUMAP(se_ast_for_pca_post_clean, dimred=paste("harmony_theta_", res, sep=""),name = paste("UMAP_theta_", res, sep=""))

}

print("Done, saving the file!")

saveRDS(se_ast_for_pca_post_clean, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.onlychanHVG_for_PCA.with_harmony.liddelow_cellranger_rest_cellbender.rds")



