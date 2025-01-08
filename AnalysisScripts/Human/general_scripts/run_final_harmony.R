library(scran)
library(scater)
library(SingleCellExperiment)
library(batchelor)
library(scran.chan)
library(dplyr)
library(harmony)
library(BiocParallel)

se_ast_for_pca_post_clean <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.CLEANED_space.onlyHVG_for_PCA.rds")

se_ast_for_pca_post_clean <- RunHarmony(se_ast_for_pca_post_clean, group.by.vars = "donor", theta=1.5, reduction.save = "harmony_theta_1.5")
se_ast_for_pca_post_clean <- runUMAP(se_ast_for_pca_post_clean, dimred="harmony_theta_1.5", name = "UMAP_theta_1.5")

print("Done with harmony, saving the file!")
saveRDS(se_ast_for_pca_post_clean, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.CLEANED_space.onlyHVG_for_PCA.with_harmony.rds")


