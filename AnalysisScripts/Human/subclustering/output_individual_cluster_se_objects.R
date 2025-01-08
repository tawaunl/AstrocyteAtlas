library(ShadowArray)
library(SingleCellExperiment)
library(scran.chan)
library(gp.sa.solo)
library(tidyr)
library(scater)

print("Loading the se file")

se_ast <- readRDS("/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/another_run_after_removing_high_mito_clusters/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.WITH_logcounts.liddelow_cellranger_rest_cellbender.rds")
print("Done!")
clusters <- readRDS("/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/another_run_after_removing_high_mito_clusters/integration/Clusters_res_harmony_2_donor_groupvar_0.5.rds")

se_ast$clusters_res_0.5 <- clusters$membership

for (i in 1:8) {
print(i)
tmp <- se_ast[,se_ast$clusters_res_0.5 == i]
assays(tmp) <- assays(tmp)[1]
print(paste("Outputing file for cluster", i, sep=" "))

saveRDS(tmp, paste("/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/another_run_after_removing_high_mito_clusters/clustering/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.Cluster_", i, "_only.rds", sep=""))

}
