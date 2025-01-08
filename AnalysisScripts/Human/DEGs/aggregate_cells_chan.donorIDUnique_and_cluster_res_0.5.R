library(ShadowArray)
library(SingleCellExperiment)
library(scran.chan)
library(gp.sa.solo)
library(tidyr)
library(scater)

print("Loading in the se file")

se_clean_full <- readRDS("/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/another_run_after_removing_high_mito_clusters/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.WITH_logcounts.liddelow_cellranger_rest_cellbender.rds")

clusters <- readRDS("/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/another_run_after_removing_high_mito_clusters/integration/Clusters_res_harmony_2_donor_groupvar_0.5.rds")

se_clean_full$clusters_0.5 <- clusters$membership

print("starting initialize matrix call")
x <- initializeSparseMatrix(assay(se_clean_full, 1), num.threads = 14)
print("Done initializing the matrix!")

print("Aggregating cells...")
agg <- aggregateAcrossCells.chan(x, list(donor=se_clean_full$donorIdUnique, cluster = se_clean_full$clusters_0.5), num.threads=14)

saveRDS(agg, "aggregated_cells.AD_MS_PD.070724.donorIdUnique_and_cluster.rds")

print("Done!")
