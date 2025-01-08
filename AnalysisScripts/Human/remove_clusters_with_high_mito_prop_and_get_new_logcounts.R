library(ShadowArray)
library(SingleCellExperiment)
library(scran.chan)
library(gp.sa.solo)
library(tidyr)
library(dplyr)

print("Loading in the se file")

se_clean_full <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes.cellbender_counts.20cell_filter.WITH_logcounts.liddelow_cellranger_rest_cellbender.rds")
clusters <- readRDS("/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/clustering/Clusters_res_0.5.rds")
donorIdUnique <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes.cellbender_counts.20cell_filter.COLDATA.donorIDUnique.rds")

print("removing high mito clusters and donors with less than 20 cells")
se_clean_full$clusters <- clusters$membership
se_clean_full$donorIdUnique <- donorIdUnique
se_clean_full <- se_clean_full[,! se_clean_full$clusters %in% c(2, 5, 8)]
donors_to_keep <- names(table(se_clean_full$donorIdUnique))[as.logical(table(se_clean_full$donorIdUnique) >= 20)]
se_clean_full <- se_clean_full[,se_clean_full$donorIdUnique %in% donors_to_keep]

### get logcounts

library(batchelor)
library(BiocParallel)
print("Getting logcounts")

out <- multiBatchNorm(se_clean_full, normalize.all = TRUE, batch = se_clean_full$studybatch, BPPARAM=MulticoreParam(workers = 14))

print("Done computing logcounts")
saveRDS(out, "AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.WITH_logcounts.liddelow_cellranger_rest_cellbender.rds")

