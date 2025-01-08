library(batchelor)
library(BiocParallel)

out_saved <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes.cellbender_counts.20cell_filter.no_logcounts.liddelow_cellranger_rest_cellbender.rds")

out <- multiBatchNorm(out_saved, normalize.all = TRUE, batch = out_saved$studybatch, BPPARAM=MulticoreParam(workers = 14))

saveRDS(out, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes.cellbender_counts.20cell_filter.WITH_logcounts.liddelow_cellranger_rest_cellbender.rds")
