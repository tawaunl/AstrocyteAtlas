library(SingleCellExperiment)
library(scran.chan)
library(dplyr)
library(scran)
library(scater)

out <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.with_logcounts.v3.liddelow_cellranger_rest_cellbender.rds")
print("done reading in combined file")

var_genes_2k_ast_post_clean <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/var_genes_2k_ast_post_clean_with_chan.liddelow_cellranger_rest_cellbender.rds")

A_logcounts <-  as(assays(out)[[2]],  "sparseMatrix")
se_ast_for_pca_post_clean <- SingleCellExperiment(assays = list(logcounts = A_logcounts))
colData(se_ast_for_pca_post_clean) <- colData(out)
se_ast_for_pca_post_clean <- se_ast_for_pca_post_clean[var_genes_2k_ast_post_clean$var_genes,]
print("logcounts matrix preped")

saveRDS(se_ast_for_pca_post_clean,"/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/se_ast_for_pca_post_clean.onlyHVGs.liddelow_cellranger_rest_cellbender.rds")
