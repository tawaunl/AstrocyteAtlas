library(SingleCellExperiment)
library(batchelor)
library(scran.chan)
library(dplyr)
library(scran)
library(scater)
library(harmony)
library(BiocParallel)

print("Getting a logcounts matrix")

out <- readRDS("/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/another_run_after_removing_high_mito_clusters/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.WITH_logcounts.liddelow_cellranger_rest_cellbender.rds")

var_genes_2k_ast_post_clean <- readRDS("/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/another_run_after_removing_high_mito_clusters/integration/var_genes_2k_ast.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.rds")

A_logcounts <-  as(assays(out[var_genes_2k_ast_post_clean$var_genes,])[[2]],  "sparseMatrix")
se_ast_for_pca_post_clean <- SingleCellExperiment(assays = list(logcounts = A_logcounts))
colData(se_ast_for_pca_post_clean) <- colData(out)
print("logcounts matrix preped")

se_ast_for_pca_post_clean <- fixedPCA(se_ast_for_pca_post_clean,BPPARAM=MulticoreParam(workers = 14))
print("done computing PCA, saving file!")

saveRDS(se_ast_for_pca_post_clean, "AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.onlychanHVG_for_PCA.liddelow_cellranger_rest_cellbender.rds")

print("Running harmony with a default resolution")

se_ast_for_pca_post_clean <- RunHarmony(se_ast_for_pca_post_clean, group.by.vars = "donorIdUnique",theta=2, reduction.save = "harmony_theta_2_donor_groupvar")
se_ast_for_pca_post_clean <- runUMAP(se_ast_for_pca_post_clean, dimred="harmony_theta_2_donor_groupvar", name = "UMAP_theta_2_donor_groupvar")


print("Done, saving the file!")

saveRDS(se_ast_for_pca_post_clean, "AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.onlychanHVG_for_PCA.with_harmony.liddelow_cellranger_rest_cellbender.rds")

print("Done saving file, starting to run clustering")
resolutions <- seq(from = 0.3, to = 0.8, by = 0.1)

for (res in resolutions) {
message(paste("Starting to run resolution", res, sep=" "))

clusters <- clusterSNNGraph.chan(t(reducedDim(se_ast_for_pca_post_clean, "harmony_theta_2_donor_groupvar")), num.threads = 14, resolution = res)
message(paste("Saving data for resolution", res, sep=" "))
saveRDS(clusters, paste("Clusters_res_harmony_2_donor_groupvar_", res, ".rds",sep=""))
}

rm(se_ast_for_pca_post_clean)

print("starting initialize matrix call")
x <- initializeSparseMatrix(assay(out, 2), num.threads = 14)
print("Done initializing the matrix!")

resolutions <- seq(from = 0.3, to = 0.8, by = 0.1)

for (res in resolutions) {

print(paste("Getting markers for ", res, sep=""))

res_tmp <- readRDS(paste("Clusters_res_harmony_2_donor_groupvar_", res, ".rds", sep=""))

m.out <- scoreMarkers.chan(x, res_tmp$membership, batch = out$studybatch)
saveRDS(m.out, paste("m.out.harmony_2_donor_groupvar.", res, ".rds", sep=""))
}

print("Done!")




