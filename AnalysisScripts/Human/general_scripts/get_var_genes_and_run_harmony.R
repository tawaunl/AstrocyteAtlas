library(SingleCellExperiment)
library(batchelor)
library(scran.chan)
library(dplyr)
library(scran)
library(scater)
library(harmony)
library(BiocParallel)
print("Reading in combined dataset file")
out <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.with_logcounts.rds")
print("Done!")

print("Getting variable genes")


var_genes_list <- list()

studies <- unique(out$studybatch)


for (i in 1:length(studies)) {
  print(paste("Running var genes for", studies[i], sep=" "))
  var_genes_list[[i]] <- getTopHVGs(out[,out$studybatch==studies[i]],n=2000)
}


for (i in 1:length(studies)) {
    colnames(var_genes_list[[i]]) <- c("var_genes", paste("var.", studies[i], sep=""))
}

saveRDS(var_genes_list, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/var_genes_list.rds")

library(plyr)
### Combine all together

all_var_genes_combined <- plyr::join_all(var_genes_list, by='var_genes', type='full')
all_var_genes_combined <- all_var_genes_combined %>% mutate(median_var=rowMeans(all_var_genes_combined[,-1], na.rm = TRUE))

counts_ast_post_clean <- data.frame(gene = all_var_genes_combined$var_genes, Freq = rowSums(! is.na(all_var_genes_combined)) - 2) ### 1 for var_genes and 1 for median_var_column 

all_var_genes_combined_final <- left_join(all_var_genes_combined,counts_ast_post_clean,by=c("var_genes"="gene")) %>% dplyr::filter(var_genes %in% rownames(out))

var_genes_2k_ast_post_clean <- all_var_genes_combined_final %>% arrange(-Freq,median_var) %>% head(2000) %>% dplyr::select(var_genes)


saveRDS(var_genes_2k_ast_post_clean, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/var_genes_2k_ast_post_clean.rds")


print("Getting a logcounts matrix")

A_logcounts <-  as(assays(out)[[2]],  "sparseMatrix")
se_ast_for_pca_post_clean <- SingleCellExperiment(assays = list(logcounts = A_logcounts))
colData(se_ast_for_pca_post_clean) <- colData(out)
se_ast_for_pca_post_clean <- se_ast_for_pca_post_clean[var_genes_2k_ast_post_clean$var_genes,]
print("logcounts matrix preped")

se_ast_for_pca_post_clean <- fixedPCA(se_ast_for_pca_post_clean,BPPARAM=MulticoreParam(workers = 14))
print("done computing PCA, saving file!")

saveRDS(se_ast_for_pca_post_clean, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.onlyHVG_for_PCA.rds")


### Let's try multiple integrations

resolutions <- c(0.5, 1, 1.5, 2, 2.5, 3)

print("Running through multiple resolutions")

for (res in resolutions) {

se_ast_for_pca_post_clean <- RunHarmony(se_ast_for_pca_post_clean, group.by.vars = "donor",theta=res, reduction.save = paste("harmony_theta_", res, sep=""))
se_ast_for_pca_post_clean <- runUMAP(se_ast_for_pca_post_clean, dimred=paste("harmony_theta_", res, sep=""),name = paste("UMAP_theta_", res, sep=""))

}

print("Done, saving the file!")

saveRDS(se_ast_for_pca_post_clean, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.onlyHVG_for_PCA.with_harmony.rds")

