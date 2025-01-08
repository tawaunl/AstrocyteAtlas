library(scran)
library(scater)
library(SingleCellExperiment)
library(batchelor)
library(scran.chan)
library(dplyr)
library(harmony)
library(BiocParallel)

#out_clean <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes.cellbender_counts.20cell_filter.no_logcounts.liddelow_cellranger_rest_cellbender.rds")

#print("Getting logcounts")
#out <- batchelor::multiBatchNorm(out_clean, normalize.all = TRUE, batch = out_clean$studybatch, BPPARAM=MulticoreParam(workers = 14))
#print("Done")

out <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes.cellbender_counts.20cell_filter.WITH_logcounts.liddelow_cellranger_rest_cellbender.rds")

get_var_genes <- function(se,n_genes) {

  x <- initializeSparseMatrix(assay(se, 1), num.threads = 14)
  is.mito <- integer(0)
  qc.metrics <- perCellQCMetrics.chan(x,subsets = list(mito = is.mito),
                                    num.threads = 4)
  qc.filters <- perCellQCFilters.chan(
    sums = qc.metrics$sums,
    detected = qc.metrics$detected,
    subsets = qc.metrics$subsets,
    nmads = 3)

  qc.discard <- qc.filters$filters$overall
  x2 <- filterCells.chan(x, qc.discard)
  lib.sizes <- qc.metrics$sums[!qc.discard]
  normed <- logNormCounts.chan(x2, lib.sizes)

  variances <- modelGeneVar.chan(normed,num.threads = 14)
  variances <- variances$statistics
  keep <- rank(-variances$residuals, ties.method = "first") <= n_genes

  var_genes <- rownames(se)[keep]
  variances <- -variances$residuals[keep]
  return(data.frame(var_genes=var_genes, variances=variances))
}

studies_long <- names(table(out$studybatch))[as.logical(table(out$studybatch) > 20000)]
studies_short <- names(table(out$studybatch))[as.logical(table(out$studybatch) < 20000)]

var_genes_list_short <- list()
var_genes_list_long <- list()

for (i in 1:length(studies_short)) {
  print(paste("Running var genes for", studies_short[i], sep=" "))
  var_genes_list_short[[i]] <- get_var_genes(out[,out$studybatch==studies_short[i]],2000)
}


for (i in 1:length(studies_long)) {
  print(paste("Running var genes for", studies_long[i], sep=" "))
  out_tmp <- out[,out$studybatch==studies_long[i]]
  var_genes_list_long[[i]] <- get_var_genes(out_tmp[,sample(colnames(out_tmp), 20000)], n=2000)
}

studies <- c(studies_short, studies_long)
var_genes_list <- c(var_genes_list_short, var_genes_list_long)

for (i in 1:length(studies)) {
    colnames(var_genes_list[[i]]) <- c("var_genes", paste("var.", studies[i], sep=""))
}

saveRDS(var_genes_list, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes.cellbender_counts.20cell_filter.no_logcounts.liddelow_cellranger_rest_cellbender.rds")

library(plyr)
### Combine all together

all_var_genes_combined <- plyr::join_all(var_genes_list, by='var_genes', type='full')
all_var_genes_combined <- all_var_genes_combined %>% mutate(median_var=rowMeans(all_var_genes_combined[,-1], na.rm = TRUE))

counts_ast_post_clean <- data.frame(gene = all_var_genes_combined$var_genes, Freq = rowSums(! is.na(all_var_genes_combined)) - 2) ### 1 for var_genes and 1 for median_var_column

all_var_genes_combined_final <- left_join(all_var_genes_combined,counts_ast_post_clean,by=c("var_genes"="gene")) %>% dplyr::filter(var_genes %in% rownames(out))

var_genes_2k_ast_post_clean <- all_var_genes_combined_final %>% arrange(-Freq,median_var) %>% head(2000) %>% dplyr::select(var_genes)


saveRDS(var_genes_2k_ast_post_clean, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/var_genes_2k_ast_post_clean.CLEANED_space.wuth_subsampling_for_large_studies.rds")




print("Getting a logcounts matrix for PCA")

A_logcounts <-  as(assays(out)[[2]],  "sparseMatrix")
se_ast_for_pca_post_clean <- SingleCellExperiment(assays = list(logcounts = A_logcounts))
colData(se_ast_for_pca_post_clean) <- colData(out)
se_ast_for_pca_post_clean <- se_ast_for_pca_post_clean[var_genes_2k_ast_post_clean$var_genes,]
print("logcounts matrix preped, running the PCA")

se_ast_for_pca_post_clean <- fixedPCA(se_ast_for_pca_post_clean,BPPARAM=MulticoreParam(workers = 14))
print("done computing PCA, saving file!")

saveRDS(se_ast_for_pca_post_clean, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.CLEANED_space.onlyHVG_for_PCA.rds")

se_ast_for_pca_post_clean <- RunHarmony(se_ast_for_pca_post_clean, group.by.vars = "donor", theta=1.5, reduction.save = paste("harmony_theta_", res, sep=""))
se_ast_for_pca_post_clean <- runUMAP(se_ast_for_pca_post_clean, dimred=paste("harmony_theta_", 1.5, sep=""),name = paste("UMAP_theta_", res, sep=""))

print("Done with harmony, saving the file!")
saveRDS(se_ast_for_pca_post_clean, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.CLEANED_space.onlyHVG_for_PCA.with_harmony.rds")



