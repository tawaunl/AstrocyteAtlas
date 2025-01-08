library(SingleCellExperiment)
library(batchelor)
library(scran.chan)
library(dplyr)
library(scran)
library(scater)
library(harmony)
library(BiocParallel)

print("Reading in combined dataset file")

out <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/another_run_after_removing_high_mito_clusters/subclustering/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.Cluster_2_only.rds")

print("Getting var genes...")

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

studies <- unique(out$studybatch)

var_genes_list <- list()

for (i in 1:length(studies)) {
  print(paste("Running var genes for", studies[i], sep=" "))
  var_genes_list[[i]] <- get_var_genes(out[,out$studybatch==studies[i]],2000)
}


for (i in 1:length(studies)) {
    colnames(var_genes_list[[i]]) <- c("var_genes", paste("var.", studies[i], sep=""))
}

saveRDS(var_genes_list, "var_genes_list_with_chan.Cluster2.rds")

library(plyr)
### Combine all together

all_var_genes_combined <- plyr::join_all(var_genes_list, by='var_genes', type='full')
all_var_genes_combined <- all_var_genes_combined %>% mutate(median_var=rowMeans(all_var_genes_combined[,-1], na.rm = TRUE))

counts_ast_post_clean <- data.frame(gene = all_var_genes_combined$var_genes, Freq = rowSums(! is.na(all_var_genes_combined)) - 2) ### 1 for var_genes and 1 for median_var_column 

all_var_genes_combined_final <- left_join(all_var_genes_combined,counts_ast_post_clean,by=c("var_genes"="gene")) %>% dplyr::filter(var_genes %in% rownames(out))
var_genes_2k_ast_post_clean <- all_var_genes_combined_final %>% arrange(-Freq,median_var) %>% head(2000) %>% dplyr::select(var_genes)

saveRDS(var_genes_2k_ast_post_clean, "var_genes_2k_ast.Cluster2.rds")

