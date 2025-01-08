library(SingleCellExperiment)
library(batchelor)
library(scran.chan)
library(dplyr)
library(scran)
library(scater)
library(harmony)
library(BiocParallel)

print("Reading in combined dataset file")

out <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes.cellbender_counts.20cell_filter.WITH_logcounts.liddelow_cellranger_rest_cellbender.rds")

clusters <- readRDS("/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/clustering/Clusters_res_0.5.rds")
donorIdUnique <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes.cellbender_counts.20cell_filter.COLDATA.donorIDUnique.rds")

print("removing high mito clusters and donors with less than 20 cells")
out$clusters <- clusters$membership
out$donorIdUnique <- donorIdUnique
out <- out[,! out$clusters %in% c(2, 5, 8)]
donors_to_keep <- names(table(out$donorIdUnique))[as.logical(table(out$donorIdUnique) >= 20)]
out <- out[,out$donorIdUnique %in% donors_to_keep]

print(dim(out))

print("Getting variable genes")

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

saveRDS(var_genes_list, "var_genes_list_with_chan.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.rds")

library(plyr)
### Combine all together

all_var_genes_combined <- plyr::join_all(var_genes_list, by='var_genes', type='full')
all_var_genes_combined <- all_var_genes_combined %>% mutate(median_var=rowMeans(all_var_genes_combined[,-1], na.rm = TRUE))

counts_ast_post_clean <- data.frame(gene = all_var_genes_combined$var_genes, Freq = rowSums(! is.na(all_var_genes_combined)) - 2) ### 1 for var_genes and 1 for median_var_column 

all_var_genes_combined_final <- left_join(all_var_genes_combined,counts_ast_post_clean,by=c("var_genes"="gene")) %>% dplyr::filter(var_genes %in% rownames(out))
var_genes_2k_ast_post_clean <- all_var_genes_combined_final %>% arrange(-Freq,median_var) %>% head(2000) %>% dplyr::select(var_genes)

saveRDS(var_genes_2k_ast_post_clean, "var_genes_2k_ast.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.rds")

