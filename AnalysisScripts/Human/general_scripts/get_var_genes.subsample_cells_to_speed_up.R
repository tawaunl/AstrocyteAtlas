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

colnames(out[,out$studybatch == "kamath"]) <- paste("cell_", 1:ncol(out[,out$studybatch == "kamath"]), sep="") ### empty colnames in kamath object for some reason
colnames(out) <- paste(out$studybatch, colnames(out), sep="_")
print("Getting variable genes")



studies_long <- names(table(out$studybatch))[as.logical(table(out$studybatch) > 20000)]
studies_short <- names(table(out$studybatch))[as.logical(table(out$studybatch) < 20000)]


for (i in 1:length(studies_short)) {
  print(paste("Running var genes for", studies_short[i], sep=" "))
  var_genes_tmp <- getTopHVGs(out[,out$studybatch==studies_short[i]], n=2000)
  saveRDS(var_genes_tmp, paste("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/", studies_short[i], ".var_genes.rds", sep=""))
print(paste("Done saving  var genes for", studies_short[i], sep=" "))
}



for (i in 1:length(studies_long)) {
  print(paste("Running var genes for", studies_long[i], sep=" "))
  out_tmp <- out[,out$studybatch==studies_long[i]]
  var_genes_tmp <- getTopHVGs(out_tmp[,sample(colnames(out_tmp), 20000)], n=2000)
  saveRDS(var_genes_tmp, paste("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/", studies_long[i], ".var_genes.rds", sep=""))
 print(paste("Done saving  var genes for", studies_long[i], sep=" "))
}
