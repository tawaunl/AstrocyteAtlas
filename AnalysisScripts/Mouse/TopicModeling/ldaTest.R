library(Seurat)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(foreach)
library(MatrixExtra)

topicModeling <- function(mat, n = 16, tol = 0.1, save.path) {
  
  # set parameters for topic modeling 
  n.topics <- as.numeric(n)
  tolerance <- as.numeric(tol)
  
  # create directories to save the results in 
  sub.dir <- paste0(n.topics, "topics_tol", tolerance)
  dir.create(file.path(save.path, sub.dir), recursive = T)
  
  fit <- fit_topic_model(mat,
                         k = n.topics, 
                         method.main = "em",
                         method.refine = "scd",
                         numiter.main = 100,
                         numiter.refine = 100,)
  saveRDS(fit,file = file.path(save.path, sub.dir, paste0('FitTopic_k', n.topics, '_tol', tolerance, '.RDS')))
  return(fit)
}


subsample = T # Sub-sample for testing
data <- readRDS("/gstore/project/neurodegen_meta/data/AstrocyteIntegration_AmbientRemoved_filtered_noneuron.RDS")
downsampled.obj <- data[, sample(colnames(data), size = 100000, replace=F)]
rownames(downsampled.obj@assays[["RNA"]]@meta.features) <- rownames(data)

counts <- GetAssayData(downsampled.obj,slot = "counts",assay = "RNA")
counts <- counts[-which(rowSums(counts)==0)]
counts <- t_deep(counts)
out.dir <- "/gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Rscripts/TopicModeling"
ks = c(4)

dir.create(out.dir, recursive = T)

for (k in ks){
  fit <- topicModeling(counts, k, tol = 0.1,out.dir)
  
  plot_progress(fit,x = "iter",add.point.every = 10,colors = "black") +
    theme_cowplot(font_size = 10)+ ggtitle(paste("Convergence Model fitting K =", k))
}
saveRDS(downsampled.obj,"/gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Rscripts/TopicModeling/object.RDS")