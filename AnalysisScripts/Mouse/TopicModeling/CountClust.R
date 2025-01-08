## Count Clust Topics
library(Seurat)
library(CountClust)
library(ggplot2)
library(cowplot)
library(MatrixExtra)


topicModeling <- function(mat, n = 16, tol = 0.1, save.path) {
  
  # set parameters for topic modeling 
  n.topics <- as.numeric(n)
  tolerance <- as.numeric(tol)
  
  # create directories to save the results in 
  sub.dir <- paste0(n.topics, "topics_tol", tolerance)
  dir.create(file.path(save.path, sub.dir), recursive = T)
  
  FitGoM(mat, K = n.topics, tol = tolerance,
         path_rda = file.path(save.path, sub.dir, paste0('FitGoM_50K_k', n.topics, '_tol', tolerance, '.rda')))
  gc()
  gc()
}

subsample = TRUE # Sub-sample for testing
data <- readRDS("/gstore/project/neurodegen_meta/data/DS_50k.rds")
mat <- GetAssayData(data,slot = "counts",assay = "RNA")

mat <- mat[-which(rowSums(mat)==0)]
mat <- as.matrix(t_deep(mat))

out.dir <- "/gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Rscripts/TopicModeling"
ks = as.numeric(commandArgs(trailingOnly = TRUE))
dir.create(out.dir, recursive = T)

for (k in ks){
  fit <- topicModeling(mat, k, tol = 0.1,out.dir)
}


