# Write Tables and Extract BIC 
library(Seurat)
library(fastTopics)
library(ggplot2)
library(CountClust)
library(MatrixExtra)

out.dir <- "/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling"
ks = c(4,6,8,10,12,14)
models <- list()
i <- 1
tol=0.1
data <- readRDS("/gstore/data/astroMetaAnalysis/data/AstrocyteIntegration_AmbientRemoved_filtered_noneuron.RDS")

for (k in ks){
  k_dir <- file.path(out.dir,"results",paste0("fasttopics_k",k,"_results"))
  if (!dir.exists(k_dir)) {
    dir.create(k_dir, recursive = TRUE)
  }
  
  rda_fname <- file.path(out.dir,"results",paste0("fasttopics_k",k,"_results.rds"))
  
  res <- readRDS(rda_fname)
  
  omega <- Topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = data$finalClusters)
  rownames(omega) <- annotation$sample_id
# Save Structureplot
  png(file=paste0(k_dir,"/StructurePlot.png"))
  print(StructureGGplot(omega = omega,
                  annotation = annotation,
                  palette = scales::hue_pal()(dim(omega)[2]),
                  yaxis_label = "Seurat Cluster",
                  order_sample = TRUE,
                  axis_tick = list(axis_ticks_length = .1,
                                   axis_ticks_lwd_y = .1,
                                   axis_ticks_lwd_x = .1,
                                   axis_label_size = 7,
                                   axis_label_face = "bold"),
                  figure_title = paste0("Topic Model K = ",k),
                  legend_key_size = 1, legend_text_size = 14) +
    guides(fill=guide_legend(title="Topic Number")))
  dev.off()
}

i=1
for (k in ks){
  k_dir <- paste0(out.dir, "/", k, "topics_tol", tol)
  rda_fname <- paste0(k_dir, "/", "FitGoM_50K_k", k, "_tol", tol , ".rda")
  
  load(rda_fname)
  models[[i]] <- Topic_clus
  
  usage <- as.data.frame(Topic_clus$omega)
  colnames(usage) <- paste0("lda_", colnames(usage))
  
  theta <- as.matrix(Topic_clus$theta)
  colnames(theta) <- paste0("lda_", colnames(theta))
  
  write.csv(usage, file.path(k_dir, "usage_50k.csv"))
  write.csv(theta, file.path(k_dir, "theta_50k.csv"))
  
  top_features_min <- ExtractTopFeatures(theta, top_features = 200, shared = T, method = "poisson", options = "min")
  write.csv(top_features_min, file.path(k_dir, "score_min_50k.csv"))
  
  i <- i + 1
}

counts <- GetAssayData(data,slot = "counts",assay = "RNA")
counts <- counts[-which(rowSums(counts)==0),]
counts <- as.matrix(t_deep(counts))

out <- compGoM(counts, models)

saveRDS(out,"/gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Rscripts/TopicModeling/comGOM_50k.rds")
n.topics <- ks
names(out) <- paste0("topic_", n.topics)
bic.plot <- sapply(names(out), function(x) out[[x]]$BIC)
saveRDS(bic.plot,"/gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Rscripts/TopicModeling/bicplot_50k.rds")

bic <- data.frame(BIC=bic.plot, k=n.topics)
write.csv(bic, file.path(out.dir, "bic_50k.csv"))