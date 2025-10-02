# Load necessary libraries
library(fastTopics)
library(ggplot2)
library(dplyr)
library(purrr)
library(CountClust)
library(Seurat) 

library(MatrixExtra)



# --- 1. Define the k values and model types to assess ---
# These should match the values used in your sbatch array job.
k_values <- c(2,3, 10,14,16,18,20,22,26)
data.dir <- "/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling"
out.dir <- "/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling/MouseResults"
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}


# --- 2. Load the data and create a results data frame ---
# Create an empty list to store all the results.
all_results <- list()
i <- 1
tol=0.1
data <- readRDS("/gstore/data/astroMetaAnalysis/data/DS_50k.rds")
models <- list()

for (k in k_values) {
  rda_fname <- file.path(data.dir,"results",paste0("GoM_k",k,"_results.rda"))
  load(rda_fname)
  if (file.exists(rda_fname)) {
    message(paste("Loading file:", rda_fname))
    load(rda_fname)
    omega <- Topic_clus$omega
    annotation <- data.frame(
      sample_id = paste0("X", c(1:NROW(omega))),
      tissue_label = data$finalClusters)
    rownames(omega) <- annotation$sample_id
    message(paste("Saving StructurePlot for k =", k))
    png(file=paste0(out.dir,"/StructurePlot_k",k,".png"))
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
}
i=1
for (k in k_values){
  message(paste("Processing k =", k))
  
  rda_fname <- file.path(data.dir,"results",paste0("GoM_k",k,"_results.rda"))
  load(rda_fname)
  
  load(rda_fname)
  models[[i]] <- Topic_clus
  
  usage <- as.data.frame(Topic_clus$omega)
  colnames(usage) <- paste0("lda_", colnames(usage))
  
  theta <- as.matrix(Topic_clus$theta)
  colnames(theta) <- paste0("lda_", colnames(theta))
  
  write.csv(usage, file.path(out.dir,paste0("usage_k",k ,".csv")))
  write.csv(theta, file.path(out.dir,paste0("theta_k",k ,".csv")))
  
  top_features_min <- ExtractTopFeatures(theta, top_features = 200, shared = T, method = "poisson", options = "min")
  write.csv(top_features_min, file.path(out.dir,paste0("score_min_k",k ,".csv"))
  )
  message(paste("Saved usage, theta, and top features for k =", k))
  
  i <- i + 1
}

counts <- GetAssayData(data,slot = "counts",assay = "RNA")
counts <- counts[rowSums(counts > 0) >= 10, ]
feat <- as.data.frame(genomitory::getFeatures("GMTY17:GRCm38/GRCm38.IGIS4.0.genes.rds@REVISION-3")) # get gene names
library(dplyr)
protein_coding <- feat %>% filter(type =="protein_coding")

counts <- counts[which(rownames(counts) %in% protein_coding$symbol),]


counts <- as.matrix(t_deep(counts))

out <- compGoM(counts, models)

saveRDS(out,file.path(out.dir, "compGoM_results.rds"))
message("compGoM results saved to:", file.path(out.dir, "compGoM_results.rds"))
n.topics <- ks
names(out) <- paste0("topic_", n.topics)
bic.plot <- sapply(names(out), function(x) out[[x]]$BIC)
saveRDS(bic.plot,file.path(out.dir, "bic_plot.rds"))
message("BIC values saved to:", file.path(out.dir, "bic_plot.rds"))

bic <- data.frame(BIC=bic.plot, k=n.topics)
write.csv(bic, file.path(out.dir, "bicValues.csv"))
message("BIC values saved to:", file.path(out.dir, "bicValues.csv"))