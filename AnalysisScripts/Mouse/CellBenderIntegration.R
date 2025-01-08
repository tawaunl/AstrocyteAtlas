### ------ Load Libraries --------------------------------------------------------
pkgs <- c("Seurat","SingleCellExperiment","ggplot2","R.utils",
          "sctransform","future","harmony","Matrix","patchwork",
          "glmGamPoi","cowplot")
for(p in pkgs){
  suppressPackageStartupMessages(library(p, quietly=TRUE,character.only=TRUE))
}
options(future.globals.maxSize= 1e+11)
### ------ Load data ----------------------------------------------------
rm(list=ls())
gc()
# Load all astrocyte data
dsids <- c("GSE140511","GSE143758","GSE150934",
           "GSE161224","GSE166261","NGS2330",
           "NGS2453","NGS2664","NGS2722",
           "NGS3033","NGS3111","NGS3971",
           "GSE150934","GSE129763","GSE129609")

data.list <-lapply(X=dsids, FUN=function(x){
  DietSeurat(readRDS(file =  paste0("/gstore/project/neurodegen_meta/data/cellbender/",x,"/Astrocytes.RDS")))
})

# name everything with the same assay name "RNA"
for (study in 1:length(data.list)) {
  data <- data.list[[study]]
  tryCatch({
   data <-  RenameAssays(data,originalexp="RNA")
  }, error=function(e){})
  data.list[[study]] <- data
}
feat <- as.data.frame(genomitory::getFeatures("GMTY17:GRCm38/GRCm38.IGIS4.0.genes.rds@REVISION-3"))

# Name everything with symbols first to calculate MitoPercent
for (study in 1:length(data.list)) {
  data <- data.list[[study]]
  counts <- GetAssayData(data, assay = "RNA")
  if(startsWith(rownames(data)[1],"ENSMUSG")){
    x <- match(rownames(counts),rownames(feat))
  }
  else{x <- match(rownames(counts),feat$symbol)}
  features <- feat[x,]
  rownames(counts) <- make.unique(features$symbol)
  if(length(which(is.na(rownames(counts))==TRUE))==0){
    new <- CreateSeuratObject(counts)
  }
  else{
    counts <- counts[-which(is.na(rownames(counts))==TRUE),]
    new <- CreateSeuratObject(counts) }
  
  new@assays[["RNA"]]@meta.features <- features
  data.list[[study]] <- new
}

### ------ Seurat Integration CCA----------------------------------------


plan("multicore", workers = 7)
features <- SelectIntegrationFeatures(object.list = data.list )
anchors <- FindIntegrationAnchors(object.list = data.list,anchor.features = features)
cca.integrated <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the
DefaultAssay(cca.integrated) <- "integrated"

plan("sequential")
# Run the standard workflow for visualization and clustering
cca.integrated <- ScaleData(cca.integrated, verbose = FALSE)
cca.integrated <- RunPCA(cca.integrated, npcs = 30, verbose = FALSE)
cca.integrated <- RunUMAP(cca.integrated, reduction = "pca", dims = 1:20)
cca.integrated <- FindNeighbors(cca.integrated, reduction = "pca", dims = 1:20)
cca.integrated <- FindClusters(cca.integrated,resolution = 1.5)
saveRDS(cca.integrated, file ="/gstore/project/neurodegen_meta/data/AstrocyteIntegration_AmbientRemoved.RDS" )


#Plotting data
# p1 <- DimPlot(cca.integrated,shuffle=TRUE, reduction = "umap",
#               group.by = "StudyName",pt.size = 1) + ggtitle("CCA Integrated")
# p2 <- DimPlot(cca.integrated, reduction = "umap", label = TRUE, repel = TRUE)
# split <- DimPlot(cca.integrated, reduction = "umap",
#                  label = FALSE, group.by = "StudyName", split.by = "StudyName",
#                  label.size = 10, ncol=4, pt.size = 1) 
# 
# split2 <- DimPlot(cca.integrated, reduction = "umap",
#                   label = FALSE, group.by = "Model",split.by = "Model",
#                   label.size = 10, ncol = 4, pt.size = .8) 
# 
# model <- DimPlot(cca.integrated, reduction = "umap",
#                  label = FALSE, group.by = "Model",
#                  label.size = 10, pt.size = .8) 
# 
# split3 <- DimPlot(cca.integrated, reduction = "umap",
#                   label = FALSE, group.by = "BrainRegion",split.by = "BrainRegion",
#                   label.size = 10,ncol = 3, pt.size = 1) 
# 
# region <- DimPlot(cca.integrated, reduction = "umap",
#                   label = FALSE, group.by = "BrainRegion",
#                   label.size = 10, pt.size = 1)
# 
# 
# pdf(file = "/gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Figures/CCA_integratedUMAP_AmbientRemoved.pdf", width=20,height = 10)
# print(p1 + p2 )
# dev.off()
# 
# pdf(file = "/gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Figures/CCA_integratedUMAP_byStudy_AmbientRemoved.pdf", width=35,height = 16)
# print(plot_grid(split, p1,ncol = 2))
# dev.off()
# 
# pdf(file = "/gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Figures/CCA_integratedUMAP_byModel_AmbientRemoved.pdf", width=35,height = 16)
# print(plot_grid(split2, odel,ncol = 2))
# dev.off()
# 
# pdf(file = "/gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Figures/CCA_integratedUMAP_byRegion_AmbientRemoved.pdf", width=35,height = 16)
# print(plot_grid(split3, region,ncol = 2))
# dev.off()
# 
