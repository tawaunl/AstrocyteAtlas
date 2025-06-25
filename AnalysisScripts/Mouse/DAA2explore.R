# Reviewewr Response

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(wesanderson)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 
library(scales)
library(ggrepel)
library(scran.chan)
library(scuttle)
library(edgeR)
source("~/Documents/scHelpers.R")
library(openxlsx)

# Load Data ====================================================================
dir <- "~/Documents/AstrocytePaper"
data <- readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))

all <- table(data$Sample)
all <- data.frame(all)
colnames(all) <- c("Sample","Total")

DAA2 <- subset(data, subset = ClusterNames=="DAA2")
library(scCustomize)
DAA2 <- Split_Layers(DAA2,split.by = "StudyID")


DAA2 <- runSeurat(DAA2)

DAA2 <- IntegrateLayers(object = DAA2, method = CCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose=FALSE,k.weight=40)

# re-join layers after integration
DAA2[["RNA"]] <- JoinLayers(DAA2[["RNA"]])

DAA2 <- FindNeighbors(DAA2, reduction = "integrated.cca", dims = 1:30)
DAA2 <- FindClusters(DAA2, resolution = c(0.1,0.15,0.20,.25,0.3,0.35,0.4,0.45,0.5))
DAA2 <- RunUMAP(DAA2,reduction = "integrated.cca",dims=1:30)
DimPlot(DAA2,group.by = "refinedClusters",label=T,label.size = 10.5,pt.size = 1)
DimPlot(DAA2,group.by = "StudyID")
DAA2 <- Add_Cell_QC_Metrics(object = DAA2, species = "mouse")
data <- Add_Cell_QC_Metrics(data,species = "mouse")
DAA2 <- subset(DAA2,subset= percent_mito <10)
data <- subset(data,subset= percent_mito <10)
DAA2$refinedClusters <- DAA2$RNA_snn_res.0.5
DAA2$refinedClusters <- gsub("\\<4\\>", "0", DAA2$refinedClusters)
DAA2$refinedClusters <- gsub("\\<6\\>", "2", DAA2$refinedClusters)
DAA2$refinedClusters <- gsub("\\<5\\>", "2", DAA2$refinedClusters)
DAA2$refinedClusters <- factor(DAA2$refinedClusters)




DotPlot(DAA2,features = c("Myoc","Crym","Meg3","Ttr"),
        group.by ="refinedClusters",cols = "RdYlBu",assay = "RNA" ) + RotatedAxis() +coord_flip() + xlab('Gene') +  ylab('Cluster') +
  theme(axis.text = element_text(size = 14),axis.title = element_text(size=18,face='bold'))  


DotPlot(DAA2,features = "Crym",group.by ="refinedClusters" )
DotPlot(DAA2,features = "Meg3",group.by ="RNA_snn_res.0.5" )
DotPlot(DAA2,features = "Ttr",group.by ="RNA_snn_res.0.5" )


coldata <- as.data.frame(DAA2@meta.data)
study_cluster_tmp_by_study <- coldata %>%
  group_by(StudyName, refinedClusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

df_normalized <- study_cluster_tmp_by_study %>%
  group_by(refinedClusters) %>%
  mutate(TotalCellsInCluster = sum(freq),
         NormalizedCells = freq / TotalCellsInCluster)
#### Plot study Composition -----------

ggplot(df_normalized, aes(fill=StudyName, y=NormalizedCells, x=refinedClusters)) +
  geom_bar(position="stack", stat="identity", alpha=0.6) + theme_void()+ theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14)) +
  xlab("Clusters") + ylab("Proportion of cells within a study") +
  theme(axis.title.y =  element_text(angle = 90))



QC_Plots_Combined_Vln(seurat_object = DAA2, pt.size = 0.1,group.by ="refinedClusters",
                      colors_use = scales::hue_pal()(7))
QC_Plots_Combined_Vln(seurat_object = data, pt.size = 0,group.by ="finalClusters",
                      colors_use = scales::hue_pal()(4))

counts <- round(GetAssayData(DAA2, slot="counts", assay="RNA"))   
counts <- initializeSparseMatrix(counts)
lognorm <- logNormCounts.chan(counts,batch = DAA2$StudyID)
markers <- scoreMarkers.chan(lognorm,
                             groups=DAA2$refinedClusters, batch=DAA2$StudyID,lfc=0)
saveRDS(markers,"~/Documents/AstrocytePaper/DAA2subclusterMarkers.rds")
# Get top 5 markers from each cluster 
topmarkers <- c()
for (cluster in 1:length(markers$statistics)) {
  clustermarkers <- data.frame(markers$statistics[[cluster]])
  clustermarkers <- clustermarkers[order(clustermarkers$logFC,decreasing = TRUE),]
  top <- rownames(clustermarkers)[1:5]
  topmarkers <- c(topmarkers,top)
}

gcSample <- list()
for (cluster in 1:length(markers$statistics)) {
  clustermarkers <- data.frame(markers$statistics[[cluster]])
  clustermarkers <- clustermarkers[order(clustermarkers$logFC,decreasing = TRUE),]
  geneList <- as.numeric(clustermarkers$logFC[which(clustermarkers$logFC > 0)])
  names(geneList) <- rownames(clustermarkers[which(clustermarkers$logFC > 0),])
  
  gcSample[[paste0("Cluster_",cluster)]] <- geneList
  
}

library(clusterProfiler)

ck <- lapply(names(gcSample),function(x){
  gseGO(gcSample[[x]],ont = "BP",
        OrgDb = org.Mm.eg.db::org.Mm.eg.db,
        keyType = "SYMBOL",eps=0, pAdjustMethod = "BH")
})
names(ck) <- names(gcSample)
gsea_results <- lapply(names(ck), function(cluster){
  ck[[cluster]]@result
})

names(gsea_results) <- names(ck)
DotPlotCompare(
  gsea_list = gsea_results,
  n = 10,
  size_col = "NES",
  color_col = "p.adjust",
  size_cutoff = NULL,      # Filter to pathways with NES >= 1.5
  color_cutoff = 0.01,
  direction = "positive" # Filter to pathways with p.adjust <= 0.05
)

topmarkers <- c(topmarkers,"Thbs4")

plot <- DotPlot(DAA2,features = c(unique(topmarkers)), assay = "RNA",
                cols = "RdYlBu",group.by = "refinedClusters",scale = TRUE)

plot + RotatedAxis() +coord_flip() + xlab('Gene') +  ylab('Cluster') +
  theme(axis.text = element_text(size = 14),axis.title = element_text(size=18,face='bold'))  

sce <- DAA2 %>% 
  Seurat::as.SingleCellExperiment() 

library(scran)
library(purrr)
library(dplyr)
library(AnnotationDbi)
features <- markers$statistics %>% 
  map(~ .x |>
        as.data.frame() |>
        arrange(desc(cohen.mean))|>
        dplyr::slice(1:15) |>
        rownames())|> 
  unlist2()

features[1:3] <- c("Gfap","Crym","Hsp90ab1")
rowData(sce)$Marker <- features[match(rownames(sce), features)] |>
  names() %>% 
  factor(levels = levels(sce$refinedClusters))

scDotPlot::scDotPlot(sce,features = unique(features),
                       group = "refinedClusters",
                       #block = "Sample",
                       scale = TRUE,
                       cluster = FALSE,
                       groupAnno = "refinedClusters",
                       featureAnno = "Marker",
                       featureLegends = FALSE,
                       annoHeight = 0.025,
                       dotColors=c("#313695FF","#FFFFBFFF","#A50026FF"),
                       annoColors = list("refinedClusters" = scales::hue_pal()(7),
                                         "Marker" = scales::hue_pal()(7)),
                       annoWidth = 0.1,fontSize = 20) 


VlnPlot(data,"Ttr",group.by = "finalClusters",slot = "counts")
Dim
## Renaming DAA2 subclusters ---------
DAA2$subclusters <- DAA2$refinedClusters
DAA2$subclusters <- gsub("\\<0\\>", "NT-Transport", DAA2$subclusters)
DAA2$subclusters <- gsub("\\<1\\>", "Fth1+", DAA2$subclusters)

DAA2$subclusters <- gsub("\\<2\\>", "OxPhos", DAA2$subclusters)


DAA2$subclusters <- gsub("\\<3\\>", "Gfap+", DAA2$subclusters)




BRumaps <- lapply(unique(DAA2$BrainRegion), function(x){
  
  plot_umap_library(sobj=DAA2,meta="BrainRegion",sample=x)
})
names(BRumaps) <- unique(DAA2$BrainRegion)

plot_grid(BRumaps$Hippocampus,BRumaps$Cortex,BRumaps$CNS,BRumaps$`Corpus Callosum`,
          ncol=2)


abundances<- table(DAA2$refinedClusters,factor(DAA2$Sample))
df_long <- GetAbundanceTable(abundances,DAA2@meta.data,samplecol = "Sample",
                             subset.cols = c("Sample","Disease","BrainRegion","StudyName"),
                             calc.col = "refinedClusters")
level_order <- c("AD_control","AD","MS_control","MS")
df_long$Disease <- factor(df_long$Disease,levels = level_order) 

fill.codes <- c("white","darkgreen","white", "purple")
color.codes <- c("darkgreen","darkgreen","purple","purple")
zone <- levels(factor(df_long$Disease))

filter <- c("Zhao 5XFAD","Choi 5XFAD","Shi P301S")

df_long <- df_long[! df_long$StudyName %in% filter,]

ggplot(df_long, aes(x=Disease, y=Percent,fill=Disease,colour = Disease)) + theme_classic() +
  geom_jitter(color="darkgrey",width = 0.4,size=3,alpha=0.6)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  facet_wrap(~Cluster,scales = "free",ncol=4) + xlab("Disease Label") +
  theme(axis.text.y = element_text(size=20)) +
  ylab("Cellularity Proportion") +  
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(values=setNames(fill.codes, zone))+
  scale_color_manual(values=setNames(color.codes, zone))

# LPS scoring ----------
lps.de <- read.csv("~/Documents/AstrocytePaper/Figure3/LPSvsSalineDE.csv",row.names = 1)
lps.genes <- rownames(lps.de[lps.de$diffexpressed=="LPS",])
score_intersect <- intersect(lps.genes,rownames(DAA2))
scores_by_cell <- colMeans(
  GetAssayData(DAA2,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
DAA2[["LPSCurrentAstroScore"]] <- scores_by_cell


scvi_coords <- get_scvi_coords(DAA2,DAA2$RNA_snn_res.0.5)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))
text <- DAA2$RNA_snn_res.0.5
text_x <- vapply(split(scvi_coords$UMAP1, text), median, FUN.VALUE=0)
text_y <- vapply(split(scvi_coords$UMAP2, text), median, FUN.VALUE=0)

ggplot(scvi_coords%>%
         arrange(LPSCurrentAstroScore),aes(x=UMAP1, y=UMAP2, colour=LPSCurrentAstroScore)) +
  geom_point() +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) +
  scale_colour_gradientn(limits= c(0,1), oob=scales::squish,
                         colours =  RColorBrewer::brewer.pal(5,"Purples"))


#AD Scoring -----------------------
AD.DE <- readRDS(file.path(dir, "res_dsl_AD_full.011824.rds"))
MS.DE <- readRDS(file.path(dir, "res_dsl_MS_full.011824.rds"))
AD.DE_plot <- generatePlotTable(AD.DE)
AD.DE_plot$diffexpressed <- "unchanged"

AD.DE_plot$diffexpressed[AD.DE_plot$dl_mu >= 0.25 &
                           AD.DE_plot$FDR<=0.05 &
                           AD.DE_plot$sig=="yes" ] <- "up"

AD.DE_plot$diffexpressed[AD.DE_plot$dl_mu <= -0.25 &
                           AD.DE_plot$FDR<=0.05 &
                           AD.DE_plot$sig=="yes" ] <- "down"

AD.genes <-AD.DE_plot$symbol[AD.DE_plot$diffexpressed=="up"]
score_intersect <- intersect(AD.genes,rownames(DAA2))
scores_by_cell <- colMeans(
  GetAssayData(DAA2,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
DAA2[["ADCurrentAstroScore"]] <- scores_by_cell


library(scales)
scvi_coords <- get_scvi_coords(DAA2,DAA2$RNA_snn_res.0.5)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))
text <- DAA2$RNA_snn_res.0.5
text_x <- vapply(split(scvi_coords$UMAP1, text), median, FUN.VALUE=0)
text_y <- vapply(split(scvi_coords$UMAP2, text), median, FUN.VALUE=0)

ggplot(scvi_coords%>%
         arrange(ADCurrentAstroScore),aes(x=UMAP1, y=UMAP2, colour=ADCurrentAstroScore)) +
  geom_point(size=3) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) +
  scale_colour_gradientn(limits= c(0,.1), oob=squish,
                         colours =  RColorBrewer::brewer.pal(9,"Purples"))

# MS Scoring --------------
MS.DE_plot <- generatePlotTable(MS.DE)
MS.DE_plot$diffexpressed <- "unchanged"

MS.DE_plot$diffexpressed[MS.DE_plot$dl_mu >= 0.25 &
                           MS.DE_plot$FDR<=0.05 &
                           MS.DE_plot$sig=="yes" ] <- "up"

MS.DE_plot$diffexpressed[MS.DE_plot$dl_mu <= -0.25 &
                           MS.DE_plot$FDR<=0.05 &
                           MS.DE_plot$sig=="yes" ] <- "down"

MS.genes <-MS.DE_plot$symbol[MS.DE_plot$diffexpressed=="up"]
score_intersect <- intersect(MS.genes,rownames(DAA2))
scores_by_cell <- colMeans(
  GetAssayData(DAA2,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
DAA2[["MSCurrentAstroScore"]] <- scores_by_cell


scvi_coords <- get_scvi_coords(DAA2,DAA2$subclusters)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))
text <- DAA2$subclusters
text_x <- vapply(split(scvi_coords$UMAP1, text), median, FUN.VALUE=0)
text_y <- vapply(split(scvi_coords$UMAP2, text), median, FUN.VALUE=0)

ggplot(scvi_coords%>%
         arrange(MSCurrentAstroScore),aes(x=UMAP1, y=UMAP2, colour=MSCurrentAstroScore)) +
  geom_point(size=3) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) +
  scale_colour_gradientn(limits= c(0,1), oob=squish,
                         colours =  RColorBrewer::brewer.pal(9,"Purples"))

VlnPlot(DAA2,features= c("LPSCurrentAstroScore","ADCurrentAstroScore",
                         "MSCurrentAstroScore"),
        group.by = "refinedClusters",add.noise = T,raster = T)

ggplot(DAA2@meta.data,aes(x=refinedClusters,y=ADCurrentAstroScore,fill = refinedClusters))+
  geom_boxplot()

# Abundance plot -----------------------

abundances<- table(DAA2$subclusters,factor(DAA2$Sample))

df_long <- GetAbundanceTable(abundances,coldata = DAA2@meta.data,samplecol = "Sample",
                             subset.cols = c("Sample","Disease","BrainRegion","StudyName"),
                             calc.col = "RNA_snn_res.0.5")
df_long$Disease <- factor(df_long$Disease,levels = level_order) 

ggplot(df_long, aes(x=Disease, y=Percent,fill=Disease,colour = Disease))+ theme_classic() +
  geom_jitter(color="darkgrey",width = 0.4,size=3,alpha=0.6)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  facet_wrap(~Cluster,scales = "free",ncol=4) + xlab("Disease Label") +
  theme(axis.text.y = element_text(size=20)) +
  ylab("Cellularity Proportion") +  
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(values=setNames(fill.codes, zone))+
  scale_color_manual(values=setNames(color.codes, zone))

df_transformed <- GetAbundanceTable(abundances,DAA2@meta.data,samplecol = "Sample",
                                    subset.cols = c("Sample","Disease","BrainRegion","StudyName"),
                                    calc.col = "subclusters",return_clr = T)

sig.clusters <- list()
for(cluster in 1:length(levels(factor(DAA2$subclusters)))){
  clusterName <- levels(factor(DAA2$subclusters))[cluster]
  dataCluster <- df_transformed[which(df_transformed$Cluster==clusterName),]
  res <- kruskal.test(Percent ~ Disease, data = dataCluster)
  # check if there is any significance
  if(res$p.value < 0.8/length(levels(factor(DAA2$subclusters)))){
    sig.clusters[[clusterName]] <- res
    # do a Welch's t-test to find which specific groups are different
    res.AD <- wilcox.test(dataCluster$Percent[which(dataCluster$Disease=="AD_control")],
                          dataCluster$Percent[which(dataCluster$Disease=="AD")])
    res.MS <- wilcox.test(dataCluster$Percent[which(dataCluster$Disease=="MS_control")],
                          dataCluster$Percent[which(dataCluster$Disease=="MS")])
    if(res.MS$p.value < 0.9){
      sig.clusters[[clusterName]][["MS"]] <- res.MS
    }
    if(res.AD$p.value < 0.9){
      sig.clusters[[clusterName]][["AD"]] <- res.AD
    }
  }
}


FeaturePlot(DAA2,"MyocScore")
saveRDS(DAA2,"~/Documents/AstrocytePaper/DAA2subcluster.rds")
