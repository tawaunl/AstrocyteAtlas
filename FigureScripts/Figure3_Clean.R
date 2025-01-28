# Script to generate all plots in Figure 3

# ======================== Load Libraries ======================================
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

# Load Data ====================================================================
dir <- "~/Documents/AstrocytePaper"
data <- readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))


# A. DE DAA1 vs DAA2=========================================================
studies <- unique(data$StudyName)
studies <- studies[!studies %in% c("Shi P301S")] #Shi only has two samples

counts <- GetAssayData(object = data,
                       slot = "counts",
                       assay="RNA")
se <- SingleCellExperiment(assays = list(counts = counts),
                           colData = data@meta.data)
se <- se[,se$StudyName%in% studies]

summed <- aggregateAcrossCells(se, 
                               id=colData(se)[,c("StudyName", "Sample","finalClusters")])

summed.filt <- summed[,summed$ncells >= 10]
current <- summed.filt
y <- DGEList(counts(current), samples=colData(current))
keep <- filterByExpr(y, group=current$finalClusters)
y <- y[keep,]
y <- calcNormFactors(y)

design <- model.matrix(~0+factor(finalClusters), y$samples)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
colnames(design) <- c("Intercept","cluster2","cluster3","cluster4")
my.contrasts <- makeContrasts(DAA1vsDAA2=cluster3-cluster2, levels = design)

res <- glmQLFTest(fit, coef=ncol(design),contrast = my.contrasts)

res <- topTags(res,n=60000,p.value = 2)
res <- res$table
res <- res %>% mutate(sig = ifelse((PValue<=0.05 | PValue <= 0.05) & abs(logFC) >= 1, "yes", "no"))


res$diffexpressed <- "unchanged"

res$diffexpressed[res$logFC >= 1 &
                    res$PValue<=0.05 &
                    res$sig=="yes" ] <- "DAA1"

res$diffexpressed[res$logFC <= -1 &
                    res$PValue<=0.05 &
                    res$sig=="yes" ] <- "DAA2"

mycolors <- c("red", "goldenrod1", "grey")
names(mycolors) <- c("DAA2", "DAA1", "unchanged")

labels <- c("Gfap","Id3","Apoe","Clu","Aldoc","Cst3","Aqp4")

# Define thresholds for labeling outliers
positive_logFC_thresh <- 2  # Threshold for upregulated genes
negative_logFC_thresh <- -1.5 # Threshold for downregulated genes
pvalue_thresh <- 1e-28       # Common p-value threshold for both groups
#### Plot A ----------------
res <- read.csv(file.path(dir,"Figure3",paste0("DAA1vsDAA2_DEtable.csv")))
res$label <- with(res, 
                     ((logFC > positive_logFC_thresh) | 
                        (FDR < pvalue_thresh) |
                        (logFC < negative_logFC_thresh ) |
                        (symbol %in% labels)))

cairo_pdf(file.path(dir,"Figure3","DAA1vsDAA2_VolcanoPlot.pdf" ),
    width = 15,height=15)
ggplot(data=res, aes(x=logFC, y=-log10(FDR), col=diffexpressed, label=symbol)) + 
  ggrastr::geom_point_rast(alpha=0.6,size=3.5) + 
  theme_classic()+ geom_text_repel(size=11,data = subset(res, label & diffexpressed %in% c("DAA1","DAA2")),
                                   show.legend = FALSE, point.padding = 2,box.padding =1.2) +
  scale_colour_manual(values = mycolors) + ggtitle("DAA1 vs DAA2")+
  theme(legend.position = "none",
        legend.title = element_text(size=18,face='bold'),
        legend.key.size = unit(1 ,'cm'),
        axis.title = element_text(size=24),
        plot.title = element_text(size=28,hjust = 0.5)) +
  xlab(bquote(log[2]("Fold Change"))) + ylab(bquote(-log[10](PValue)))
dev.off()

#write.csv(res, file.path(dir,"Figure3",paste0("DAA1vsDAA2_DEtable.csv")))

# B. Pathways ================================================================
library(clusterProfiler)
res <- read.csv(file.path(dir,"Figure3",paste0("DAA1vsDAA2_DEtable.csv")))

DAA1.genes <- res$symbol[which(res1$diffexpressed=="DAA1")]

DAA2.genes <- res$symbol[which(res1$diffexpressed=="DAA2")]

ck <- compareCluster(list(DAA1=DAA1.genes,DAA2=DAA2.genes),
                     fun = "enrichGO",keyType="SYMBOL",
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db,ont="BP")

####  Plot B ----------- 
cairo_pdf(file.path(dir,"Figure3",paste0("DAA1vsDAA2_Pathways.pdf")),
    width = 4.5,height=7)
dotplot(ck)
dev.off()

# C. DAAs vs Homeostatic ------------------------------------
## Read DE ----------
DAA1 <- readRDS(file.path(dir,"DAA1vsHomeostatic_metaDE.rds"))
DAA2 <- readRDS(file.path(dir,"DAA2vsHomeostatic_metaDE.rds"))

DAA1 <- generatePlotTable(DAA1)
DAA2 <- generatePlotTable(DAA2)
DAA1$dl_mu <- DAA1$dl_mu *-1
DAA2$dl_mu <- DAA2$dl_mu *-1

DAA1$diffexpressed <- "unchanged"

DAA1$diffexpressed[DAA1$dl_mu   >= 0.5 &
                     DAA1$PValue <= 0.05 &
                     DAA1$sig    == "yes" ] <- "DAA1"

DAA1$diffexpressed[DAA1$dl_mu   <= -0.5 &
                     DAA1$PValue <= 0.05 &
                     DAA1$sig    == "yes" ] <- "Homeostatic"

DAA2$diffexpressed <- "unchanged"

DAA2$diffexpressed[DAA2$dl_mu   >= 0.25 &
                     DAA2$PValue <= 0.05 &
                     DAA2$sig    == "yes" ] <- "DAA2"

DAA2$diffexpressed[DAA2$dl_mu   <= -0.5 &
                     DAA2$PValue <= 0.05 &
                     DAA2$sig   == "yes" ] <- "Homeostatic"

## Pathways --------------------
DAA1.gsea <- setNames(DAA1$dl_mu,DAA1$ID)
DAA2.gsea <- setNames(DAA2$dl_mu,DAA2$ID)


DAA1.gsea <- DAA1.gsea[order(DAA1.gsea,decreasing = T)]
DAA2.gsea <- DAA2.gsea[order(DAA2.gsea,decreasing = T)]
ck <- compareCluster(list(DAA1=DAA1.gsea,DAA2=DAA2.gsea),
                     fun = "gseGO",OrgDb = org.Mm.eg.db::org.Mm.eg.db,ont="BP",eps=0)


DAA1.gsea <- DAA1.gsea[complete.cases(DAA1.gsea)]
DAA2.gsea <- DAA2.gsea[complete.cases(DAA2.gsea)]

D1 <- gseGO(DAA1.gsea,OrgDb = org.Mm.eg.db::org.Mm.eg.db,ont="BP",eps=0,keyType = "ENSEMBL")
D2 <- gseGO(DAA2.gsea,OrgDb = org.Mm.eg.db::org.Mm.eg.db,ont="BP",eps=0,keyType = "ENSEMBL")

go_DAA1 <- dotplot(D1, x="NES",showCategory=5) +
  ggtitle(str_wrap("GO:BP Pathways Enriched in DAA1 vs Homeo",width=30)) + 
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=24))
D2 <- filter(D2, NES>0)
go_DAA2 <- dotplot(D2, x="NES",showCategory=5) +
  ggtitle(str_wrap("GO:BP Pathways Enriched in DAA2 vs Homeo.",width=30)) + 
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size=10)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=24))


cairo_pdf("~/Documents/AstrocytePaper/Figure3/indvPAthwaysHOMEO_DAA1vsDAA2.pdf",
          width=11,height=5)
go_DAA1+go_DAA2
dev.off()

mycolors <- c("dodgerblue", "goldenrod1", "grey")
names(mycolors) <- c("Homeostatic", "DAA1", "unchanged")
### Plot C Left -------
positive_logFC_thresh <- 2  # Threshold for upregulated genes
negative_logFC_thresh <- -1 # Threshold for downregulated genes
pvalue_thresh <- 1e-28  
DAA1$label <- with(DAA1, 
                  ((dl_mu > positive_logFC_thresh) | 
                     (FDR < pvalue_thresh) |
                     (dl_mu < negative_logFC_thresh ) |
                     (symbol %in% labels)))

lfc_max <- max(abs(DAA1$dl_mu),na.rm = T) * 1.1  # Adding 10% for padding
x_limits <- c(-lfc_max+3, lfc_max)

cairo_pdf(file.path(dir,"Figure3","DAA1vsHomeostatic_VolcanoPlot.pdf" ),
          width = 15,height=15)
ggplot(data=DAA1, aes(x=dl_mu, y=-log10(FDR), col=diffexpressed, label=symbol)) + 
  ggrastr::geom_point_rast(alpha=0.6,size=3.5) + 
  theme_classic()+ geom_text_repel(size=11,
                                   data = subset(DAA1,label & diffexpressed %in% c("DAA1","Homeostatic")),
                                   show.legend = FALSE,
                                   box.padding = 1,
                                   max.overlaps = 15,
                                   point.padding = 1,
                                   min.segment.length = 0.1) +
  scale_colour_manual(values = mycolors) + ggtitle("DAA1 vs Homeostatic")+
  theme(legend.position = "none",
        legend.title = element_text(size=18,face='bold'),
        legend.key.size = unit(1 ,'cm'),
        axis.title = element_text(size=24),
        plot.title = element_text(size=28,hjust = 0.5)) +
  xlab("Avg. LogFoldChange") + ylab(bquote(bold(-log[10](PValue)))) +
  scale_x_continuous(limits = x_limits) 
dev.off()


###   Plot C Middle ---------------------------

# Define thresholds for labeling outliers
positive_logFC_thresh <- 1  # Threshold for upregulated genes
negative_logFC_thresh <- -3 # Threshold for downregulated genes
pvalue_thresh <- 1e-60         # Common p-value threshold for both groups

# Create a logical vector for labeling outlier points
DAA2$label <- with(DAA2, 
                   (dl_mu > positive_logFC_thresh | 
                      dl_mu < negative_logFC_thresh) | 
                     FDR < pvalue_thresh)

mycolors <- c("dodgerblue", "red", "grey")

names(mycolors) <- c("Homeostatic", "DAA2", "unchanged")
cairo_pdf(file.path(dir,"Figure3","DAA2vsHomeostatic_VolcanoPlot.pdf" ),
    width = 15,height=15)
ggplot(data=DAA2, aes(x=dl_mu, y=-log10(FDR), col=diffexpressed, label=symbol)) + 
  ggrastr::geom_point_rast(alpha=0.6,size=3.5) + 
  theme_classic()+ geom_text_repel(size=11,
                                   data = subset(DAA2,label),
                                   show.legend = FALSE,
                                   box.padding = 1,
                                   max.overlaps = 15,
                                   point.padding = 1,
                                   min.segment.length = 0.1) +
  scale_colour_manual(values = mycolors) + ggtitle("DAA2 vs Homeostatic")+
  theme(legend.position = "none",
        legend.title = element_text(size=18,face='bold'),
        legend.key.size = unit(1 ,'cm'),
        axis.title = element_text(size=24),
        plot.title = element_text(size=28,hjust = 0.5)) +
  xlab("Avg. LogFoldChange") + ylab(bquote(bold(-log[10](PValue))))
dev.off()

###  Plot C Right -----------

cairo_pdf(file.path(dir,"Figure3","DAA1_DAA2vsHomeo_Pathways.pdf"),
          width=6.5,height=10)
dotplot(ck_up, by="NES")  + 
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=20))+xlab("")

dev.off()

write.csv(ck@result,file.path(dir,"Figure3","DAAs_vs_HomeoPathways.csv"))

# D LPS Data UMAPS---------------------------------
## Recluster LPS Astros =====================================
lps.data <- readRDS(file.path(dir,"Figure3","LPS_astrocytes.rds"))
library(scCustomize)
library(scDblFinder)

sce <- as.SingleCellExperiment(lps.data)
dbl <- scDblFinder(sce,clusters = sce$seurat_clusters,samples = sce$orig.ident)


counts <- GetAssayData(lps.data,layer="counts")
meta <- lps.data@meta.data
astros <- CreateSeuratObject(counts = counts,meta.data = meta)
astros <- runSeurat(astros)
# remove bad clusters
toRemove <- c(9,10)
astros <- subset(astros,subset=seurat_clusters %in% toRemove,invert=T)
set.seed(824)
astros <- runSeurat(astros)
astros <- FindClusters(astros,resolution = 0.15)
#rename clusters
astros$seurat_clusters <- gsub("0", "A", astros$seurat_clusters)
astros$seurat_clusters <- gsub("1", "B", astros$seurat_clusters)
astros$seurat_clusters <- gsub("2", "C", astros$seurat_clusters)
astros$seurat_clusters <- gsub("3", "D", astros$seurat_clusters)
astros$seurat_clusters <- gsub("4", "E", astros$seurat_clusters)

astros$seurat_clusters <- gsub("A", "1", astros$seurat_clusters)
astros$seurat_clusters <- gsub("B", "2", astros$seurat_clusters)
astros$seurat_clusters <- gsub("C", "3", astros$seurat_clusters)
astros$seurat_clusters <- gsub("D", "4", astros$seurat_clusters)
astros$seurat_clusters <- gsub("E", "5", astros$seurat_clusters)
### Plot UMAPs =======================================
scvi_coords <- get_scvi_coords(astros,astros$seurat_clusters)

png(file.path(dir,"Figure3,UMAP_clusters.png"),
    width = 1500,height=1500,res = 300)
ggplot(data=scvi_coords, aes(x=UMAP1, y=UMAP2, colour  = seurat_clusters)) + 
  geom_point(size=2,alpha = 0.9) +theme_void()
dev.off()
p1 <- ggplot(data=scvi_coords, aes(x=UMAP1, y=UMAP2)) +
  geom_point(color="grey", size=2,alpha = 0.7) +
  geom_point(data = scvi_coords %>% dplyr::filter(TREATMENT_NAME=="wt_wt|LPS"), color=hue_pal()(2)[2], size=1, alpha=0.8)+
  theme_classic() + ggtitle("LPS") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 24, face = "bold"))
p2 <- ggplot(data=scvi_coords, aes(x=UMAP1, y=UMAP2)) +
  geom_point(color="grey", size=2,alpha = 0.7) +
  geom_point(data = scvi_coords %>% dplyr::filter(TREATMENT_NAME=="wt_wt|Saline"), color=hue_pal()(2)[1], size=1, alpha=0.8)+
  theme_classic() + ggtitle("Saline") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 24, face = "bold"))

png(file.path(dir,"Figure3","UMAP_Treatment.png"),
    width = 3200,height=1500,res = 300)
plot_grid(p1,p2)
dev.off()

saveRDS(astros, file.path(dir,"Figure3","ProcessedLPS_astros.rds"))
# E. LPS Abundance Plot ----------------------------
astros <- readRDS(file.path(dir,"Figure3","ProcessedLPS_astros.rds"))
sce <- as.SingleCellExperiment(astros)
dbl <- scDblFinder(sce,clusters = sce$seurat_clusters,samples = sce$orig.ident)
astros@meta.data <- as.data.frame(colData(dbl))
DimPlot(astros,group.by = "scDblFinder.class")

library(edgeR)
library(RColorBrewer)
library(viridis)
abundances <- table(astros$seurat_clusters,astros$Sample) 
abundances <- unclass(abundances) 

extra.info <- astros@meta.data[match(colnames(abundances), astros$Sample),]
d <- DGEList(abundances, samples=extra.info)
d = calcNormFactors(d)
d= estimateCommonDisp(d, verbose=TRUE)



norm_counts <- as.data.frame(t(d$counts)) 
colnames(norm_counts) <- paste("Cluster ", colnames(norm_counts), sep="")
norm_counts <- norm_counts %>% mutate(Sample = rownames(.))
coldata <- as.data.frame(astros@meta.data)
coldata_short <- coldata %>% dplyr::select(Sample,Dose,TREATMENT_NAME,Mouse) %>% unique


df_long_final <- left_join(norm_counts,coldata_short, by="Sample") 

percentages <- data.frame(matrix(nrow = dim(df_long_final)[1],ncol = dim(df_long_final)[2]))
sums <- data.frame(clustersums=rowSums(df_long_final[,-c(dim(df_long_final)[2]-3,
                                                         dim(df_long_final)[2]-2,
                                                         dim(df_long_final)[2]-1,
                                                         dim(df_long_final)[2])]),ident=df_long_final$Sample)

for (clust in 1:length(levels(factor(astros$seurat_clusters)))) {
  for (sample in 1:length(df_long_final[,clust])) {
    percent <- df_long_final[sample,clust]/sums[sample,1]
    percentages[sample,clust]<- percent
  }
}

percentages <- percentages +min(percentages[percentages>0],na.rm =T ) # value to add to offset zeros

#add col data back in 
percentages[,dim(df_long_final)[2]-3]<- df_long_final[,dim(df_long_final)[2]-3]
percentages[,dim(df_long_final)[2]-2]<- df_long_final[,dim(df_long_final)[2]-2]
percentages[,dim(df_long_final)[2]-1]<- df_long_final[,dim(df_long_final)[2]-1]
percentages[,dim(df_long_final)[2]]<- df_long_final[,dim(df_long_final)[2]]


colnames(percentages) <-  colnames(df_long_final)
## Clr transformation for statistics
clr <- percentages
for (sample in 1:dim(percentages)[1]) {
  tran <- log((percentages[sample,1:4]/exp(mean(as.numeric(log( percentages[sample,1:4]))))))
  #tran <- clr(percentages[sample,1:4])
  clr[sample,1:4]<- tran
}

df_long <- percentages %>% 
  pivot_longer(c(c(paste0("Cluster ",1:(dim(df_long_final)[2]-4)))))
df_transformed <- clr %>% 
  pivot_longer(c(c(paste0("Cluster ",1:(dim(df_long_final)[2]-4)))))
colnames(df_long)[5:6] <- c("Cluster","Percent")
colnames(df_transformed)[5:6] <- c("Cluster","Percent")

level_order <- c("Saline","LPS")
df_long$Dose <- factor(df_long$Dose,levels = level_order) 
df_long$Cluster <- factor(df_long$Cluster,levels =c(paste0("Cluster ",1:(dim(df_long_final)[2]-4))))
df_transformed$Dose <- factor(df_transformed$Dose,levels = level_order) 
df_transformed$Cluster <- factor(df_transformed$Cluster,levels =c(paste0("Cluster ",1:(dim(df_long_final)[2]-4))))

fill.codes <- hue_pal()(2)
color.codes <- hue_pal()(2)
zone <- levels(factor(df_long$Dose))

####Plot E ---------
cairo_pdf(file.path(dir,"Figure3","E_AbundancePlot.pdf"),
    width=14,height = 4)
ggplot(df_long, aes(x=Dose, y=Percent,fill=Dose,colour = Dose)) +
  theme_classic() +
  geom_jitter(color="darkgrey",width = 0.4,size=3,alpha=0.6)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  facet_wrap(~Cluster,scales = "free",ncol=5) + xlab("Dose") +
  theme(axis.text.y = element_text(size=20)) +
  ylab("Cellularity Proportion") +  
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values=setNames(fill.codes, zone))+
  scale_color_manual(values=setNames(color.codes, zone))
dev.off()

## Significant testing ==================
#Sig test on CLR transformed values
sig.clusters <- list()

for(cluster in 1:length(levels(factor(astros$seurat_clusters)))){
  clusterName <- paste0("Cluster ",levels(factor(astros$seurat_clusters))[cluster])
  dataCluster <- df_transformed[which(df_transformed$Cluster==clusterName),]
  res <- aov(Percent ~ Dose, data = dataCluster)
  # check if there is any significance
  if(summary(res)[[1]][["Pr(>F)"]][1] < 0.05/length(levels(factor(astros$seurat_clusters)))){
    sig.clusters[[clusterName]] <- res
    # do a Welch's t-test to find which specific groups are different
    LPS <- t.test(dataCluster$Percent[which(dataCluster$Dose=="Saline")],
                          dataCluster$Percent[which(dataCluster$Dose=="LPS")])
  
    if(LPS$p.value < 0.05){
      sig.clusters[[clusterName]][["LPS"]] <- LPS
    }
  }
}

write.csv(df_transformed1,"~/Documents/AstrocytePaper/2024-11-26 lucast3/Fig2D_AbundancePlotDataforStatistics_Mouse.csv")


# F. DE Scoring -----------------------
counts <- GetAssayData(object = astros,
                       slot = "counts",
                       assay="RNA")
se <- SingleCellExperiment(assays = list(counts = counts),
                           colData = astros@meta.data)

summed <- aggregateAcrossCells(se, 
                               id=colData(se)[,c("Sample","TREATMENT_NAME")])

summed.filt <- summed[,summed$ncells >= 10]
current <- summed.filt
y <- DGEList(counts(current), samples=colData(current))
keep <- filterByExpr(y, group=current$TREATMENT_NAME)
y <- y[keep,]
y <- calcNormFactors(y)

design <- model.matrix(~0+factor(TREATMENT_NAME), y$samples)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
colnames(design) <- c("Saline","LPS")
my.contrasts <- makeContrasts(LPSvsSaline=LPS-Saline, levels = design)

res <- glmQLFTest(fit, coef=ncol(design),contrast = my.contrasts)

res <- topTags(res,n=60000,p.value = 2)
res <- res$table

res <- res %>% mutate(sig = ifelse((FDR <= 0.05 | FDR >= 0.05) & abs(logFC) >= 1, "yes", "no"))


res$diffexpressed <- "unchanged"

res$diffexpressed[res$logFC >= 0.5 &
                    res$PValue <=0.05 &
                    res$sig == "yes" ] <- "LPS"

res$diffexpressed[res$logFC <= -0.5 &
                    res$PValue<=0.05 &
                    res$sig=="yes" ] <- "Saline"
lps.de <- res
## score cells -----------
lps.genes <- rownames(lps.de[lps.de$diffexpressed=="LPS",])
score_intersect <- intersect(lps.genes,rownames(data))
scores_by_cell <- colMeans(
  GetAssayData(data,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
data[["LPSCurrentAstroScore"]] <- scores_by_cell

data$ClusterNames <- factor(data$finalClusters)
levels(data$ClusterNames) <-c("Homeostatic","DAA1","DAA2","Synapse")

scvi_coords <- get_scvi_coords(data,data$finalClusters)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))
text <- data$ClusterNames
text_x <- vapply(split(scvi_coords$UMAP1, text), median, FUN.VALUE=0)
text_y <- vapply(split(scvi_coords$UMAP2, text), median, FUN.VALUE=0)
## Plot F -----------------------
png(file.path(dir,"Figure3","E_LPSCurrentAstroScoringUMAP.png"),
    width = 1500,height=1500,res=300)
ggplot(scvi_coords%>%
         arrange(LPSCurrentAstroScore),aes(x=UMAP1, y=UMAP2, colour=LPSCurrentAstroScore)) +
  geom_point() +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) +
  scale_colour_gradientn(limits= c(0,1), oob=squish,
                         colours =  RColorBrewer::brewer.pal(5,"Purples"))
dev.off()

# G. Sankey Scoring -----------------------

library(networkD3)

counts = GetAssayData(object = data,
                      slot = "counts",
                      layer ="RNA")

integrated.astros <- CreateSeuratObject(counts = counts,
                                        meta.data = data@meta.data) # Recreate object for simplicity


integrated.astros <- FindVariableFeatures(integrated.astros,nfeatures = 5000)
integrated.astros <- NormalizeData(integrated.astros)
integrated.astros <- ScaleData(integrated.astros)
integrated.astros <- RunPCA(integrated.astros)
integrated.astros@reductions$umap <- data@reductions$umap

anchors <- FindTransferAnchors(reference = integrated.astros,
                               query = astros, dims = 1:30,
                               reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors,
                            refdata = integrated.astros$finalClusters,
                            dims = 1:30)

colnames(predictions) <- str_c("LPS_", colnames(predictions))
astros <- AddMetaData(astros, metadata = predictions)

## project into integrated reference UMAP ====================================
integrated.astros <- RunUMAP(integrated.astros, dims = 1:30, return.model = TRUE)
integrated.astros@reductions[["umap"]]@cell.embeddings <- x


astros <- MapQuery(anchorset = anchors, reference = integrated.astros, query = astros,
                   #refdata = list(celltype = "celltype"),
                   reference.reduction = "pca", reduction.model = "umap")

scvi_coords <- get_scvi_coords(astros,astros$LPS_predicted.id)
color.codes <- c("dodgerblue","goldenrod1","red", "darkgreen")
zone <- levels(factor(astros$LPS_predicted.id))
astros <- subset(astros,ident=(c(0,1,2,3,4)))
astros$finalClusters <- astros$seurat_clusters
astros$finalClusters <- factor(astros$finalClusters)
levels(astros$finalClusters) <- c("1","2","3","4","5")
links <- data.frame(
  source = astros$finalClusters,
  target = astros$LPS_predicted.id,
  value  = astros$LPS_prediction.score.max,
  group  = astros$TREATMENT_NAME
)
links$target <- factor(links$target)
levels(links$target) <- c("Homeostatic","DAA1","DAA2","Synapse")

nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
my_color <- paste('d3.scaleOrdinal()',
                  '.domain(["wt_wt|LPS", "wt_wt|Saline",',
                  '"1", "2", "4", "3", "5", "Homeostatic",',
                  '"DAA1", "DAA2", "Synapse"])',
                  ' .range(["#00BFC4", "#F8766D","#F8766D", "#B79F00",',
                  '"#00BA38", "#619CFF", "#F564E3", "dodgerblue",',
                  '"gold", "red", "darkgreen" ])')
 ## Plot G ------------------------
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name",
              sinksRight=FALSE, LinkGroup="group",colourScale=my_color,
              nodeWidth=40, fontSize=13, nodePadding=50)