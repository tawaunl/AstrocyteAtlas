# supplemental Figure #2

# Write Tables and Extract BIC 
library(Seurat)
library(fastTopics)
library(ggplot2)
library(MatrixExtra)
library(grid)
library(gridExtra)

data.dir <- "/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Rscripts/TopicModeling/14topics_tol0.1"
fig.dir <- "~/Documents/AstrocytePaper/Supplemental/SuppFig2"
load("/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Rscripts/TopicModeling/14topics_tol0.1/FitGoM_50K_k14_tol0.1.rda")
data <- readRDS("/gstore/data/astroMetaAnalysis/data/DS_50k.rds")
usage <- read.csv(file.path(data.dir,"usage_50k.csv"))
theta <- read.csv(file.path(data.dir,"theta_50k.csv"))
scores <- read.csv(file.path(data.dir,"score_min_50k.csv"))

colnames(usage)[1] <- "cell_id"
colnames(theta)[1] <- "gene"
colnames(scores)[1] <- "topic"
n_genes = 25
n_features = 200
prolif_topic = 'lda_14'
prolif_cutoff = 0.08
topic_range <- paste0("lda_",c(1,2,3,5,6,9,12,14))
# join Usage to Seurat object
for (topic in topic_range) {
  data[[paste0("lda_",topic)]] <- usage[[paste0("lda_",topic)]]
  
}

# Prepare list for scores
topic_score_dict = list()
plot_score_dict = list()

for (topic in topic_range) {
  t <- strsplit(topic,"_")[[1]][2]
  topic_scores = scores[t,]
  
  topic_dict = list()
  for (i in 1:n_features) {
    gene_index <- as.numeric(topic_scores[paste0("indices.",i)])
    gene <- theta$gene[gene_index]
    score = as.numeric(topic_scores[paste0("scores.",i)])
    topic_dict[gene] = score
  }
  topic_score_dict[[topic]] = topic_dict
  
  plot_series = lapply(topic_dict,sort,decreasing=TRUE)
  plot_score_dict[[topic]] = plot_series
}

#Plot
topic_plots <- list()

for(topic in topic_range){
  df <- data.frame(gene=names(topic_score_dict[[topic]]),score=unlist(topic_score_dict[[topic]]))
  topic_plots[[topic]] <- ggplot(df[1:10,],aes(x=reorder(gene,score),y=score)) + geom_bar(stat="identity") +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0),labels = function(x) format(x, scientific = TRUE),breaks = range(df$score)) + 
    xlab("Gene") + ggtitle(topic) +
    theme_light()+ theme(plot.title = element_text(hjust = 0.5)) 
  x=x+1
}

pdf(file.path(fig.dir,"TopGenesinTopics_condensed.pdf"),width = 12,height=8)
do.call("grid.arrange",list(grobs=topic_plots,top=textGrob("Gene scores for each topic"),ncol=4))
dev.off()
# write gene loadings to CSV -------------
write.csv(theta,file = file.path(fig.dir,"TopicGeneLoadings.csv"))
write.csv(usage,file = file.path(fig.dir,"TopicCellScores.csv"))

topic_plots_umap <- list()

for(topic in topic_range){
  
  topic_plots_umap[[topic]] <- FeaturePlot(data,topic,pt.size = 1.5) +
    ggtitle(topic) +
    theme_void() + theme(plot.title = element_text(hjust = 0.5,face = "bold"), legend.position = "none") 
  
}

png(file.path(fig.dir,"UMAPscoringTopics_condensed.png"),
    width = 1000, height=500)
do.call("grid.arrange",list(grobs=topic_plots_umap,
                            top=textGrob(expression(bold(underline("UMAP embedding for each topic")))),ncol=4))
dev.off()


# GSEA on top 50 Topic loadings
library(clusterProfiler)
# First get top 50 loadings from each topic
geneList <- list()
for(topic in 2:15){
  lda_data <- theta[,c(1,topic)]
  lda_data <- lda_data[order(lda_data[,2],decreasing = TRUE),]
  genes <- lda_data[,2]
  names(genes) <-lda_data[,1]
  genes <- genes[order(genes,decreasing = T)]
  geneList[[colnames(lda_data)[2]]] <- genes
}


ck <- lapply(geneList, function(x){
  gseGO(x,OrgDb = org.Mm.eg.db::org.Mm.eg.db,
        ont="BP",keyType="SYMBOL",eps=0,
        pvalueCutoff = 0.1)
})

saveRDS(ck, file.path(data.dir,"GSEA_results.rds"))
ck <- readRDS(file.path(data.dir,"GSEA_results.rds"))

gsea_results <- lapply(names(ck), function(cluster){
  ck[[cluster]]@result
})
names(gsea_results) <- names(ck)
# Assuming gsea_results is a list of GSEA results per cluster
gsea_results <- gsea_results[topic_range]
source("~/Documents/scHelpers.R")

pdf(file.path(fig.dir,"GSEAonTopics_condensed.pdf"),width=8,height = 8)
DotPlotCompare(
  gsea_list = gsea_results,
  n = 3,
  size_col = "NES",
  color_col = "p.adjust",
  size_cutoff = NULL,      # Filter to pathways with NES >= 1.5
  color_cutoff = 0.01,
  direction = "positive" # Filter to pathways with p.adjust <= 0.05
)
dev.off()

# DAA subclusters ------------

## UMAP ----------
DAA1<- readRDS("~/Documents/AstrocytePaper/DAA1subcluster.rds")
scvi_coords <- get_scvi_coords(DAA1,DAA1$subclusters)
library(cowplot)


cluster_centroids <- aggregate(cbind(UMAP1, UMAP2) ~ subclusters, scvi_coords, mean)
scvi_coords$subclusters <- factor(scvi_coords$subclusters,levels = c("Myoc+","Reactive1",
                                                                    "Reactive2","Interferon-Responsive",
                                                                    "Meg3+"))

cairo_pdf(file.path(fig.dir,"DAA1_subClustersUMAP.pdf"),
    width = 10,height=10)
ggplot(data=scvi_coords %>% arrange(subclusters) , aes(x=UMAP1, y=UMAP2, colour  = subclusters)) + 
  ggrastr::geom_point_rast(size=2,alpha = 0.9) +
  geom_text(
    data = cluster_centroids,                         # Add cluster labels at centroids
    aes(x = UMAP1, y = UMAP2, label = subclusters),
    size = 6, fontface = "bold",family="Arial", color = "black"
  ) + theme_nothing()
dev.off()

## Markers ---------
library(scran.chan)

counts <- round(GetAssayData(DAA1, slot="counts", assay="RNA"))   
counts <- initializeSparseMatrix(counts)
lognorm <- logNormCounts.chan(counts,batch = DAA1$StudyID)
markers <- scoreMarkers.chan(lognorm,
                             groups=DAA1$subclusters, batch=DAA1$StudyID,lfc=0)
saveRDS(markers,"~/Documents/AstrocytePaper/DAA1subclustermarkers.rds")
# Get top 5 markers from each cluster 
topmarkers <- c()
for (cluster in 1:length(markers$statistics)) {
  clustermarkers <- data.frame(markers$statistics[[cluster]])
  clustermarkers <- clustermarkers[order(clustermarkers$logFC,decreasing = TRUE),]
  top <- rownames(clustermarkers)[1:15]
  topmarkers <- c(topmarkers,top)
}


sce <- DAA1 %>% 
  Seurat::as.SingleCellExperiment() 

library(scran)
library(purrr)
library(dplyr)
library(AnnotationDbi)


sce$subclusters <- factor(sce$subclusters,levels = c("Myoc+","Reactive1",
                                                     "Reactive2","Interferon-Responsive",
                                                     "Meg3+"))

features <- c("Myoc","Id1","Ctnna2","Wdr17","Gpm6b",
              "Acsl3","Aqp4","Ckb","Kcnn2","Gja1",
              "Mdga2","Ncam2","Ank2","Sorbs1","C4b",
              "Bst2","Ifit3","Oasl2","B2m","Stat1",
              "Meg3","Nrg3","Opcml","Ptprd","Lrrtm4")
features <- setNames(features,rep(levels(sce$subclusters), each = 5))



out_for_plotting <- sce[features,]

rowData(out_for_plotting)$Marker <- features[match(rownames(out_for_plotting), features)] |>
  names() %>% 
  factor(levels =levels(sce$subclusters))

cairo_pdf(file.path(fig.dir,"DAA1subClustersMarkerDotPlot.pdf"))
out_for_plotting %>% 
  scDotPlot::scDotPlot(features = unique(features),
                       group = "subclusters",
                       #block = "Sample",
                       scale = TRUE,
                       cluster = FALSE,
                       groupAnno = "subclusters",
                       featureAnno = "Marker",
                       featureLegends = FALSE,
                       annoHeight = 0.025,
                       annoColors = list("subclusters" = scales::hue_pal()(5),
                                         "Marker" = scales::hue_pal()(5)),
                       annoWidth = 0.05,fontSize = 20,
                       dotColors=c( "blue","#FFFFBF", "red")) 
dev.off()
## Cellularity ---------------------------
dir <- "~/Documents/AstrocytePaper"
data <- readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))
calc.col = "subclusters"
subset.cols = c("Sample","Disease","BrainRegion","StudyName")
samplecol = "Sample"
library(edgeR)
abundances<- table(DAA1$subclusters,factor(DAA1$Sample))

abundances <- unclass(abundances)

extra.info <- DAA1@meta.data[match(colnames(abundances),DAA1@meta.data[[samplecol]]),]
d <- DGEList(abundances, samples=extra.info)

d = calcNormFactors(d)
d= estimateCommonDisp(d, verbose=TRUE)

norm_counts <- as.data.frame(t(d$counts)) 
norm_counts <- norm_counts %>% mutate( !! samplecol := rownames(.))
coldata <- as.data.frame(DAA1@meta.data)
coldata_short <- coldata %>% dplyr::select(all_of(subset.cols)) %>% unique

df_long_final <- left_join(norm_counts,coldata_short, by=samplecol)
totsums <- data.frame(table(data$Sample))
#totsums <- totsums %>% filter(Var1 == "DAA1")
exclude <- (dim(df_long_final)[2]-dim(coldata_short)[2]+1):dim(df_long_final)[2]
rownames(totsums) <- totsums$Var1
df_long_final[[samplecol]] <- factor(df_long_final[[samplecol]])
percentages <- data.frame(matrix(nrow = dim(df_long_final)[1],ncol = dim(df_long_final)[2]))

sums <- data.frame(clustersums=rowSums(df_long_final[,-exclude]),ident=df_long_final[[samplecol]])
totsums <- totsums[match(df_long_final$Sample,totsums$Var1),]
for (clust in 1:length(levels(factor(coldata[[calc.col]])))) {
  for (sample in 1:length(df_long_final[,clust])) {
    sam <- df_long_final[sample,6]
    percent <- (df_long_final[sample,clust]/totsums[sample,2])
    percentages[sample,clust]<- percent
  }
}


percentages[,exclude] <- df_long_final[,exclude]

colnames(percentages) <-  colnames(df_long_final)


df_long <- percentages %>% 
  pivot_longer(unique(coldata[[calc.col]]),values_to = "Percent",names_to = "Cluster")


level_order <- c("AD_control","AD","MS_control","MS")
df_long$Disease <- factor(df_long$Disease,levels = level_order) 

fill.codes <- c("white","darkgreen","white", "purple")
color.codes <- c("darkgreen","darkgreen","purple","purple")
zone <- levels(factor(df_long$Disease))

filter <- c("Zhao 5XFAD","Choi 5XFAD","Shi P301S")

df_long <- df_long[! df_long$StudyName %in% filter,]

cairo_pdf(file.path(fig.dir, "DAA1subclusterCellularity.pdf"),width=10)
ggplot(df_long, aes(x=Disease, y=Percent,fill=Disease,colour = Disease))+ theme_classic() +
  geom_jitter(color="darkgrey",width = 0.4,size=3,alpha=0.6)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  facet_wrap(~Cluster,scales = "free",ncol=3) + xlab("Disease Label") +
  theme(axis.text.y = element_text(size=20)) +
  ylab("Cellularity Proportion") +  
  theme(strip.text.x = element_text(size = 16,face = "bold"),
        axis.title=element_text(size=12,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(values=setNames(fill.codes, zone))+
  scale_color_manual(values=setNames(color.codes, zone))
dev.off()


sig.clusters <- list()
for(cluster in 1:length(levels(factor(DAA1$subclusters)))){
  clusterName <- levels(factor(DAA1$subclusters))[cluster]
  dataCluster <- df_transformed[which(df_transformed$Cluster==clusterName),]
  res <- kruskal.test(Percent ~ Disease, data = dataCluster)
  # check if there is any significance
  if(res$p.value < 0.8/length(levels(factor(DAA1$subclusters)))){
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



#DAA2 -----------
DAA2 <- readRDS("~/Documents/AstrocytePaper/DAA2subcluster.rds")


## UMAP ----------
scvi_coords <- get_scvi_coords(DAA2,DAA2$subclusters)
library(cowplot)


cluster_centroids <- aggregate(cbind(UMAP1, UMAP2) ~ subclusters, scvi_coords, mean)
scvi_coords$subclusters <- factor(scvi_coords$subclusters,levels = c("NT-Transport","Fth1+",
                                                                     "OxPhos","Gfap+"))

cairo_pdf(file.path(fig.dir,"DAA2_subClustersUMAP.pdf"),
          width = 10,height=10)
ggplot(data=scvi_coords %>% arrange(subclusters) , aes(x=UMAP1, y=UMAP2, colour  = subclusters)) + 
  ggrastr::geom_point_rast(size=3,alpha = 0.9) +
  geom_text(
    data = cluster_centroids,                         # Add cluster labels at centroids
    aes(x = UMAP1, y = UMAP2, label = subclusters),
    size = 12, fontface = "bold",family="Arial", color = "black"
  ) + theme_nothing()
dev.off()

## Markers -------------
sce <- DAA2 %>% 
  Seurat::as.SingleCellExperiment() 
sce$subclusters <- factor(scvi_coords$subclusters,levels = c("NT-Transport","Fth1+",
                                                             "OxPhos","Gfap+"))

features <- c("Atp1a2","Slc4a4","Slc6a1","Atp1b2","Slc6a11",
              "Fth1","Cpe","Mt3","Gpm6b","Cst3",
              "Eef1b2","H3f3b","Clta","Cbr3","Prex2",
              "Tmsb4x","Gfap","C4b","Pcsk1n","Celf4")

features <- setNames(features,rep(levels(sce$subclusters), each = 5))



out_for_plotting <- sce[features,]
out_for_plotting$subclusters <- factor(out_for_plotting$subclusters)
rowData(out_for_plotting)$Marker <- features[match(rownames(out_for_plotting), features)] |>
  names() %>% 
  factor(levels =levels(out_for_plotting$subclusters))


cairo_pdf(file.path(fig.dir,"DAA2subClustersMarkerDotPlot.pdf"))
out_for_plotting %>% 
  scDotPlot::scDotPlot(features = unique(features),
                       group = "subclusters",
                       #block = "Sample",
                       scale = TRUE,
                       cluster = FALSE,
                       groupAnno = "subclusters",
                       featureAnno = "Marker",
                       featureLegends = FALSE,
                       annoHeight = 0.025,
                       annoColors = list("subclusters" = scales::hue_pal()(4),
                                         "Marker" = scales::hue_pal()(4)),
                       annoWidth = 0.05,fontSize = 20,
                       dotColors=c( "blue","#FFFFBF", "red")) 
dev.off()
## Cellularity ---------------
calc.col = "subclusters"
subset.cols = c("Sample","Disease","BrainRegion","StudyName")
samplecol = "Sample"
library(edgeR)
abundances<- table(DAA2$subclusters,factor(DAA2$Sample))

abundances <- unclass(abundances)

extra.info <- DAA2@meta.data[match(colnames(abundances),DAA2@meta.data[[samplecol]]),]
d <- DGEList(abundances, samples=extra.info)

d = calcNormFactors(d)
d= estimateCommonDisp(d, verbose=TRUE)

norm_counts <- as.data.frame(t(d$counts)) 
norm_counts <- norm_counts %>% mutate( !! samplecol := rownames(.))
coldata <- as.data.frame(DAA2@meta.data)
coldata_short <- coldata %>% dplyr::select(all_of(subset.cols)) %>% unique

df_long_final <- left_join(norm_counts,coldata_short, by=samplecol)
totsums <- data.frame(table(data$Sample))
exclude <- (dim(df_long_final)[2]-dim(coldata_short)[2]+1):dim(df_long_final)[2]
rownames(totsums) <- totsums$Var1
df_long_final[[samplecol]] <- factor(df_long_final[[samplecol]])
percentages <- data.frame(matrix(nrow = dim(df_long_final)[1],ncol = dim(df_long_final)[2]))

sums <- data.frame(clustersums=rowSums(df_long_final[,-exclude]),ident=df_long_final[[samplecol]])
totsums <- totsums[match(df_long_final$Sample,totsums$Var1),]
for (clust in 1:length(levels(factor(coldata[[calc.col]])))) {
  for (sample in 1:length(df_long_final[,clust])) {
    sam <- df_long_final[sample,6]
    percent <- (df_long_final[sample,clust]/totsums[sample,2])
    percentages[sample,clust]<- percent
  }
}


percentages[,exclude] <- df_long_final[,exclude]

colnames(percentages) <-  colnames(df_long_final)


df_long <- percentages %>% 
  pivot_longer(unique(coldata[[calc.col]]),values_to = "Percent",names_to = "Cluster")


level_order <- c("AD_control","AD","MS_control","MS")
df_long$Disease <- factor(df_long$Disease,levels = level_order) 

fill.codes <- c("white","darkgreen","white", "purple")
color.codes <- c("darkgreen","darkgreen","purple","purple")
zone <- levels(factor(df_long$Disease))

filter <- c("Zhao 5XFAD","Choi 5XFAD","Shi P301S")

df_long <- df_long[! df_long$StudyName %in% filter,]

df_long$Cluster <- factor(df_long$Cluster,levels = c("NT-Transport","Fth1+",
                                                             "OxPhos","Gfap+"))
cairo_pdf(file.path(fig.dir, "DAA2subclusterCellularity.pdf"),height = 6,width = 13)
ggplot(df_long, aes(x=Disease, y=Percent,fill=Disease,colour = Disease))+ theme_classic() +
  geom_jitter(color="darkgrey",width = 0.4,size=3,alpha=0.6)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  facet_wrap(~Cluster,scales = "free",ncol=4) + xlab("Disease Label") +
  theme(axis.text.y = element_text(size=20)) +
  ylab("Cellularity Proportion") +  
  theme(strip.text.x = element_text(size = 16,face = "bold"),
        axis.title=element_text(size=12,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(values=setNames(fill.codes, zone))+
  scale_color_manual(values=setNames(color.codes, zone))
dev.off()
