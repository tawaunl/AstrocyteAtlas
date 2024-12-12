# Script to generate all plots in Figure 2
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
source("~/Documents/scHelpers.R")
#install for scDotPlot
#remotes::install_git("ssh://git@ssh.code.roche.com/omni-bioinfo/packages/scDotPlot.git",git = "external")

# Load Data ====================================================================
dir <- "~/Documents/AstrocytePaper"
data <- readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))


# A. UMAP Clusters ============================================================
scvi_coords <- get_scvi_coords(data,data$finalClusters)
color.codes <- c("dodgerblue","goldenrod1","red", "darkgreen")
zone <- levels(factor(data$finalClusters))
png(file.path(dir,"Figure2","A_ClustersUMAP.png"),
    width = 500,height=500)
ggplot(data=scvi_coords , aes(x=UMAP1, y=UMAP2, colour  = finalClusters)) + 
  geom_point(size=1,alpha = 0.6) +
  scale_colour_manual(values=setNames(color.codes, zone))+ theme_nothing()
dev.off()

##  Gfap Violin Plot=================
png(file.path(dir,"Figure2","A_GfapViolin.png"),
    width = 2500,height=1800,res=300)
VlnPlot(data, features = "Gfap", cols = color.codes, group.by = "finalClusters",
        pt.size = 0) + geom_hline(yintercept = 1,linetype='dotted',size=2)
dev.off()


# B Bar Plots-------------------------------
cells <- data.frame(table(data$finalClusters))
colnames(cells) <- c("Cluster","Value")
#### Plot B -------------
png(file.path(dir,"Figure2","B_CellsBarPlot.png"),
    width = 300,height=500)
ggplot(cells, aes(fill=Cluster, y=log10(Value), x= Cluster )) +
  geom_bar(stat="identity",fill= color.codes) + theme_nothing()+ theme(axis.title = element_text(size = 20),
                                                                       axis.title.y = element_text(angle = 90)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14,face="bold"))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14,face="bold")) +
  ylab("Log10 Cell Count") + xlab("Cluster") +
  scale_colour_manual(values=setNames(color.codes, zone))
dev.off()

#C. Study Bar Plot ==============
coldata <- as.data.frame(data@meta.data)
study_cluster_tmp_by_study <- coldata %>%
  group_by(StudyName, finalClusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

cluster1 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==1),]
cluster2 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==2),]
cluster3 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==3),]
cluster4 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==4),]

cluster1$scaled <- cluster1$freq/sum(cluster1$freq)
cluster2$scaled <- cluster2$freq/sum(cluster2$freq)
cluster3$scaled <- cluster3$freq/sum(cluster3$freq)
cluster4$scaled <- cluster4$freq/sum(cluster4$freq)


study_cluster_tmp_by_study_scaled <- rbind(cluster1,cluster2,cluster3,cluster4)
#### Plot C -----------
png(file.path(dir,"Figure2","C_StudyBarPlot.png"),
    width = 550,height=500)
ggplot(study_cluster_tmp_by_study_scaled, aes(fill=StudyName, y=scaled, x=finalClusters)) +
  geom_bar(position="stack", stat="identity", alpha=0.6) + theme_void()+ theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14)) +
  xlab("Clusters") + ylab("Proportion of cells within a study") +
  theme(axis.title.y =  element_text(angle = 90))
dev.off()

# D. Abundance Plots ===========================================================
library(edgeR)
library(RColorBrewer)
library(viridis)
abundances <- table(data$finalClusters,data$Sample) 
abundances <- unclass(abundances) 

extra.info <- data@meta.data[match(colnames(abundances), data$Sample),]
d <- DGEList(abundances, samples=extra.info)
d = calcNormFactors(d)
d= estimateCommonDisp(d, verbose=TRUE)



norm_counts <- as.data.frame(t(d$counts)) 
colnames(norm_counts) <- paste("Cluster ", colnames(norm_counts), sep="")
norm_counts <- norm_counts %>% mutate(Sample = rownames(.))
coldata <- as.data.frame(data@meta.data)
coldata_short <- coldata %>% dplyr::select(Sample,Disease,StudyName,Model) %>% unique


df_long_final <- left_join(norm_counts,coldata_short, by="Sample") 

percentages <- data.frame(matrix(nrow = dim(df_long_final)[1],ncol = dim(df_long_final)[2]))
sums <- data.frame(clustersums=rowSums(df_long_final[,-c(dim(df_long_final)[2]-3,
                                                         dim(df_long_final)[2]-2,
                                                         dim(df_long_final)[2]-1,
                                                         dim(df_long_final)[2])]),ident=df_long_final$Sample)

for (clust in 1:length(levels(factor(data$finalClusters)))) {
  for (sample in 1:length(df_long_final[,clust])) {
    percent <- df_long_final[sample,clust]/sums[sample,1]
    percentages[sample,clust]<- percent+0.0005737235 # value to add to offset zeros
  }
}

#dd col data back in 
percentages[,dim(df_long_final)[2]-3]<- df_long_final[,dim(df_long_final)[2]-3]
percentages[,dim(df_long_final)[2]-2]<- df_long_final[,dim(df_long_final)[2]-2]
percentages[,dim(df_long_final)[2]-1]<- df_long_final[,dim(df_long_final)[2]-1]
percentages[,dim(df_long_final)[2]]<- df_long_final[,dim(df_long_final)[2]]


colnames(percentages) <-  colnames(df_long_final)
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
colnames(df_long) <- c("Sample","DiseaseLabel","StudyName","Model","Cluster","Percent")
colnames(df_transformed) <- c("Sample","DiseaseLabel","StudyName","Model","Cluster","Percent")

level_order <- c("AD_control","AD","MS_control","MS")
df_long$DiseaseLabel <- factor(df_long$DiseaseLabel,levels = level_order) 
df_long$Cluster <- factor(df_long$Cluster,levels =c(paste0("Cluster ",1:(dim(df_long_final)[2]-4))))
df_transformed$DiseaseLabel <- factor(df_transformed$DiseaseLabel,levels = level_order) 
df_transformed$Cluster <- factor(df_transformed$Cluster,levels =c(paste0("Cluster ",1:(dim(df_long_final)[2]-4))))

fill.codes <- c("white","darkgreen","white", "purple")
color.codes <- c("darkgreen","darkgreen","purple","purple")
zone <- levels(factor(df_long$DiseaseLabel))

filter <- c("Zhao 5XFAD","Choi 5XFAD","Shi P301S")

df_transformed1 <- df_transformed[!df_transformed$StudyName%in% filter,]

####Plot D ---------
png(file.path(dir,"Figure2","D_AbundancePlot.png"),
    width = 4500,height=1000,res=300)
ggplot(df_long, aes(x=DiseaseLabel, y=Percent,fill=DiseaseLabel,colour = DiseaseLabel)) + theme_classic() +
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
        axis.text.x = element_blank()) +
  scale_fill_manual(values=setNames(fill.codes, zone))+
  scale_color_manual(values=setNames(color.codes, zone))
dev.off()

## Significant testing ==================
#Sig test on CLR transformed values
sig.clusters <- list()
for(cluster in 1:length(levels(factor(data$finalClusters)))){
  clusterName <- paste0("Cluster ",levels(factor(data$finalClusters))[cluster])
  dataCluster <- df_transformed1[which(df_transformed1$Cluster==clusterName),]
  res <- kruskal.test(Percent ~ DiseaseLabel, data = dataCluster)
  # check if there is any significance
  if(res$p.value < 0.8/length(levels(factor(data$finalClusters)))){
    sig.clusters[[clusterName]] <- res
    # do a Welch's t-test to find which specific groups are different
    res.AD <- wilcox.test(dataCluster$Percent[which(dataCluster$DiseaseLabel=="AD_control")],
                          dataCluster$Percent[which(dataCluster$DiseaseLabel=="AD")])
    res.MS <- wilcox.test(dataCluster$Percent[which(dataCluster$DiseaseLabel=="MS_control")],
                          dataCluster$Percent[which(dataCluster$DiseaseLabel=="MS")])
    if(res.MS$p.value < 0.9){
      sig.clusters[[clusterName]][["MS"]] <- res.MS
    }
    if(res.AD$p.value < 0.9){
      sig.clusters[[clusterName]][["AD"]] <- res.AD
    }
  }
}

write.csv(df_transformed1,"~/Documents/AstrocytePaper/2024-11-26 lucast3/Fig2D_AbundancePlotDataforStatistics_Mouse.csv")

# E. Markers Dot Plot ==========================================================
counts <- round(GetAssayData(data, slot="counts", assay="RNA"))   
counts <- initializeSparseMatrix(counts)
lognorm <- logNormCounts.chan(counts,batch = data$StudyID)
markers <- scoreMarkers.chan(lognorm,
                             groups=data$finalClusters, batch=data$StudyID,lfc=0)
# Get top 5 markers from each cluster 
topmarkers <- c()
for (cluster in 1:length(markers$statistics)) {
  clustermarkers <- data.frame(markers$statistics[[cluster]])
  clustermarkers <- clustermarkers[order(clustermarkers$logFC,decreasing = TRUE),]
  top <- rownames(clustermarkers)[1:5]
  topmarkers <- c(topmarkers,top)
}

topmarkers[5] <- "Nrxn1"
topmarkers[9:10] <- c("Id4","Clu")

plot <- DotPlot(data,features = c(unique(topmarkers)), assay = "RNA",col.min = 0,
                cols = "RdYlBu",group.by = "finalClusters",scale = TRUE)

png(file.path(dir,"Figure2","C_DotPlotMarkers.png"),
    width = 2500,height=2500,res=400)
plot + RotatedAxis() +coord_flip() + xlab('Gene') +  ylab('Cluster') +
  theme(axis.text = element_text(size = 14),axis.title = element_text(size=18,face='bold'))  
dev.off()



# E. Markers Dot Plot ==========================================================
subtypeList <- list("1" = "Homeostatic",
                    "2" = "DAA1",
                    "3" = "DAA2",
                    "4" = "Synapse")
sce <- data %>% 
  Seurat::as.SingleCellExperiment() 

colData(sce) <- data.frame(colData(sce)) %>% 
  dplyr::mutate(Cluster = dplyr::recode_factor(finalClusters,
                                               !!!subtypeList)) %>% DataFrame()

features <- c("Nrxn3", "Nrg3", "Nkain2", "Dlg2", "Meg3", "mt-Co3", "Aldoc", "Mt1",
              "Apoe", "Cst3", "Clu", "Id4", "Mt2", "Id3", "Gfap", "Nrxn1", "Wdr17",
              "Ntm", "Lsamp", "Gpc5") %>% 
  purrr::set_names(rep(rev(subtypeList), each = 5))

rowData(sce)$Marker <- features[match(rownames(sce), features)] |>
  names() %>% 
  factor(levels = unlist(subtypeList))



pdf(file.path(dir,"Figure2","C_Updataed_DotPlotMarkers.pdf"))
sce %>% 
  scDotPlot::scDotPlot(features = features,
                       group = "Cluster",
                       #block = "Sample",
                       scale = TRUE,
                       cluster = FALSE,
                       groupAnno = "Cluster",
                       featureAnno = "Marker",
                       annoColors = list("Cluster" = c("dodgerblue", "goldenrod1",
                                                       "red", "darkgreen"),
                                         "Marker" = c("dodgerblue", "goldenrod1",
                                                      "red", "darkgreen")),
                       featureLegends = FALSE,
                       annoHeight = 0.025,
                       annoWidth = 0.1,fontSize = 14,
                       dotColors=c( "blue","#FFFFBF", "red"))
dev.off()



# F. Pathway Analysis Go and KEGG  =============================================
# Cluster profiler needs to be up to date to get KEGG to work
#remotes::install_github("GuangchuangYu/GOSemSim")
#devtools::install_github("YuLab-SMU/clusterProfiler")

library(clusterProfiler)
GoPlots <- list()
gcSample <- list()
for (cluster in 1:length(markers$statistics)) {
  clustermarkers <- data.frame(markers$statistics[[cluster]])
  clustermarkers <- clustermarkers[order(clustermarkers$logFC,decreasing = TRUE),]
  geneList <- as.numeric(clustermarkers$logFC[which(clustermarkers$logFC > 0)])
  names(geneList) <- rownames(clustermarkers[which(clustermarkers$logFC > 0),])
  
  gene.df <- bitr( names(geneList), fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org.Mm.eg.db::org.Mm.eg.db)
  x <- match(gene.df$SYMBOL,names(geneList))
  geneList1 <- geneList[x]
  names(geneList1) <- gene.df$ENTREZID
  gcSample[[paste0("Cluster_",cluster)]] <- geneList1
  
}

ck.kegg <- compareCluster(geneClusters = gcSample, fun = "gseKEGG",
                          organism = "mmu",nPermSimple = 10000,scoreType = "pos")


ck.GO <- compareCluster(geneClusters = gcSample, fun = "gseGO",nPermSimple = 10000,
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db,ont="BP",scoreType = "pos",eps=0)

p1 <- dotplot(ck.kegg,by="NES") + ggtitle("KEGG Pathways") + 
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=20))
p2 <- dotplot(ck.GO,by="NES") + ggtitle("GO:BP Pathways") + 
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=20))

png(file.path(dir,"Figure2","D_DotPlotPathways.png"),
    width = 3200,height=2750,res=300)
plot_grid(p1,p2)
dev.off()

for(name in names(gcSample)){
  go <- gseGO(gcSample[[name]],OrgDb = org.Mm.eg.db::org.Mm.eg.db,ont="BP")
  write.csv(data.frame(go@result),file = file.path(dir,"Figure2",paste0(name,"_GOresult.csv")))
  edox2 <- pairwise_termsim(go)
  
  GoPlots[[name]] <-  treeplot(edox2, hclust_method = "average",showCategory=20)
}
for(name in names(gcSample)){
  go <- gseKEGG(gcSample[[name]],organism = "mmu",nPermSimple = 10000,scoreType = "pos")
  write.csv(data.frame(go@result),file = file.path(dir,"Figure2",paste0(name,"_KEGGresult.csv")))
  
}
for(name in names(GoPlots)){
  pdf(file.path(dir,"Figure2",paste0(name,"_GODotPlot_tree.pdf")),width=10,height=10)
  print(GoPlots[[name]])
  dev.off()
}


for(name in names(gcSample)){
  go <- gseGO(gcSample[[name]],OrgDb = org.Mm.eg.db::org.Mm.eg.db,ont="BP")
  GoPlots[[name]] <-  dotplot(go,showCategory=5)
}
for(name in names(GoPlots)){
  pdf(file.path(dir,"Figure2",paste0(name,"_GODotPlot.pdf")),width=10,height=10)
  print(GoPlots[[name]])
  dev.off()
}