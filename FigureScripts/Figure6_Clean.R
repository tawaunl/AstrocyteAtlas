# Figure 6

# Load libraries -----------------
library(ggplot2)
library(clusterProfiler)
library(rstatix)
library(ggrepel)


# Load Data  ----------------------------
source("~/Documents/scHelpers.R")
dir <- "~/Documents/AstrocytePaper/Human"
human.AD <- readRDS(file.path(dir,"res_dsl.AD.cain_fixed.072924.rds"))

human.MS <- readRDS(file.path(dir,"res_dsl.MS.080124.rds"))
human.PD <- readRDS(file.path(dir,"res_dsl.PD.072924.rds"))
fig.dir <- "~/Documents/AstrocytePaper/Figure6"
# Process Data -------------------------
feat <- readRDS(file.path(dir,"humanFeatures.rds"))

human.AD <- generatePlotTable(human.AD)
human.MS <- generatePlotTable(human.MS)
human.PD <- generatePlotTable(human.PD)

human.AD$diffexpressed <- ifelse(human.AD$FDR<=0.05 & human.AD$dl_mu>=0.5 & human.AD$n_up>=3, "up",
                                 ifelse(human.AD$FDR<=0.05 & human.AD$dl_mu<= -0.5 & human.AD$n_down>=3, "down","unchanged"))
human.MS$diffexpressed <- ifelse(human.MS$FDR<=0.05 & human.MS$dl_mu>=0.5 & human.MS$n_up>2, "up",
                                 ifelse(human.MS$FDR<=0.05 & human.MS$dl_mu<= -0.5 & human.MS$n_down>=3, "down","unchanged"))
human.PD$diffexpressed <- ifelse(human.PD$FDR<=0.05 & human.PD$dl_mu>=0.5 & human.PD$n_up>=2, "up",
                                 ifelse(human.PD$FDR<=0.05 & human.PD$dl_mu<= -0.5 & human.PD$n_down>=2, "down","unchanged"))
protein_coding <- feat %>% dplyr::filter(type == "protein_coding")
human.AD <- human.AD %>% dplyr::filter(symbol %in% protein_coding$symbol)
human.MS <- human.MS %>% dplyr::filter(symbol %in% protein_coding$symbol)
human.PD <- human.PD %>% dplyr::filter(symbol %in% protein_coding$symbol)

positive_logFC_thresh <- 1.3 # Threshold for upregulated genes
negative_logFC_thresh <- -1 # Threshold for downregulated genes
pvalue_thresh <- 1e-3         

# A. Diseae Volcano Plots ------------------
## AD -------------
genes <- c("CP","SLC5A3","COL27A1","COL8A1","FBXO2","FBXO32","S100A6","MT1F","C3")

human.AD$label <- with(human.AD, 
                     ((dl_mu > positive_logFC_thresh & PValue <= pvalue_thresh) | 
                        (dl_mu < negative_logFC_thresh & PValue <= pvalue_thresh) |
                        symbol %in% genes | (dl_mu > 0.5 & PValue < 1e-10)))
AD_volcano <- ggplot(data=human.AD, aes(x=dl_mu, y=-log10(PValue), col=diffexpressed, label=symbol)) + 
  ggrastr::geom_point_rast(alpha=0.6,size=3.5) + 
  theme_classic() +
  geom_text_repel(size=8,data = subset(human.AD,label),show.legend = FALSE,box.padding = 0.5,
                  max.overlaps = 10) + 
  scale_color_manual(values=c("blue", "grey", "firebrick3"), name = "")+
  theme(legend.position = "none",
        axis.title = element_text(size=20),
        plot.title = element_text(size=28,hjust = 0.5)) +
  xlab("Meta-LogFoldChange") + ylab(bquote(bold(-log[10](Meta-PValue))))


## MS ---------------
human.MS$label <- with(human.MS, 
                       ((dl_mu > positive_logFC_thresh & PValue <= pvalue_thresh) | 
                          (dl_mu < -0.8 & PValue <= pvalue_thresh) |
                          symbol %in% genes | (dl_mu > 0.8 & PValue < 1e-6)))
MS_volcano <- ggplot(data=human.MS, aes(x=dl_mu, y=-log10(PValue), col=diffexpressed, label=symbol)) + 
  ggrastr::geom_point_rast(alpha=0.6,size=3.5) + 
  theme_classic() +
  geom_text_repel(size=8,data = subset(human.MS,label & diffexpressed %in% c("up","down")),show.legend = FALSE,box.padding = 0.5,
                  max.overlaps = 10) + 
  scale_color_manual(values=c("blue", "grey", "firebrick3"), name = "")+
  theme(legend.position = "none",
        axis.title = element_text(size=20),
        plot.title = element_text(size=28,hjust = 0.5)) +
  xlab("Meta-LogFoldChange") + ylab(bquote(bold(-log[10](Meta-PValue))))

## PD ---------------
pvalue_thresh <- 1e-5
human.PD$label <- with(human.PD, 
                       ((dl_mu > positive_logFC_thresh & PValue <= pvalue_thresh) | 
                          (dl_mu < negative_logFC_thresh & PValue <= pvalue_thresh) |
                          symbol %in% genes | (abs(dl_mu) > 0.8 & PValue < 1e-6)))
PD_volcano <- ggplot(data=human.PD, aes(x=dl_mu, y=-log10(PValue), col=diffexpressed, label=symbol)) + 
  ggrastr::geom_point_rast(alpha=0.6,size=3.5) + 
  theme_classic() +
  geom_text_repel(size=8,
                  data = subset(human.PD,label & diffexpressed %in% c("up","down")),
                  show.legend = FALSE, box.padding = 0.5,
                  max.overlaps = 15) + 
  scale_color_manual(values=c("blue", "grey", "firebrick3"), name = "")+
  theme(legend.position = "none",
        axis.title = element_text(size=20),
        plot.title = element_text(size=28,hjust = 0.5)) +
  xlab("Meta-LogFoldChange") + ylab(bquote(bold(-log[10](Meta-PValue))))
### Plot A --------------------
cairo_pdf(file.path(fig.dir,"A_VolcanoPlots.pdf"),
    width = 30,height = 10)
cowplot::plot_grid(AD_volcano,MS_volcano,PD_volcano,ncol=3) 
dev.off()

#B. DE gene up Scores by Study ------------------
se_agg <-  readRDS(file.path(data.dir,"aggregated_cells.AD_MS_PD.070724.DonoruniqueID.SE.rds"))

## AD------------------------
up_candidates <- human.AD %>% dplyr::filter(diffexpressed == "up")
down_candidates <- human.AD %>% dplyr::filter(diffexpressed == "down")

se_agg$AD_up_score <- colMeans(assay(se_agg, 'logCPM')[up_candidates$ID,],na.rm = TRUE)
se_agg$AD_down_score <- colMeans(assay(se_agg, 'logCPM')[down_candidates$ID,],na.rm = TRUE)

se_agg$diagnosis_harmonized_by_disease <- 
  ifelse(se_agg$studybatch %in% c("jakel", "schirmer", "absinta", "bryois") &
           se_agg$diagnosis_harmonized == "Control", "MS_Control",
         
         ifelse(se_agg$studybatch %in% c("smajic","wang", "kamath") &
                  se_agg$diagnosis_harmonized == "Control", "PD_Control",
                
                ifelse(se_agg$studybatch %in% c("cain", "morabito","gerrits", "smith", "SEA-AD", "liddelow") &
                         se_agg$diagnosis_harmonized == "Control", "AD_Control",
                       se_agg$diagnosis_harmonized)))

se_agg$diagnosis_harmonized_by_disease <- ifelse(se_agg$diagnosis_harmonized == "RRMS", "MS",
                                                 se_agg$diagnosis_harmonized_by_disease)

coldata_se_agg <- as.data.frame(colData(se_agg)) %>%
  dplyr::filter(diagnosis_harmonized_by_disease %in% c("AD", "AD_Control"))

coldata_se_agg$diagnosis_harmonized <- factor(coldata_se_agg$diagnosis_harmonized, c("Control","AD"))

fill.codes <- c("white","darkgreen","white", "purple")
color.codes <- c("darkgreen","darkgreen","purple","purple")
zone <- levels(coldata_se_agg$diagnosis_harmonized)
p_AD_up_score_by_study <- ggplot(coldata_se_agg,
                                 aes(x=diagnosis_harmonized,
                                     y=AD_up_score,
                                     fill=diagnosis_harmonized,
                                     color =diagnosis_harmonized )) +
  theme_classic() +  geom_jitter(color="darkgrey",width = 0.1,size=3,alpha=0.6)+ 
  geom_boxplot(outlier.shape=NA,alpha=0.5) +
  facet_wrap(~studybatch, ncol= 3) +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.5)) + 
  scale_fill_manual(values=setNames(fill.codes, zone), name = "Diagnosis") +
  scale_color_manual(values=setNames(color.codes, zone)) +
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())
### Plot AD ----------
pdf(file.path(fig.dir,"Human_AD_DEG_up_score_per_study.pdf"))
p_AD_up_score_by_study
dev.off()

## MS ------------------

up_candidates <- human.MS %>% dplyr::filter(diffexpressed == "up")
down_candidates <- human.MS %>% dplyr::filter(diffexpressed == "down")

se_agg$MS_up_score <- colMeans(assay(se_agg, 'logCPM')[up_candidates$ID,],na.rm = TRUE)
se_agg$MS_down_score <- colMeans(assay(se_agg, 'logCPM')[down_candidates$ID,],na.rm = TRUE)

coldata_se_agg <- as.data.frame(colData(se_agg)) %>% dplyr::filter(diagnosis_harmonized_by_disease %in% c("MS", "MS_Control"))
coldata_se_agg$diagnosis_harmonized <- factor(coldata_se_agg$diagnosis_harmonized, c("Control","MS"))

fill.codes <- c("white", "purple")
color.codes <- c("purple","purple")
zone <- levels(coldata_se_agg$diagnosis_harmonized)
p_MS_up_score_by_study <- ggplot(coldata_se_agg,
                                 aes(x=diagnosis_harmonized, y=MS_up_score,
                                     fill=diagnosis_harmonized, color =diagnosis_harmonized )) +
  theme_classic() +  geom_jitter(color="darkgrey",width = 0.1,size=3,alpha=0.6)+ 
  geom_boxplot(outlier.shape=NA,alpha=0.5) +
  facet_wrap(~studybatch, ncol= 2) +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.5)) + 
  scale_fill_manual(values=setNames(fill.codes, zone), name = "Diagnosis") +
  scale_color_manual(values=setNames(color.codes, zone)) +
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

### Plot MS ----------------
pdf(file.path(fig.dir,"Figure6","Human_MS_DEG_up_score_per_study.pdf"))
p_MS_up_score_by_study
dev.off()

## PD ------------------

up_candidates <- human.PD %>% dplyr::filter(diffexpressed == "up")
down_candidates <- human.PD %>% dplyr::filter(diffexpressed == "down")

se_agg$PD_up_score <- colMeans(assay(se_agg, 'logCPM')[up_candidates$ID,],na.rm = TRUE)
se_agg$PD_down_score <- colMeans(assay(se_agg, 'logCPM')[down_candidates$ID,],na.rm = TRUE)

coldata_se_agg <- as.data.frame(colData(se_agg)) %>% dplyr::filter(diagnosis_harmonized_by_disease %in% c("PD", "PD_Control"))
coldata_se_agg$diagnosis_harmonized <- factor(coldata_se_agg$diagnosis_harmonized, c("Control","PD"))

fill.codes <- c("white", "mediumblue")
color.codes <- c("mediumblue","mediumblue")
zone <- levels(coldata_se_agg$diagnosis_harmonized)
p_PD_up_score_by_study <- ggplot(coldata_se_agg,
                                 aes(x=diagnosis_harmonized, y=PD_up_score,
                                     fill=diagnosis_harmonized, color =diagnosis_harmonized )) +
  theme_classic() +  geom_jitter(color="darkgrey",width = 0.1,size=3,alpha=0.6)+ 
  geom_boxplot(outlier.shape=NA,alpha=0.5) +
  facet_wrap(~studybatch, ncol= 2) +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.5)) + 
  scale_fill_manual(values=setNames(fill.codes, zone), name = "Diagnosis") +
  scale_color_manual(values=setNames(color.codes, zone)) +
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())
###Plot PD -------------
pdf(file.path(fig.dir,"Human_PD_DEG_up_score_per_study.pdf"))
p_PD_up_score_by_study
dev.off()


# C. Compared Pathways dot Plot ---------------
library(clusterProfiler)
## Set-up -------
###AD -------
human.AD <- human.AD[order(human.AD$dl_mu,decreasing = TRUE),]
AD.genes <- human.AD$dl_mu
names(AD.genes) <- human.AD$ID

gene.df <- bitr( names(AD.genes), fromType = "ENSEMBL",
                 toType = c("SYMBOL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db::org.Hs.eg.db)
x <- match(gene.df$ENSEMBL,names(AD.genes))
AD.gsea <- AD.genes[x]
names(AD.gsea) <- gene.df$ENTREZID

AD.gsea <- AD.gsea[-which(is.na(AD.gsea)==TRUE)]
###MS --------
human.MS <- human.MS[order(human.MS$dl_mu,decreasing = TRUE),]
MS.genes <- human.MS$dl_mu
names(MS.genes) <- human.MS$ID

gene.df <- bitr( names(MS.genes), fromType = "ENSEMBL",
                 toType = c("SYMBOL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db::org.Hs.eg.db)
x <- match(gene.df$ENSEMBL,names(MS.genes))
MS.gsea <- MS.genes[x]
names(MS.gsea) <- gene.df$ENTREZID
MS.gsea <- MS.gsea[-which(is.na(MS.gsea)==TRUE)]

### PD ------------
human.PD <- human.PD[order(human.PD$dl_mu,decreasing = TRUE),]
PD.genes <- human.PD$dl_mu
names(PD.genes) <- human.PD$ID

gene.df <- bitr( names(PD.genes), fromType = "ENSEMBL",
                 toType = c("SYMBOL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db::org.Hs.eg.db)
x <- match(gene.df$ENSEMBL,names(PD.genes))
PD.gsea <- PD.genes[x]
names(PD.gsea) <- gene.df$ENTREZID
PD.gsea <- PD.gsea[-which(is.na(PD.gsea)==TRUE)]


## Individual Pathways============

## GO Pathways =====================
### AD -------
AD.paths <- gseGO(geneList     = AD.gsea,
                  OrgDb        = org.Hs.eg.db::org.Hs.eg.db,
                  ont          = "All",
                  minGSSize    = 50,
                  maxGSSize    = 500,
                  eps          = 0,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)

### MS --------
MS.paths <- gseGO(geneList     = MS.gsea,
                  OrgDb        = org.Hs.eg.db::org.Hs.eg.db,
                  ont          = "All",
                  minGSSize    = 50,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  eps          =0,
                  verbose      = FALSE)

###PD --------------------------
PD.paths <- gseGO(geneList     = PD.gsea,
                  OrgDb        = org.Hs.eg.db::org.Hs.eg.db,
                  ont          = "All",
                  minGSSize    = 50,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  eps          = 0,
                  verbose      = FALSE)

##Save Pathways--------
DiseasePaths <- list(AD_GSEAGO=AD.paths ,MS_GSEAGO =MS.paths,PD_GSEAGO=PD.paths)
saveRDS(DiseasePaths,file.path(dir,"HumanDiseasePathwaysGO.rds"))


## Grouped Pathways =================================
ck <- compareCluster(list(AD=AD.gsea,MS=MS.gsea,PD=PD.gsea),
                     fun = "gseGO",OrgDb = org.Hs.eg.db::org.Hs.eg.db ,ont="All",
                     eps=0)

ck1 <- filter(ck,NES>1,)
ck2 <- filter(ck,NES<1)

go_up <- dotplot(ck1, by="NES") + ggtitle("GO Pathways Upregulated") + 
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=20)) +xlab("Disease")
go_down <- dotplot(ck2, by="NES") + ggtitle("GO Pathways Downregulated") + 
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=20))+xlab("Disease")

###Plot C ---------------------
pdf(file.path(fig.dir,"DiseasePathwasyDotPlot.pdf"))
cowplot::plot_grid(go_up,go_down,align = "v",ncol= 2)
dev.off()

#D. Abundance Plot -------------------------------
library(edgeR)
library(RColorBrewer)
library(tidyverse)

coldata <- readRDS(file.path(dir,
                             "AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.COLDATA.correct_diagnosis_harmonized_by_disease.Cain_donors_subsetted.rds"))

abundances <- table(coldata$clusters_res_0.5,coldata$donorIdUnique)
abundances <- unclass(abundances) 

extra.info <- coldata[match(colnames(abundances), coldata$donorIdUnique),]
d <- DGEList(abundances, samples=extra.info)
d = calcNormFactors(d)
d= estimateCommonDisp(d, verbose=TRUE)

norm_counts <- as.data.frame(t(d$counts)) 
colnames(norm_counts) <- paste("Cluster ", colnames(norm_counts), sep="")
norm_counts <- norm_counts %>% mutate(donorIdUnique = rownames(.))

coldata_short <- as.data.frame(coldata) %>% dplyr::select(donorIdUnique,diagnosis_harmonized_by_disease,studybatch,apoe) %>% unique


df_long_final <- left_join(norm_counts,coldata_short, by="donorIdUnique") 
percentages <- data.frame(matrix(nrow = dim(df_long_final)[1],ncol = dim(df_long_final)[2]))
sums <- data.frame(clustersums=rowSums(df_long_final[,1:dim(norm_counts)[2]-1]),ident=df_long_final$donorIdUnique)



for (clust in 1:length(levels(factor(coldata$clusters_res_0.5)))) {
  for (sample in 1:length(df_long_final[,clust])) {
    percent <- (df_long_final[sample,clust]/sums[sample,1])*100
    percentages[sample,clust]<- percent
  }
}


cols2replace <- dim(norm_counts)[2]:(dim(norm_counts)[2]+dim(coldata_short)[2]-1)
percentages[,cols2replace]<- df_long_final[,cols2replace]



colnames(percentages) <-  colnames(df_long_final)
df_long <- percentages %>% 
  pivot_longer(c(c(paste0("Cluster ",1:(dim(norm_counts)[2]-1)))))

level_order <- c("AD_Control","AD","MS_Control","MS","PD_Control","PD")
df_long$diagnosis_harmonized_by_disease <- factor(df_long$diagnosis_harmonized_by_disease,levels = level_order) 

df_long$name <- factor(df_long$name,levels =c(paste0("Cluster ",1:(dim(norm_counts)[2]-1))))

fill.codes <- c("white","darkgreen","white", "purple","white","mediumblue")
color.codes <- c("darkgreen","darkgreen","purple","purple","mediumblue","mediumblue")

zone <- levels(factor(df_long$diagnosis_harmonized_by_disease))

###Plot D -----------------------------
png(file.path(fig.dir,"Human_CellularityPlot_with_percentages.png"),
    width = 6000,height=1400,res=300)
ggplot(df_long, aes(x=diagnosis_harmonized_by_disease, y=value,fill=diagnosis_harmonized_by_disease,colour = diagnosis_harmonized_by_disease)) +
  theme_classic() +
  geom_jitter(color="darkgrey",width = 0.2,size=1,alpha=.5)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5)+ 
  
  facet_wrap(~name,scales = "free", ncol=8) + xlab("Disease Label") +
  theme(axis.text.x = element_text(angle=75,vjust = 0.5),axis.text.y = element_text(size=20)) +
  ylab("Cellularity Percent") +  
  theme(strip.text.x = element_text(size = 14,face = "italic"),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50,75))+
  scale_fill_manual(values=setNames(fill.codes, zone))+
  scale_color_manual(values=setNames(color.codes, zone))

dev.off()

## Sig testing --------------
coldata_short <- as.data.frame(coldata) %>% dplyr::select(donorIdUnique,diagnosis_harmonized_by_disease,studybatch,apoe) %>% unique

proportions <- t(apply(norm_counts[,1:8], 1, function(x) (x/sum(x)))) %>% as.data.frame()
value_to_add_to_offset_zeros <- min(proportions[proportions > 0])
proportions <- proportions + value_to_add_to_offset_zeros

CLR_transf_proportions <- t(apply(proportions, 1, function(x) log((x/exp(mean(log(x)))) + 1))) %>% as.data.frame()

CLR_transf_proportions$donorIdUnique <- rownames(CLR_transf_proportions)

CLR_transf_proportions_full <- left_join(CLR_transf_proportions, coldata_short %>% dplyr::select(donorIdUnique, diagnosis_harmonized_by_disease, studybatch, apoe) %>% unique, by = "donorIdUnique")

df_long <- CLR_transf_proportions_full %>% 
  pivot_longer(c(c(paste0("Cluster ",1:(dim(norm_counts)[2]-1)))))

level_order <- c("AD_Control","AD","MS_Control","MS","PD_Control","PD")
df_long$diagnosis_harmonized_by_disease <- factor(df_long$diagnosis_harmonized_by_disease,levels = level_order) 

df_long$name <- factor(df_long$name,levels =c(paste0("Cluster ",1:(dim(norm_counts)[2]-1))))

sig.clusters <- list()

## slipping cluster 4 because it's unidentified:


for(cluster in 1:length(levels(factor(coldata$clusters_res_0.5)))){
  clusterName <- paste0("Cluster ",levels(factor(coldata$clusters_res_0.5))[cluster])
  dataCluster <- df_long[which(df_long$name==clusterName),]
  res <- kruskal.test(value ~ diagnosis_harmonized_by_disease, data = dataCluster)
  # check if there is any significance
  if(res$p.value < 0.05/length(levels(factor(coldata$clusters_res_0.5)))){ # bonferroni correction
    sig.clusters[[clusterName]] <- res
    
    # do a wilcox test to find which specific groups are different within cluster
    res.AD <- wilcox.test(dataCluster$value[which(dataCluster$diagnosis_harmonized_by_disease=="AD_Control")],
                          dataCluster$value[which(dataCluster$diagnosis_harmonized_by_disease=="AD")])
    res.MS <- wilcox.test(dataCluster$value[which(dataCluster$diagnosis_harmonized_by_disease=="MS_Control")],
                          dataCluster$value[which(dataCluster$diagnosis_harmonized_by_disease=="MS")])
    res.PD <- wilcox.test(dataCluster$value[which(dataCluster$diagnosis_harmonized_by_disease=="PD_Control")],
                          dataCluster$value[which(dataCluster$diagnosis_harmonized_by_disease=="PD")])
    
    sig.clusters[["AD"]][[clusterName]] <- res.AD
    sig.clusters[["MS"]][[clusterName]] <- res.MS    
    sig.clusters[["PD"]][[clusterName]] <- res.PD
    
  }
  
}

write.csv(df_long,"~/Documents/AstrocytePaper/2024-11-26 lucast3/Fig6D_AbundancePlotDataforStatistics_Human.csv")
#With CLR, the following are significant:
#AD - clusters 3, 5 and 8
#MS - clusters 1, 5, 6 and 8





