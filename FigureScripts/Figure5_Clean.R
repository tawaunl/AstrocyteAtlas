# Figure 5
# ======================== Load Libraries ======================================
source("~/Documents/scHelpers.R")
library(dplyr)
library(ggplot2)
library(scater)
library(ggfun)

load("~/Documents/AstrocytePaper/humanFeatures.Rdata")
# Load Data ---------------------------------------
dir <- "~/Documents/AstrocytePaper/Human"
coldata <- readRDS(file.path(dir,
                             "AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.COLDATA.correct_diagnosis_harmonized_by_disease.rds"))

dims <- readRDS(file.path(dir,
                          "AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.onlychanHVG_for_PCA.with_harmony.liddelow_cellranger_rest_cellbender.rds"))
clusters <- readRDS(file.path(dir,
                              "Clusters_res_harmony_2_donor_groupvar_0.5.rds"))
# A. Human UMAP clusters -------------------------------------
dims$clusters <- clusters$membership
harmony_coords <- get_harmony_coords(dims,3)

color.codes <- c(adjustcolor('#31a354', alpha.f = 0.9),
                 '#8856a7', adjustcolor('#fa9fb5', 0.9),
                 '#2b8cbe', adjustcolor('#fec44f', alpha.f = 0.9),
                 adjustcolor('#c51b8a', 0.9),
                 adjustcolor('#d95f0e', alpha.f = 0.9),
                 adjustcolor('turquoise1', alpha.f = 0.9))

zone <- levels(factor(harmony_coords$clusters))
png(file.path(dir,"Figure5","A_ClustersUMAP.png"),
    width = 500,height=500)
ggplot(data=harmony_coords , aes(x=UMAP1, y=UMAP2, colour  = clusters)) + 
  geom_point(size=1,alpha = 0.6) +
  scale_colour_manual(values=setNames(color.codes, zone))+ theme_nothing()
dev.off()


#B. Cells Bar Plot ---------------------
cells <- data.frame(table(dims$clusters))
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
coldata <- as.data.frame(colData(dims))
study_cluster_tmp_by_study <- coldata %>%
  group_by(studybatch, clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

cluster1 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$clusters==1),]
cluster2 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$clusters==2),]
cluster3 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$clusters==3),]
cluster4 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$clusters==4),]
cluster5 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$clusters==5),]
cluster6 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$clusters==6),]
cluster7 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$clusters==7),]
cluster8 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$clusters==8),]

cluster1$scaled <- cluster1$freq/sum(cluster1$freq)
cluster2$scaled <- cluster2$freq/sum(cluster2$freq)
cluster3$scaled <- cluster3$freq/sum(cluster3$freq)
cluster4$scaled <- cluster4$freq/sum(cluster4$freq)
cluster5$scaled <- cluster5$freq/sum(cluster5$freq)
cluster6$scaled <- cluster6$freq/sum(cluster6$freq)
cluster7$scaled <- cluster7$freq/sum(cluster7$freq)
cluster8$scaled <- cluster8$freq/sum(cluster8$freq)
study_cluster_tmp_by_study_scaled <- rbind(cluster1,cluster2,cluster3,cluster4,
                                           cluster5,cluster6,cluster7,cluster8)

#### Plot C -----------
png(file.path(dir,"Figure5","C_StudyBarPlot.png"),
    width = 550,height=500)
ggplot(study_cluster_tmp_by_study_scaled, aes(fill=studybatch, y=scaled, x=clusters)) +
  geom_bar(position="stack", stat="identity", alpha=1) + theme_void()+ theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14)) +
  xlab("Clusters") + ylab("Proportion of cells within a study") +
  theme(axis.title.y =  element_text(angle = 90))
dev.off()

#D. Markers Dotplot -------------------------
markers <- readRDS(file.path(dir,"m.out.harmony_2_donor_groupvar.0.5.rds"))


markers <- markers$statistics
  
m.out.with_symbols <- list()

for (i in 1:length(markers)){
  tmp <- markers[[i]]
  tmp$ID <- rownames(tmp)
  m.out.with_symbols[[i]] <- left_join(tmp, feat, by = "ID")
}

markers_for_plotting_10_each_symbol <- list()
markers_for_plotting_10_each_ID <- list()

for (i in 1:length(m.out.with_symbols)) {
  temp <- as.data.frame(m.out.with_symbols[[i]]) %>% dplyr::filter(cohen.mean!="Inf") %>% arrange(-cohen.mean)
  markers_for_plotting_10_each_symbol[[i]] <- temp$symbol %>% head(10) 
  markers_for_plotting_10_each_ID[[i]] <- temp$ID %>% head(10) 
}

subtypeList <- list("1" = "NEAT1-hi",
                    "2" = "GFAP-hi",
                    "3" = "DST-hi",
                    "4" = "BCYRN-hi",
                    "5" = "ADGRV1-hi",
                    "6" = "NRXN1-hi",
                    "7" = "APOE-hi",
                    "8" = "SLC1A2-hi")
features <- unlist(markers_for_plotting_10_each_symbol)%>% 
  purrr::set_names(rep(subtypeList, each = 10)) 

markers_for_plotting_10_each_ID <- unique(unlist(markers_for_plotting_10_each_ID))
markers_for_plotting_10_each_symbol <- unique(unlist(markers_for_plotting_10_each_symbol))

out_for_plotting <- out[markers_for_plotting_10_each_ID,]
out_for_plotting$clusters_0.5 <- clusters$membership
colData(out_for_plotting) <- data.frame(colData(out_for_plotting)) %>% 
  dplyr::mutate(Cluster = dplyr::recode_factor(clusters_0.5,
                                               !!!subtypeList)) %>% DataFrame()

rowData(out_for_plotting)$Marker <- features[match(rownames(out_for_plotting), features)] |>
  names() %>% 
  factor(levels = unlist(subtypeList))

colors <- c(adjustcolor('#31a354', alpha.f = 0.9), '#8856a7',
            adjustcolor('#fa9fb5', 0.9),  '#2b8cbe',
            adjustcolor('#fec44f', alpha.f = 0.9),
            adjustcolor('#c51b8a', 0.9),
            adjustcolor('#d95f0e', alpha.f = 0.9),
            adjustcolor('turquoise1', alpha.f = 0.9))
rownames(out_for_plotting) <- markers_for_plotting_10_each_symbol
png(file.path(dir,"Figure5","D_MarkersDotPlot.png"),
    width = 2000, height = 2100, res = 300)
plotDots(out_for_plotting, factor(rownames(out_for_plotting),markers_for_plotting_10_each_symbol),
         group="Cluster", scale=TRUE, center=TRUE) + xlab("Cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(colour="black", size = 5.5), axis.ticks.x  = element_blank()) + ylab("Cluster markers")
dev.off()

pdf(file.path(dir,"Figure5","D_UpdatedMarkersDotPlot.png"),width=16,height=10)
out_for_plotting %>% 
  scDotPlot::scDotPlot(features = factor(rownames(out_for_plotting),markers_for_plotting_10_each_symbol),
                       group = "Cluster",
                       #block = "Sample",
                       scale = TRUE,
                       cluster = FALSE,
                       groupAnno = "Cluster",
                       featureAnno = "Marker",
                       annoColors = list("Cluster" = colors,
                                         "Marker" = colors),
                       featureLegends = FALSE,
                       annoHeight = 0.025,
                       annoWidth = 0.1,flipPlot=T,
                       dotColors=c("blue","#FFFFBF", "red"))
dev.off()
