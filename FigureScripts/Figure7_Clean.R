# Clean Version Figure 7

#Load Libraries ----------------------------------
source("~/Documents/scHelpers.R")
library(wesanderson)

# Load data ---------------------------------------
dir <- "~/Documents/AstrocytePaper/Human"
human.AD <- readRDS(file.path(dir,"res_dsl.AD.cain_fixed.072924.rds"))

dims <- readRDS(file.path(dir,
                          "AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.onlychanHVG_for_PCA.with_harmony.liddelow_cellranger_rest_cellbender.rds"))
AD_up_score_for_each_cell <- readRDS("~/Documents/AstrocytePaper/Human/AD_up_score_for_each_cell.rds")
dims$AD_up <- AD_up_score_for_each_cell
coldata <- readRDS(file.path(dir,
                             "AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.COLDATA.correct_diagnosis_harmonized_by_disease.Cain_donors_subsetted.rds"))

# A. AD up Score -------------------------

##Plot A --------------------
png(file.path(dir,"Figure7","A_UMAP_ADScore.png"))
plotReducedDim(dims, dimred = "UMAP_theta_2_donor_groupvar", colour_by = "AD_up",
               point_size = 0.5, point_alpha = 0.5) +
  scale_colour_gradientn(colours = wes_palette("Zissou1", n = 5),
                         name = "AD \nup \nscore \nlogcounts") +
  theme_classic() + theme(legend.title = element_text(size = 8))
dev.off()


# B Violin plot ----------------------------
df_mic <- data.frame(AD_up=AD_up_score_for_each_cell, cluster=dims$clusters)

colors <- c(adjustcolor('#31a354', alpha.f = 0.9), '#8856a7',
            adjustcolor('#fa9fb5', 0.9),  '#2b8cbe',
            adjustcolor('#fec44f', alpha.f = 0.9),
            adjustcolor('#c51b8a', 0.9),  adjustcolor('#d95f0e', alpha.f = 0.9),
            adjustcolor('turquoise1', alpha.f = 0.9))

p <- ggplot(df_mic, aes(cluster,log(AD_up), fill=cluster)) +
  scale_fill_manual(values = colors) +
  geom_violin() +
  theme(legend.position = "none") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  ylab("logCPM") 

#C Cluster 2 Subcluster Sub cluster UMAP --------------------

# D  Cluster 2 SubCluster AD Scoring -----------
c2 <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/another_run_after_removing_high_mito_clusters/subclustering/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.Cluster_2_only.rds")
harmony_gfap <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/another_run_after_removing_high_mito_clusters/subclustering/Cluster2_results/Cluster_2_astrocytes.onlychanHVG_for_PCA.with_harmony.rds")
reducedDims(c2) <- reducedDims(harmony_gfap)

clusters <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/another_run_after_removing_high_mito_clusters/subclustering/Cluster2_results/Cluster_2_astrocytes.Clusters_res_harmony_2_donor_groupvar_0.3.rds")

c2$subcluster_res_0.3 <- clusters$membership


library(batchelor)
library(BiocParallel)
out_gfap <- multiBatchNorm(c2, normalize.all = TRUE,
                           batch = c2$studybatch,
                           BPPARAM=MulticoreParam(workers = 3))

out_gfap$AD_up <- colMeans(assay(out_gfap, 'logcounts')[up_candidates$ID,],na.rm = TRUE)
library(wesanderson)

# Plot D --------------------------
pdf("~/HumanCluster2ADScoring.pdf")
plotReducedDim(out_gfap, dimred = "UMAP_theta_2_donor_groupvar", 
               colour_by = "AD_up", point_size = 0.5, point_alpha = 0.5) +
  scale_colour_gradientn(colours = wes_palette("Zissou1", n = 5), name = "AD \nup \nscore \nlogcounts") + theme_classic() + theme(legend.title = element_text(size = 8))
dev.off()


# E. Cluster 2 Subcluster Abundance -----------------------------
coldata <- as.data.frame(colData(se_gfap))

agg_cl <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/another_run_after_removing_high_mito_clusters/DEGs/aggregated_cells.AD_MS_PD.070724.DonoruniqueID_and_cluster.SE.rds") ## cain donors are fixed here

cain_donors <- as.data.frame(colData(agg_cl)) %>% dplyr::filter(studybatch == "cain") %>% dplyr::select(donor) %>% unique

coldata <- coldata %>% dplyr::filter(! donor %in% cain_donors$donor)


coldata$diagnosis_harmonized_by_disease <-
  ifelse(coldata$studybatch %in% c("jakel", "schirmer", "absinta", "bryois") & coldata$diagnosis_harmonized == "Control", "MS_Control",
         
         ifelse(coldata$studybatch %in% c("smajic","wang", "kamath") &
                  coldata$diagnosis_harmonized == "Control", "PD_Control",
                
                ifelse(coldata$studybatch %in% c("cain", "morabito","gerrits", "smith", "SEA-AD", "liddelow") &
                         coldata$diagnosis_harmonized == "Control", "AD_Control",
                       coldata$diagnosis_harmonized)))
coldata$diagnosis_harmonized_by_disease <- ifelse(coldata$diagnosis_harmonized == "RRMS", "MS", coldata$diagnosis_harmonized_by_disease)

abundances <- table(coldata$subcluster_res_0.3,coldata$donorIdUnique)
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



for (clust in 1:length(levels(factor(coldata$subcluster_res_0.3)))) {
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

level_order <- c("AD_Control","AD")
df_long <- df_long %>% dplyr::filter(diagnosis_harmonized_by_disease %in% c("AD_Control","AD"))

df_long$diagnosis_harmonized_by_disease <- factor(df_long$diagnosis_harmonized_by_disease,levels = level_order)

df_long$name <- factor(df_long$name,levels =c(paste0("Cluster ",1:(dim(norm_counts)[2]-1))))
fill.codes <- c("white","darkgreen")
color.codes <- c("darkgreen","darkgreen")

zone <- levels(factor(df_long$diagnosis_harmonized_by_disease))

### Plot E -----------------
png(file.path(dir,"Figure7","Cluster2Subcluster_Human_CellularityPlot_with_percentages.png"),
    width = 2700,height=800,res=300)
ggplot(df_long, aes(x=diagnosis_harmonized_by_disease, y=value,fill=diagnosis_harmonized_by_disease,colour = diagnosis_harmonized_by_disease)) +
  theme_classic() +
  geom_jitter(color="darkgrey",width = 0.2,size=1,alpha=.5)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5)  +
  facet_wrap(~factor(name),scales = "free", ncol=8) + xlab("Disease Label") +
  theme(axis.text.x = element_text(angle=75,vjust = 0.5),axis.text.y = element_text(size=20)) +
  ylab("Cellularity Percent") +
  theme(strip.text.x = element_text(size = 12,face = "italic"),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14)) +
  scale_y_continuous(limits = c(0, 120), breaks = c(0, 50, 100))+
  scale_fill_manual(values=setNames(fill.codes, zone))+
  scale_color_manual(values=setNames(color.codes, zone))
dev.off()











