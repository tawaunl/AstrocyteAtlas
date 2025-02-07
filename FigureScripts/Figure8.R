# Figure 8 
library(ggplot2)
library(gg4way)
source("~/Documents/scHelpers.R")

mouse.AD <- readRDS("~/Documents/AstrocytePaper/res_dsl_AD_full.060723.rds")
human.AD <- readRDS(paste0("~/Documents/AstrocytePaper/Human/",
                           "res_dsl.AD.cain_fixed.072924.rds"))
human.AD <- generatePlotTable(human.AD)
mouse.AD <- generatePlotTable(mouse.AD)

mouse.AD$symbol <- toupper(mouse.AD$symbol)
mouse.AD$ID <- mouse.AD$symbol
human.AD$ID <- human.AD$symbol
mouse.AD <- as.data.frame(mouse.AD) %>% dplyr::mutate(diffexpressed=ifelse(dl_mu>=0.5 & FDR<=0.05 & n_up>= max(n_tested)/2, "up", 
                                                                   ifelse(dl_mu<=-0.5 & FDR<=0.05 & n_down>= max(n_tested)/2, "down", "unchanged"))) 

human.AD <- as.data.frame(human.AD) %>% dplyr::mutate(diffexpressed=ifelse(dl_mu>=0.5 & FDR<=0.05 & n_up>= max(n_tested)/2, "up", 
                                                                   ifelse(dl_mu<=-0.5 & FDR<=0.05 & n_down>= max(n_tested)/2, "down", "unchanged")))                                    

AD <- mouse.AD$symbol[which(mouse.AD$diffexpressed=="up")]


AD_list <- list("mAD vs Control"=mouse.AD,"hAD vs Control"=human.AD)


library(gg4way)
# Process AD ------------------
AD_plot <- gg4way(DGEdata = AD_list,
                  y = "mAD vs Control",FDR = "PValue",
                  x = "hAD vs Control", sep = " vs ",
                  logFC = "dl_mu",logFCcutoff = 0.5, label=T, textSize =10,
                  colorVector = c("grey80", "goldenrod1", "red", "mediumblue")) +
  ggplot2::theme(legend.title = ggplot2::element_text(size = 14),
                 legend.text = ggplot2::element_text(size = 14),axis.title = element_text(size=14))


ff <- AD_plot |> getTotals()

data <- AD_plot$data
shared_up <- data %>% filter(`mAD vs Control Direction` =="Up" &
                               `hAD vs Control Direction`=="Up" &
                               Significant=="Significant in Both")
human_up <- data %>% filter(`hAD vs Control Direction`=="Up" &
                              Significant == "Significant in hAD vs Control")
mouse_up <- data %>% filter(`mAD vs Control Direction`=="Up" &
                              Significant == "Significant in mAD vs Control")


## Pathway analysis ---------------------
library(clusterProfiler)
library(org.Hs.eg.db)
pathData <- list("AD Shared"=shared_up,"AD Human Distinct"=human_up,"AD Mouse Distinct"=mouse_up)

enrichPaths <- lapply(pathData,function(x){
  enrichGO(x$symbol,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont="BP",readable = T)
})

enrichPlots <- list()
for(set in 1:length(enrichPaths)){
  name <- names(enrichPaths)[set]
  enrichPlots[[name]] <- dotplot(enrichPaths[[set]]) +
    ggtitle(paste("GO:BP",name,"Pathways")) +
    theme(axis.title.x = element_text(size = 18,face="bold"),
          axis.text.x = element_text(angle = 90,size=14),
          axis.text.y = element_text(size=14)) +
    theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
}

# A. AD 4 way ----------------------------

cairo_pdf("~/Documents/AstrocytePaper/Figure8/MouseADvsHumanAD_4way.pdf",
    width=15,height = 12)
plot_4way( human.AD,mouse.AD,  "Human AD","Mouse AD models")
dev.off()

png(filename ="/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Figures/MouseADvsHumanAD_pathways.png",
    width=3000,height=1000,res=150)
gridExtra::grid.arrange(grobs = enrichPlots,ncol=3)
dev.off()
saveRDS(enrichPaths,"/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Rscripts/Mouse_HumanComparison/ADComparisonPathways.rds")

# C. Plot AD Comparison Plots ------------------
enrichPaths<- readRDS("~/Documents/AstrocytePaper/Figure8/ADComparisonPathways.rds")

topPaths <- lapply(enrichPaths,function(x){
  data<- x@result
  data <- data[order(data$p.adjust),]
})
names(topPaths) <- c("AD Shared Genes","AD Human Distinct", "AD Mouse Distinct")
df <- bind_rows(topPaths, .id = "Cluster")

df <- df %>% 
  separate(GeneRatio,c("count","total"),sep="/") %>% 
  mutate(GeneRatio = as.double(count)/as.double(total))
## Function to filter first N GO term under specified p.adjust cutoff.
filterCategory <- function(df, N, cutoff){
  df_topn <- NULL
  for(i in levels(factor(df$Cluster))){
    df_slice <- df %>% filter(Cluster == i) %>% 
      arrange(p.adjust) %>% 
      filter(p.adjust< cutoff) %>% 
      head(N)
    df_topn <- bind_rows(df_topn, df_slice)
  }
  return(df_topn)
}
df_top10 <- filterCategory(df,5,0.05)      # Filter the results to show only top 10 GO terms less than p.adjust 0.05 in all cluster.
descriptions <- df_top10$ID

full <- df %>% filter(ID %in% descriptions)
full1 <- filterCategory(full,10,0.05) # filter top from all paths
library(stringr)
full1$Name <- str_wrap(full1$Description, width = 30) # wrap axis labels
# reorder to plot how i want
full1$Name <- factor(full1$Name,
                     levels = c(full1$Name[2],full1$Name[3],full1$Name[1],full1$Name[8],
                                full1$Name[7],full1$Name[6],full1$Name[4],full1$Name[5],
                                full1$Name[9],full1$Name[11],
                                full1$Name[13],full1$Name[16],full1$Name[15])
)


full1$Cluster <- factor(full1$Cluster,levels = c("AD Shared Genes","AD Mouse Distinct",
                                                 "AD Human Distinct"))
## Plot C -------------------------------
pdf("~/Documents/AstrocytePaper/Figure8/ADPaths.pdf",
    height=12,width=8)
ggplot(full1, aes(x=Cluster, y=Name)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_colour_gradient(limits=c(min(full1$p.adjust), max(full1$p.adjust)), low="red",high="blue") +
  theme_bw() + ggtitle("AD Comparison Pathways") +
  theme(axis.title = element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90,size=14),
        axis.text.y = element_text(size=20)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))+ 
  scale_size(range = c(3,11))
dev.off()


# B. MS 4-way ------------------

mouse.MS <- readRDS("~/Documents/AstrocytePaper/res_dsl_MS_full.060723.rds")
human.MS <- readRDS(paste0("~/Documents/AstrocytePaper/Human/",
                           "res_dsl.MS.080124.rds"))
human.MS <- generatePlotTable(human.MS)
mouse.MS <- generatePlotTable(mouse.MS)

mouse.MS$symbol <- toupper(mouse.MS$symbol)
mouse.MS$ID <- mouse.MS$symbol
human.MS$ID <- human.MS$symbol
MS_list <- list("mMS vs Control"=mouse.MS,"hMS vs Control"=human.MS)

MS_plot <- gg4way(DGEdata = MS_list,
                  x = "mMS vs Control",FDR = "PValue",
                  y = "hMS vs Control", sep = " vs ",
                  logFC = "LogFC",logFCcutoff = 0.5, label=T, textSize =10,
                  colorVector = c("grey80", "goldenrod1", "red", "mediumblue")) +
  ggplot2::theme(legend.title = ggplot2::element_text(size = 14),
                 legend.text = ggplot2::element_text(size = 14),axis.title = element_text(size=14))


ff <- MS_plot |> getTotals()

data <- MS_plot$data
shared_up <- data %>% filter(`mMS vs Control Direction` =="Up" &
                               `hMS vs Control Direction`=="Up" &
                               Significant=="Significant in Both")
human_up <- data %>% filter(`hMS vs Control Direction`=="Up" &
                              Significant == "Significant in hMS vs Control")
mouse_up <- data %>% filter(`mMS vs Control Direction`=="Up" &
                              Significant == "Significant in mMS vs Control")

## Pathway analysis ---------------------

pathData <- list(shared=shared_up,human=human_up,mouse=mouse_up)

enrichPaths <- lapply(pathData,function(x){
  enrichGO(x$symbol,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont="BP",readable = T)
})

enrichPlots <- list()
for(set in 1:length(enrichPaths)){
  name <- names(enrichPaths)[set]
  enrichPlots[[name]] <- dotplot(enrichPaths[[set]]) +
    ggtitle(paste("GO:BP",name,"Pathways")) +
    theme(axis.title.x = element_text(size = 18,face="bold"),
          axis.text.x = element_text(angle = 90,size=14),
          axis.text.y = element_text(size=14)) +
    theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
}

label <- MS_plot$data |>
  dplyr::filter(Significant == "Significant in Both") |>
  dplyr::pull(symbol)
## Plot B ----------
cairo_pdf("~/Documents/AstrocytePaper/Figure8/MouseMSvsHumanMS.pdf",
    width=15,height=12)
plot_4way( human.MS,mouse.MS,  "Human MS","Mouse MS models")
dev.off()

png(filename ="/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Figures/MouseADvsHumanMS_pathways.png",width=3000,height=1000,res=150)
gridExtra::grid.arrange(grobs = enrichPlots,ncol=3)
dev.off()
saveRDS(enrichPaths,"/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Rscripts/Mouse_HumanComparison/MSComparisonPathways.rds")

#D. Compare MS ----------------------
enrichPaths <- readRDS("~/Documents/AstrocytePaper/Figure8/MSComparisonPathways.rds")
topPaths <- lapply(enrichPaths,function(x){
  data<- x@result
  data <- data[order(data$p.adjust),]
})
names(topPaths) <- c("MS Shared Genes","MS Human Distinct", "MS Mouse Distinct")
df <- bind_rows(topPaths, .id = "Cluster")

df <- df %>% 
  separate(GeneRatio,c("count","total"),sep="/") %>% 
  mutate(GeneRatio = as.double(count)/as.double(total))

## Function to filter first N GO term under specified p.adjust cutoff.
filterCategory <- function(df, N, cutoff){
  df_topn <- NULL
  for(i in levels(factor(df$Cluster))){
    df_slice <- df %>% filter(Cluster == i) %>% 
      arrange(p.adjust) %>% 
      filter(p.adjust< cutoff) %>% 
      head(N)
    df_topn <- bind_rows(df_topn, df_slice)
  }
  return(df_topn)
}
df_top10 <- filterCategory(df,5,0.05)      # Filter the results to show only top 10 GO terms less than p.adjust 0.05 in all cluster.
descriptions <- df_top10$ID

full <- df %>% filter(ID %in% descriptions)

## Order the GO term in a style similar to dotplot.

full1 <- filterCategory(full,10,0.05)
full1$Name <- str_wrap(full1$Description, width = 30) # wrap axis labels
full1 <- full1 %>% filter(Cluster != "MS Shared Genes") # not enough shared to genes to be meaningful

full1$Name <- factor(full1$Name,
                     levels = c(full1$Name[5],full1$Name[4],full1$Name[2],full1$Name[3],
                                full1$Name[1],full1$Name[9],full1$Name[10],full1$Name[6],
                                full1$Name[7],full1$Name[8]))


full1$Cluster <- factor(full1$Cluster,levels = c("MS Mouse Distinct",
                                                 "MS Human Distinct"))
## Plot D -------------------
pdf("~/Documents/AstrocytePaper/Figure8/MSPaths.pdf",
    width=8,height=12)
ggplot(full1, aes(x=Cluster, y=Name)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_colour_gradient(limits=c(min(full1$p.adjust), max(full1$p.adjust)), low="red",high="blue") +
  theme_bw() + ggtitle("MS Comparison Pathways") +
  theme(axis.title = element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90,size=14),
        axis.text.y = element_text(size=20)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))+ 
  scale_size(range = c(3,11))
dev.off()


# E. MetaNeighbor ----------------------
library(SingleCellExperiment)
library(Seurat)
library(MetaNeighbor)

# ========================== *** Load Data *** =================================
mouse.data <- readRDS("/gstore/project/neurodegen_meta/data/AstrocyteIntegration_AmbientRemoved_filtered_noneuron.RDS")
human.data <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/another_run_after_removing_high_mito_clusters/integration/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.WITH_logcounts.liddelow_cellranger_rest_cellbender.rds")


mouse.data$Species <- "Mouse"
mouse.counts <- GetAssayData(mouse.data,slot="counts",assay = "RNA")
rownames(mouse.counts) <- mouse.data@assays$RNA@meta.features$ID
rownames(mouse.counts) <- toupper(rownames(mouse.counts)) # loose conversion to human

## Prepare Human Data ==========================================================
human.data$Species <- "Human"
human.data$Disease <- factor(human.data$diagnosis_harmonized)

human.counts <- counts(human.data)
rownames(human.counts) <- rowData(human.data)$symbol

### Subset counts tables ====================
# find overlapping genes
x <- rownames(human.counts) %in% rownames(mouse.counts)
f.human.counts <- human.counts[x,]
x <- rownames(f.mouse.counts) %in% rownames(f.human.counts)
f.mouse.counts <- f.mouse.counts[x,]
f.mouse.counts <- f.mouse.counts[!duplicated(rownames(f.mouse.counts)),]
f.human.counts <- as(f.human.counts,"dgCMatrix")


f.human.counts <- f.human.counts[order(rownames(f.human.counts)),]
f.mouse.counts <- f.mouse.counts[order(rownames(f.mouse.counts)),]
combined.counts <- cbind(f.human.counts,f.mouse.counts)

x <- factor(mouse.data$finalClusters)
levels(x) = paste0("Mouse_Cluster_",levels(mouse.data$finalClusters))
mouse.data$combined.clusters <- x
x <- factor(human.data$clusters_res_0.5)
levels(x) = paste0("Human_Cluster_",levels(human.data$clusters_res_0.5))
human.data$combined.clusters <- x


meta <- data.frame(Species=c(human.data$Species,mouse.data$Species),
                   Disease=c(as.character(human.data$Disease),as.character(mouse.data$Disease)),
                   Clusters=c(human.data$combined.clusters,mouse.data$combined.clusters),
                   Sample=c(human.data$donor_for_DE,mouse.data$Sample),
                   StudyName=c(as.character(human.data$studybatch),mouse.data$StudyName))


combined <- SingleCellExperiment(assays = list(counts = combined.counts),colData=meta)
saveRDS(combined,"/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Rscripts/Mouse_HumanComparison/combinedSE.rds")

var_genes = variableGenes(dat = combined, exp_labels =combined$Species)

celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = combined,
                             study_id = combined$Species,
                             cell_type = combined$Clusters,
                             fast_version = TRUE)
saveRDS(celltype_NV, "/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Rscripts/Mouse_HumanComparison/MetaNeighbor.rds")


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

my_sample_col <- data.frame(Species = rep(c("Human", "Mouse"), c(8,4)))
row.names(my_sample_col) <- colnames(celltype_NV)
## Plot E -------------------------------
pheatmap::pheatmap(celltype_NV, color=cols,breaks=breaks,
                   annotation_col = my_sample_col,treeheight_row = 20,
                   treeheight_col = 20,annotation_row = my_sample_col,
                   filename="/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Rscripts/Mouse_HumanComparison/MetaNeighborHeatmap_new.png",
                   width=10,height=8
)

