library(Seurat)
library(ggplot2)
source("~/Documents/scHelpers.R")
library(dplyr)
library(scales)

#A. C2 Marker Dot ----------
## Load Libraries -------
library(scater)
library(SingleCellExperiment)
library(scran.chan)
get_cluster_marker_symbols <- function(res, feat) {
  
  m.out <- readRDS(paste("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/another_run_after_removing_high_mito_clusters/subclustering/Cluster2_results/Cluster_2_astrocytes.m.out.harmony_2_donor_groupvar.", res, ".rds", sep=""))
  m.out <- m.out$statistics
  
  m.out.with_symbols <- list()
  
  for (i in 1:length(m.out)){
    tmp <- m.out[[i]]
    tmp$ID <- rownames(tmp)
    m.out.with_symbols[[i]] <- left_join(tmp, feat, by = "ID")
  }
  
  return(m.out.with_symbols)  
}

# Load Data ---------
c2 <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/another_run_after_removing_high_mito_clusters/subclustering/AD_MS_PD_astrocytes.CLEAN_from_non_astrocytes_AND_clusters_high_in_mito.cellbender_counts.20cell_filter.Cluster_2_only.rds")
clusters <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/another_run_after_removing_high_mito_clusters/subclustering/Cluster2_results/Cluster_2_astrocytes.Clusters_res_harmony_2_donor_groupvar_0.3.rds")
dims <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/another_run_after_removing_high_mito_clusters/subclustering/Cluster2_results/Cluster_2_astrocytes.onlychanHVG_for_PCA.with_harmony.rds")
reducedDims(c2) <- reducedDims(dims)
logcounts(c2) <- normalizeCounts(c2)

library(dplyr)
library(tidyseurat)
library(tidyr)
feat <- as.data.frame(genomitory::getFeatures("GMTY17:GRCh38/GRCh38.IGIS4.0.genes.rds@REVISION-3")) %>% dplyr::select(ID, symbol)

m.out <- markers$statistics

m.out.with_symbols <- list()

for (i in 1:length(m.out)){
  tmp <- m.out[[i]]
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
subtypeList <- list("1" = "HSP90AA1-hi",
                    "2" = "NEAT1-hi",
                    "3" = "CTNND2-hi",
                    "4" = "DPP10-hi")
features <- unlist(markers_for_plotting_10_each_symbol)%>% 
  purrr::set_names(rep(subtypeList, each = 10)) 

markers_for_plotting_10_each_ID <- unique(unlist(markers_for_plotting_10_each_ID))
markers_for_plotting_10_each_symbol <- unique(unlist(markers_for_plotting_10_each_symbol))

out_for_plotting <- c2[markers_for_plotting_10_each_ID,]
rownames(out_for_plotting) <- rowData(out_for_plotting)$symbol

out_for_plotting$clusters_0.5 <- clusters$membership
colData(out_for_plotting) <- data.frame(colData(out_for_plotting)) %>% 
  dplyr::mutate(Cluster = dplyr::recode_factor(clusters_0.5,
                                               !!!subtypeList)) %>% DataFrame()

rowData(out_for_plotting)$Marker <- features[match(rownames(out_for_plotting), features)] |>
  names() %>% 
  factor(levels = unlist(subtypeList))
colors <- c(adjustcolor('darkgrey', alpha.f = 0.5), adjustcolor('aquamarine3', 0.8),  adjustcolor('chocolate2', alpha.f = 0.8), adjustcolor('cornflowerblue', 0.8))
features <- features[-which(duplicated(features)==TRUE)]
pdf("/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Figures/HumanCluster2_markers.top_10_markers.res_0.3_SCDOT.pdf")
out_for_plotting %>% 
  scDotPlot::scDotPlot(features = features,
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
                       annoWidth = 0.1)
dev.off()
#B. C2 Pathways ----------------------


clusters_gfap_res_0.3_markers.with_symbols <- get_cluster_marker_symbols(0.3, feat)
gcSample <- list()

gcSample <- lapply(clusters_gfap_res_0.3_markers.with_symbols,function(x){
  x <- x[order(x$cohen.mean,decreasing = T),]
  setNames(x$cohen.mean,x$symbol)
  
})

names(gcSample) <- c("HSP90AA1-hi",
                     "NEAT1-hi",
                     "CTNND2-hi",
                     "DPP10-hi")

ck.GO <- compareCluster(geneClusters = gcSample, fun = "gseGO",nPermSimple = 10000,
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont="BP",eps=0,keyType="SYMBOL")
#C. Snap Programs --------
dir <- "~/Documents/AstrocytePaper"
data<- readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))

sav.dir <- "~/Documents/AstrocytePaper/Supplemental/SuppFigure6"
snapA <- readRDS("~/Documents/AstrocytePaper/Snap-aProgram.rds")

counts <- GetAssayData(data,assay = "RNA",
                       layer = "counts")
counts <- round(counts)
counts <- counts[-which(rowMeans(counts)==0),]
se <- SingleCellExperiment::SingleCellExperiment(list(counts=counts))
colData(se) <- DataFrame(data@meta.data)

library(batchelor)
out_gfap <- multiBatchNorm(se, normalize.all = TRUE, batch = se$StudyName)

score_intersect <- intersect(snapA,rownames(data))
scores_by_cell <- colMeans(
  GetAssayData(data,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
data[["SnapAProgram"]] <- scores_by_cell

scvi_coords <- get_scvi_coords(data,data$finalClusters)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))
text <- data$finalClusters
text_x <- vapply(split(scvi_coords$UMAP1, text), median, FUN.VALUE=0)
text_y <- vapply(split(scvi_coords$UMAP2, text), median, FUN.VALUE=0)
## Plot C -----------------------
cairo_pdf(file.path(sav.dir,"SnapAScoringUMAP.pdf"),
          width = 8,height=6)
ggplot(scvi_coords%>%
         arrange(SnapAProgram),aes(x=UMAP1, y=UMAP2, colour=SnapAProgram)) +
  ggrastr::geom_point_rast(size=1.5) +
  theme_classic() +
  theme(plot.title = element_text(size=28,hjust = 0.5),
        legend.position = "right", 
        legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) + ggtitle("SNAP-A Program")+
  scale_colour_gradientn(colours =  RColorBrewer::brewer.pal(5,"Purples"))
dev.off()



