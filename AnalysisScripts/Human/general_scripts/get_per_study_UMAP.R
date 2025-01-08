library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(gridExtra)
library(dplyr)
out <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.with_logcounts.v3.liddelow_cellranger_rest_cellbender.rds")


print("Reading harmony output file")

se_ast_for_pca_post_clean <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.onlychanHVG_for_PCA.with_harmony.liddelow_cellranger_rest_cellbender.rds")

reducedDims(out) <- reducedDims(se_ast_for_pca_post_clean)

get_umap_plot <- function(scvi_coords_final,dataset) {
  ### can be used for dataset specific plotting, but also for clusters
  p <- ggplot(data=scvi_coords_final, aes(x=UMAP1, y=UMAP2)) +
    geom_point(color="grey", size=2,alpha = 1) +
    geom_point(data = scvi_coords_final %>% dplyr::filter((!!sym(dataset))==1), color="dodgerblue2", size=0.5, alpha=0.5)+
    theme_classic() + ggtitle(dataset) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size = 24, face = "bold"))
  return(p)
}


get_scvi_coords <- function(se) {
  scvi_coords <- reducedDims(se)[["UMAP_theta_1.5"]]
  colnames(scvi_coords) <- c("UMAP1","UMAP2")
  scvi_coords <- as.data.frame(cbind(scvi_coords,colData(se)))
  
  results <- fastDummies::dummy_cols(se$studybatch)
  colnames(results) <- unlist(lapply(colnames(results),function(x) strsplit(x,".data_")[[1]][2]))
  results <- results[,-1]
  scvi_coords <- cbind(scvi_coords, results)
  
  return(scvi_coords)
  
}

scvi_coords <- get_scvi_coords(out)


list_of_plots <- list()

studies <- unique(out$studybatch)

for (i in 1:length(studies)) {
  p <- get_umap_plot(scvi_coords,studies[i])
  list_of_plots[[i]] <- p
}


png("/gne/web/dev/apache/htdocs/people/novikovg/Astrocytes_meta/Harmony_AD_MS_PD_integration_2024/AD_PD_MS_studies.UMAPs.png",width = 2500, height = 1500)
do.call("grid.arrange", c(list_of_plots, ncol=6))
dev.off()
