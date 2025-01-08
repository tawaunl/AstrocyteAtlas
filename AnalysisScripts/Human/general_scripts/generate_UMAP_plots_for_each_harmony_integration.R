library(SingleCellExperiment)
library(scater)
print("Reading combined file")
out <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.with_logcounts.v3.liddelow_cellranger_rest_cellbender.rds")


print("Reading harmony output file")

se_ast_for_pca_post_clean <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.onlychanHVG_for_PCA.with_harmony.liddelow_cellranger_rest_cellbender.rds")

reducedDims(out) <- reducedDims(se_ast_for_pca_post_clean)

dimnames <- reducedDimNames(se_ast_for_pca_post_clean)[grepl("UMAP",reducedDimNames(se_ast_for_pca_post_clean))]

print("plotting UMAPs")
for (dim in dimnames) {

png(paste("/gne/web/dev/apache/htdocs/people/novikovg/Astrocytes_meta/Harmony_AD_MS_PD_integration_2024/", dim, ".png", sep=""),width = 2000, height = 1500)
print(plotReducedDim(out, dimred = dim))
dev.off()
}

