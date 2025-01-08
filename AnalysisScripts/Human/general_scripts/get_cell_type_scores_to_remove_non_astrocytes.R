library(SingleCellExperiment)
library(scater)
library(scran)
library(scran.chan)
library(metaGP)
library(dplyr)

print("Reading combined file")
se <- readRDS("/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.with_logcounts.v3.liddelow_cellranger_rest_cellbender.rds")

print("Loading gene sets")


gac <- installedGeneAnnotationCollection("Neuroinflammation")
mouseHipp <- gac$`Mouse Hippocampus scRNASeq Markers`
mouseHipp4 <- mapAnnotation(mouseHipp, "GMTY17:GRCh38/GRCh38.IGIS4.0.genes@REVISION-3")
ga.ctm <- geneSetList(mouseHipp4) # "cell type markers"
ga.ctm <- lapply(ga.ctm, unlist) %>% lapply(unname)
ga.ctm <- ga.ctm[-length(ga.ctm)]

gaM <- gac$`Mouse Myeloid Activation Coarse Clusters`
gaM4 <- mapAnnotation(gaM, "GMTY17:GRCh38/GRCh38.IGIS4.0.genes@REVISION-3")
gslM4 <- geneSetList(gaM4)$geneSet
names(gslM4)[6] <- "Neutrophil_Monocyte"

gene_sets_gatcm <- names(ga.ctm)
gene_sets_gslm4 <- names(gslM4)

all_scores <- list()

print("Computing gene_sets_gatcm scores")

for (gene_set in gene_sets_gatcm) {
  score_gene_set <- ga.ctm[[gene_set]]
  score_gene_set <- intersect(score_gene_set,rownames(se))
  scores_by_cell <- colMeans(assay(se, 'logcounts')[score_gene_set,],na.rm = TRUE)
  all_scores[[gene_set]] <- scores_by_cell
  }

print("Computing gene_sets_gslm4 scores ")
for (gene_set in gene_sets_gslm4) {
  score_gene_set <- gslM4[[gene_set]]
  score_gene_set <- intersect(score_gene_set,rownames(se))
  scores_by_cell <- colMeans(assay(se, 'logcounts')[score_gene_set,],na.rm = TRUE)
  all_scores[[gene_set]] <- scores_by_cell
  }

saveRDS(all_scores, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/Cell_types_scores_for_space_cleaning.rds")
