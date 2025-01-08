library(SingleCellExperiment)
library(batchelor)
library(scran.chan)
library(dplyr)
library(scran)
library(scater)


path <- "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/scbasic_runs_for_all_datasets/"

smajic <- readRDS(paste(path,"DS000015518_Smajic.scbasic_out.only_astrocytes.rds", sep=""))
absinta <- readRDS(paste(path,"DS000015618_Absinta.scbasic_out.only_astrocytes.rds", sep=""))
jakel <-  readRDS(paste(path, "DS000015650_Jakel.scbasic_out.only_astrocytes.rds", sep=""))
schirmer <-  readRDS(paste(path, "DS000015651_Schirmer.scbasic_out.only_astrocytes.rds", sep=""))
cain <- readRDS(paste(path, "DS000016300_Cain.scbasic_out.only_astrocytes.rds", sep=""))
morabito <- readRDS(paste(path, "DS000016387_Morabito.scbasic_out.only_astrocytes.rds", sep=""))
liddelow <- readRDS(paste(path, "DS000016463_Liddelow.scbasic_out.only_astrocytes.rds", sep=""))
wang <- readRDS(paste(path, "DS000016644_Wang.scbasic_out.only_astrocytes.rds", sep=""))
gerrits <-  readRDS(paste(path, "DS000016819_Gerrits.scbasic_out.only_astrocytes.rds", sep=""))
smith <- readRDS(paste(path, "DS000016915_Smith.scbasic_out.only_astrocytes.rds", sep=""))


all_inputs <- list(smajic, absinta, jakel, schirmer, cain, morabito, liddelow, wang, gerrits, smith)


common <- lapply(all.inputs, rownames)
common <- Reduce(intersect, common)
print(length(common)) 

all.inputs <- lapply(all.inputs, function(x) x[common,])

### do prevent the error about mismatched rowData
for (i in 1:length(all.inputs)) {
  rowData(all.inputs[[i]]) <- NULL
}

rse <- do.call(combineCols, unname(all.inputs))
sce_filtered <- as(rse, "SingleCellExperiment")

meta_useful <- c("age","diagnosis_harmonized","batch","sex","studybatch","apoe", "donor")
colData(sce_filtered) <- colData(sce_filtered)[,meta_useful]

metadata(sce_filtered) <- list()
assays(sce_filtered) <- list(counts = assays(sce_filtered)[[1]])
reducedDims(sce_filtered) <- list()


donors_to_keep <- names(table(sce_filtered$donor))[as.logical(table(sce_filtered$donor) >= 20)]
coldata_final <- as.data.frame(colData(sce_filtered))
cells_to_keep <- rownames(coldata_final %>% dplyr::filter(donor %in% donors_to_keep))
sce_filtered.more_than_20_cells <- sce_filtered[,cells_to_keep]


out <- multiBatchNorm(sce_filtered.more_than_20_cells, normalize.all = TRUE, batch = sce_filtered.more_than_20_cells$studybatch, BPPARAM = SnowParam(workers = 14))

saveRDS(out, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/AD_MS_PD_astrocytes.cellbender_counts.no_allen_bryois_kamath.20cell_filter.with_logcounts.rds")


