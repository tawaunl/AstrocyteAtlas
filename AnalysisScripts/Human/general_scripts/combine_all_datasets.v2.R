library(SingleCellExperiment)
library(batchelor)
library(scran.chan)
library(dplyr)
library(scran)
library(scater)
library(harmony)
library(BiocParallel)

path <- "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/scbasic_runs_for_all_datasets/"

smajic <- readRDS(paste(path,"DS000015518_Smajic.scbasic_out.only_astrocytes.rds", sep=""))
print("Done loading smajic")

absinta <- readRDS(paste(path,"DS000015618_Absinta.scbasic_out.only_astrocytes.rds", sep=""))
print("Done loading absinta")

jakel <-  readRDS(paste(path, "DS000015650_Jakel.scbasic_out.only_astrocytes.rds", sep=""))
print("Done loading jakel")

schirmer <-  readRDS(paste(path, "DS000015651_Schirmer.scbasic_out.only_astrocytes.rds", sep=""))
print("Done loading schirmer")

cain <- readRDS(paste(path, "DS000016300_Cain.scbasic_out.only_astrocytes.rds", sep=""))
print("Done loading cain")

morabito <- readRDS(paste(path, "DS000016387_Morabito.scbasic_out.only_astrocytes.rds", sep=""))
print("Done loading morabito")

liddelow <- readRDS(paste(path, "DS000016463_Liddelow.scbasic_out.only_astrocytes.rds", sep=""))
print("Done loading liddelow")

wang <- readRDS(paste(path, "DS000016644_Wang.scbasic_out.only_astrocytes.rds", sep=""))
print("Done loading wang")

gerrits <-  readRDS(paste(path, "DS000016819_Gerrits.scbasic_out.only_astrocytes.rds", sep=""))
print("Done loading gerrits")

smith <- readRDS(paste(path, "DS000016915_Smith.scbasic_out.only_astrocytes.rds", sep=""))
print("Done loading smith")

seaad_1  <- readRDS(paste(path, "sea_ad_1_sce_05-01-24.scbasic_out.no_logcounts.only_astrocytes.rds", sep=""))
seaad_2 <- readRDS(paste(path, "sea_ad_2_sce_05-01-24.scbasic_out.no_logcounts.only_astrocytes.rds", sep=""))
seaad_3 <- readRDS(paste(path, "sea_ad_3_sce_05-01-24.scbasic_out.no_logcounts.only_astrocytes.rds", sep=""))
seaad_4 <- readRDS(paste(path, "sea_ad_4_sce_05-01-24.scbasic_out.no_logcounts.only_astrocytes.rds", sep=""))
print("Done loading all seattle datasets")

bryois_1 <- readRDS(paste(path, "Bryois_part1.scbasic_out.only_astrocytes.no_logcounts.rds", sep=""))
bryois_2 <- readRDS(paste(path, "Bryois_part2.scbasic_out.only_astrocytes.no_logcounts.rds", sep=""))
print("Done loading all bryois")

kamath <- readRDS(paste(path, "Kamath.scbasic_out.no_logcounts.only_astrocytes.rds", sep=""))
print("Done loading kamath")

colnames(kamath) <- colData(kamath)$barcode

rowdata_bryois1 <- DataFrame(symbol = rownames(bryois_1))
rowdata_bryois2 <- DataFrame(symbol = rownames(bryois_2))
rowdata_kamath <- DataFrame(symbol = rownames(kamath))


feat <- as.data.frame(genomitory::getFeatures("GMTY17:GRCh38/GRCh38.IGIS4.0.genes.rds@REVISION-3"))

rowdata_bryois1 <- left_join(as.data.frame(rowdata_bryois1), as.data.frame(feat), by = "symbol", multiple="first") %>% dplyr::select(symbol, ID)
rowdata_bryois2 <- left_join(as.data.frame(rowdata_bryois2), as.data.frame(feat), by = "symbol", multiple="first") %>% dplyr::select(symbol, ID)
rowdata_kamath <- left_join(as.data.frame(rowdata_kamath), as.data.frame(feat), by = "symbol", multiple="first") %>% dplyr::select(symbol, ID)

rowData(bryois_1) <- rowdata_bryois1
rowData(bryois_2) <- rowdata_bryois2
rowData(kamath) <- rowdata_kamath


bryois_1 <- bryois_1[!is.na(rowData(bryois_1)$ID),]
bryois_2 <- bryois_2[!is.na(rowData(bryois_2)$ID),]
kamath <- kamath[!is.na(rowData(kamath)$ID),]

rownames(bryois_1) <- rowData(bryois_1)$ID
rownames(bryois_2) <- rowData(bryois_2)$ID
rownames(kamath) <- rowData(kamath)$ID

all.inputs <- list(smajic, absinta, jakel, schirmer, cain, morabito, liddelow, wang, gerrits, smith, kamath, seaad_1, seaad_2, seaad_3, seaad_4, bryois_1, bryois_2)


common <- lapply(all.inputs, rownames)
common <- Reduce(intersect, common)
print(length(common)) 

all.inputs <- lapply(all.inputs, function(x) x[common,])

### do prevent the error about mismatched rowData
#for (i in 1:length(all.inputs)) {
 # rowData(all.inputs[[i]]) <- NULL
#}

rse <- do.call(S4Vectors::combineCols, unname(all.inputs))
sce_filtered <- as(rse, "SingleCellExperiment")

meta_useful <- c("age","diagnosis_harmonized","batch","sex","studybatch","apoe", "donor", "barcode")
colData(sce_filtered) <- colData(sce_filtered)[,meta_useful]

metadata(sce_filtered) <- list()
assays(sce_filtered) <- list(counts = assays(sce_filtered)[[1]])
reducedDims(sce_filtered) <- list()

donors_to_keep <- names(table(sce_filtered$donor))[as.logical(table(sce_filtered$donor) >= 20)]
sce_filtered.more_than_20_cells <- sce_filtered[,sce_filtered$donor %in% donors_to_keep]

print("Starting to run multiBatchNorm")

out <- multiBatchNorm(sce_filtered.more_than_20_cells, normalize.all = TRUE, batch = sce_filtered.more_than_20_cells$studybatch, BPPARAM=MulticoreParam(workers = 14))

print("Saving combine dataset file")
saveRDS(out, "/gstore/data/omni/neuroscience/AD/Astrocytes_AD_MS_PD_meta/AD_MS_PD_cellbender_integration/harmony_integration/integration_with_cellbender_counts/AD_MS_PD_astrocytes.cellbender_counts.20cell_filter.with_logcounts.v2.rds")
print("Done!")

