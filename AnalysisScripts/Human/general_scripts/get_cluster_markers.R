library(ShadowArray)
library(SingleCellExperiment)
library(scran)

print("Loading in the se file")

se_clean_full <- readRDS("AD_MS_PD_astrocytes.clean.withlogcounts.062822.rds")
se_clean_full$studybatch_gerrits_merged <- ifelse(se_clean_full$studybatch %in% c("gerrits_otc","gerrits_oc"), "gerrits", se_clean_full$studybatch)

louvains <- readRDS("All_louvains_from_scvi.070522.rds")


for (n in 1:ncol(louvains)) {

print(paste("Getting markers for ", colnames(louvains)[n], sep=""))

m.out <- scoreMarkers(se_clean_full, louvains[,n], block=se_clean_full$studybatch_gerrits_merged,assay.type="logcounts",lfc=0)
saveRDS(m.out, paste("m.out.", colnames(louvains)[n], ".rds", sep=""))
}


print("Done!")

