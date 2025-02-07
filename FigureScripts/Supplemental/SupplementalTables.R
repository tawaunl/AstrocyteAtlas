library(openxlsx)
source("~/Documents/scHelpers.R")



# Table1 . Mouse DE by Disease and Pathways -----------------------------------------------
dir <- "~/Documents/AstrocytePaper"
sav.dir <- "~/Documents/AstrocytePaper/Supplemental"
DiseaseDE <- readRDS(file.path(dir,"Figure1","DE_Disease"))
DiseasePaths <- readRDS(file.path(dir,"Figure1","DiseasePathways.rds"))
dsBook <- c(DiseaseDE,DiseasePaths)
openxlsx::write.xlsx(dsBook, file =file.path(sav.dir,"SuppTable1_DiseaseDEandPathways_Mouse.xlsx"),
                     asTable = TRUE, tableStyle = "TableStyleLight9")

# Table 2 -------
# Mouse Cluster Markers ---------------------------------------------
markers <- readRDS(file.path(dir,"FinalIntegratedMarkers.rds"))

wb <- openxlsx::createWorkbook()

for (cluster in 1:length(markers$statistics)) {
  cMarkers <- markers$statistics[[cluster]]
  cMarkers <- cMarkers[order(cMarkers$cohen.mean,decreasing = T),]
  openxlsx::addWorksheet(wb,sheetName =  paste0("Cluster_",cluster))
  openxlsx::writeDataTable(wb, sheet = paste0("Cluster_",cluster), x = cMarkers, colNames = TRUE, rowNames = TRUE, tableStyle = "TableStyleLight9")
}

# Mouse Cluster Pathways --------------------------------------------
files <- list.files(file.path(dir,"Figure2"),pattern="*GOresult")     
paths <- lapply(files, function(x){read.csv(file.path(dir,"Figure2",x))})
names(paths) <- paste0("Cluster_",1:4,"_GSEAGO")

for(path in 1:length(paths)) {
  openxlsx::addWorksheet(wb,sheetName =  names(paths)[path])
  openxlsx::writeDataTable(wb, sheet = names(paths)[path], x = paths[[path]], colNames = TRUE, rowNames = FALSE, tableStyle = "TableStyleLight9")
}

## Save markers and Pathways -----------------------------------------------
openxlsx::saveWorkbook(wb, file.path(sav.dir,"SuppTable2_ClusterMarkersandPathways_Mouse.xlsx"),overwrite = TRUE)


# Human ----------------------------

## markers --------
dir <- "~/Documents/AstrocytePaper/Human"
markers <- readRDS("~/Documents/AstrocytePaper/Human/m.out.harmony_2_donor_groupvar.0.5.rds")
subtypeList <- list("1" = "NEAT1-hi",
                    "2" = "GFAP-hi",
                    "3" = "DST-hi",
                    "4" = "BCYRN-hi",
                    "5" = "ADGRV1-hi",
                    "6" = "NRXN1-hi",
                    "7" = "APOE-hi",
                    "8" = "SLC1A2-hi")
wb <- openxlsx::createWorkbook()

for (cluster in 1:length(markers$statistics)) {
  cMarkers <- markers$statistics[[cluster]]
  cMarkers <- cMarkers[order(cMarkers$cohen.mean,decreasing = T),]
  openxlsx::addWorksheet(wb,sheetName =  paste0("Cluster_",cluster,subtypeList[[cluster]]))
  openxlsx::writeDataTable(wb, sheet = paste0("Cluster_",cluster,subtypeList[[cluster]]), x = cMarkers, colNames = TRUE, rowNames = TRUE, tableStyle = "TableStyleLight9")
}


human.AD <- readRDS(file.path(dir,"res_dsl.AD.cain_fixed.072924.rds"))

human.MS <- readRDS(file.path(dir,"res_dsl.MS.080124.rds"))
human.PD <- readRDS(file.path(dir,"res_dsl.PD.072924.rds"))
# Process Data -------------------------
feat <- readRDS(file.path(dir,"humanFeatures.rds"))

human.AD <- generatePlotTable(human.AD)
human.MS <- generatePlotTable(human.MS)
human.PD <- generatePlotTable(human.PD)

human.AD$diffexpressed <- ifelse(human.AD$FDR<=0.05 & human.AD$dl_mu>=0.5 & human.AD$n_up>=3, "up",
                                 ifelse(human.AD$FDR<=0.05 & human.AD$dl_mu<= -0.5 & human.AD$n_down>=3, "down","unchanged"))
human.MS$diffexpressed <- ifelse(human.MS$FDR<=0.05 & human.MS$dl_mu>=0.5 & human.MS$n_up>2, "up",
                                 ifelse(human.MS$FDR<=0.05 & human.MS$dl_mu<= -0.5 & human.MS$n_down>=3, "down","unchanged"))
human.PD$diffexpressed <- ifelse(human.PD$FDR<=0.05 & human.PD$dl_mu>=0.5 & human.PD$n_up>=2, "up",
                                 ifelse(human.PD$FDR<=0.05 & human.PD$dl_mu<= -0.5 & human.PD$n_down>=2, "down","unchanged"))
protein_coding <- feat %>% dplyr::filter(type == "protein_coding")
human.AD <- human.AD %>% dplyr::filter(symbol %in% protein_coding$symbol)
human.MS <- human.MS %>% dplyr::filter(symbol %in% protein_coding$symbol)
human.PD <- human.PD %>% dplyr::filter(symbol %in% protein_coding$symbol)

openxlsx::addWorksheet(wb,sheetName =  "AD_DE_Human")
openxlsx::writeDataTable(wb, sheet = "AD_DE_Human", x = human.AD, colNames = TRUE, rowNames = TRUE, tableStyle = "TableStyleLight9")

openxlsx::addWorksheet(wb,sheetName =  "MS_DE_Human")
openxlsx::writeDataTable(wb, sheet = "MS_DE_Human", x = human.AD, colNames = TRUE, rowNames = TRUE, tableStyle = "TableStyleLight9")

openxlsx::addWorksheet(wb,sheetName =  "PD_DE_Human")
openxlsx::writeDataTable(wb, sheet = "PD_DE_Human", x = human.AD, colNames = TRUE, rowNames = TRUE, tableStyle = "TableStyleLight9")





paths <- readRDS("~/Documents/AstrocytePaper/Human/HumanDiseasePathwaysGO.rds")

for (disease in 1:length(paths)) {
  res <- paths[[disease]]@result
  res <- res[order(res$NES,decreasing = T),]
  openxlsx::addWorksheet(wb,sheetName =  names(paths)[disease])
  openxlsx::writeDataTable(wb, sheet = names(paths)[disease], x = res, colNames = TRUE, rowNames = TRUE, tableStyle = "TableStyleLight9")
}

openxlsx::saveWorkbook(wb, file.path(sav.dir,"SuppTable3_HumanClusterMarkersandDiseaseDEandPathways.xlsx"),overwrite = TRUE)


# Human GFAP subcluster -------
c2Markers <- readRDS("~/Documents/AstrocytePaper/Human/Cluster_2_astrocytes.m.out.harmony_2_donor_groupvar.0.3.rds")

subtypeList <- list("1" = "HSP90AA1-hi",
                    "2" = "NEAT1-hi",
                    "3" = "CTNND2-hi",
                    "4" = "DPP10-hi")
wb <- openxlsx::createWorkbook()

for (cluster in 1:length(c2Markers$statistics)) {
  cMarkers <- markers$statistics[[cluster]]
  cMarkers <- cMarkers[order(cMarkers$cohen.mean,decreasing = T),]
  openxlsx::addWorksheet(wb,sheetName =  paste0("SubCluster_",cluster,subtypeList[[cluster]]))
  openxlsx::writeDataTable(wb, sheet = paste0("SubCluster_",cluster,subtypeList[[cluster]]), x = cMarkers, colNames = TRUE, rowNames = TRUE, tableStyle = "TableStyleLight9")
}

c2Paths <- readRDS("~/Documents/AstrocytePaper/Human/C2MarkerPathways.rds")

openxlsx::addWorksheet(wb,sheetName =  "SubClusterPathways")
openxlsx::writeDataTable(wb, sheet = "SubClusterPathways", x = c2Paths@compareClusterResult, colNames = TRUE, rowNames = TRUE, tableStyle = "TableStyleLight9")

openxlsx::saveWorkbook(wb, file.path(sav.dir,"SuppTable4_GFAPhisubclusterMarkersandPathways.xlsx"),overwrite = TRUE)






