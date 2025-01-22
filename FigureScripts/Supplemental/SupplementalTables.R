# Table1 . Mouse DE by Disease and Pathways -----------------------------------------------
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
