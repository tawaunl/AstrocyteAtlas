# Label Astrocytes in cell bender counts datasets

# ===================== Load Libraries ====================================

pkgs <- c("Seurat","SingleCellExperiment","ggplot2","R.utils",
          "sctransform","future","harmony","Matrix","patchwork",
          "glmGamPoi","cowplot","scater","gridExtra","RColorBrewer","dplyr",
          "gtable","ggpubr","ggrepel","wesanderson","edgeR","tidyr","scran","scater",
          "gp.sa.solo","SingleR")
for(p in pkgs){
  suppressPackageStartupMessages(library(p, quietly=TRUE,character.only=TRUE))
}
options(future.globals.maxSize= 1e+11)
data.dir <- "/gstore/data/project/neurodegen_adapt_immune/02_tcell_meta_analysis/00_data/00_cellbender/"
main.dir <- "/gstore/data/project/neurodegen_meta/data/cellbender/"


# ==================== Load reference dataset ==================================

datasets <- c("ngs2330","ngs2453","ngs2664","ngs2722","ngs3033" ,"ngs3111","ngs3971")

#load reference dataset
ref <- DataSetDB::getDatasetAsSE("DS000002655") #NGS2722 labels
names(assays(ref)) <- c("counts","logcounts")

ref <- ref[,colSums(counts(ref)) > 0]
ref <- ref[,!is.na(ref$biologicalInterest)]
features <- rownames(ref)
symbols <- ref@rowRanges@elementMetadata@listData[["symbol"]]
ref <- as(ref,"SummarizedExperiment")
rownames(ref) <- symbols


# ==================== Label cells and subset on astrocytes ====================
library(BiocParallel)

for (dataset in datasets) {
  if(file.exists(paste0(data.dir,dataset,"/",dataset,"_merged_cellbender_may30.rds"))){
    data <- readRDS(paste0(data.dir,dataset,"/",dataset,"_merged_cellbender_may30.rds"))
  }
  else{ data <- readRDS(paste0(data.dir,dataset,"/",dataset,"_merged_cellbender_may23.rds")) }
  
  cat(paste0("*****Loading Dataset ",dataset,"********\n\n\n\n"))
  data <- data[,colSums(counts(data)) > 0]
  data.se <- as(data,"SummarizedExperiment")
  rownames(data.se) <- data@rowRanges@elementMetadata@listData[["symbol"]]
  
  pred <- SingleR(test = data.se, ref = ref,
                  labels = ref$biologicalInterest, de.method="wilcox", BPPARAM = MulticoreParam() )
  # anno <- gp.sa.solo::runSingleCellAnnotate(data.se,ref = ref,
  #                               label.field = "biologicalInterest",
  #                               assay.ref = "logcounts",
  #                               commit = "never",
  #                               lognorm.ref = TRUE,
  #                               save.all = FALSE,
  #                               ncpu = 5,
  #                               fname = paste0("~/MetaAnalysis/Annotation.Rmd"))
  
  
  data$FineLabels <- pred$labels
  data$PrunedLabels <- pred$pruned.labels
  dir.create(file.path(main.dir,toupper(dataset)),showWarnings=FALSE)
  saveRDS(data, paste0(main.dir,toupper(dataset),"/",
                       toupper(dataset), "_merged_labeled.rds"))
  astros <- subset(data, ,FineLabels == "astrocyte")
  saveRDS(astros, paste0(main.dir,toupper(dataset),"/",
                         toupper(dataset), "_Astrocytes_merged_labeled.rds"))
}