# Load necessary libraries
library(CountClust)
library(maptpx)
library(Matrix)
library(tools)
library(Seurat)
library(cowplot)
library(MatrixExtra)

# --- Define Command-Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide the number of topics (k) as a command-line argument.")
}
k <- as.numeric(args[1])
if (is.na(k) || k <= 0) {
  stop("The number of topics (k) must be a positive integer.")
}
message(paste("Running fastTopics with k =", k, "topics."))

# --- Load and Prepare Data ---
# Replace 'your_data.mtx' and 'genes.tsv' with your actual file paths.
# The matrix should be a sparse count matrix (e.g., in Matrix Market format).
# The genes file should contain the gene names.
tryCatch({
  # Load the count matrix. Assuming it's in Matrix Market format.
  data <- readRDS("/gstore/data/astroMetaAnalysis/data/DS_50k.rds")
  counts <- GetAssayData(data,slot = "counts",assay = "RNA")
  rm(data)
  gc()
  # Load gene names if available
  # genes <- read.delim("genes.tsv", header = FALSE, stringsAsFactors = FALSE)$V1
  # rownames(counts) <- genes
  
  # --- Recommended Pre-processing ---
  # Filter out lowly expressed genes or cells if needed.
  #remove genes expressed in fewer than 10 cells.
  counts <- counts[rowSums(counts > 0) >= 10, ]
  
}, error = function(e) {
  stop("Error loading data files. Please check the file paths and format.\n", e$message)
})
feat <- as.data.frame(genomitory::getFeatures("GMTY17:GRCm38/GRCm38.IGIS4.0.genes.rds@REVISION-3")) # get gene names
library(dplyr)
protein_coding <- feat %>% filter(type =="protein_coding")

filt.counts <- counts[which(rownames(counts) %in% protein_coding$symbol),]


# --- Fit the Model ---
# The results will be saved in a list.
fit_results <- list()

# 1. Poisson NMF Fit
message("Starting Poisson GOM fit...")
results.dir <- "/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling/results"

FitGoM(as.matrix(t_deep(filt.counts)), K = k, tol = 0.1,
       path_rda = file.path(results.dir,paste0("GoM_k", k, "_results.rda")))

message("Script finished successfully.")