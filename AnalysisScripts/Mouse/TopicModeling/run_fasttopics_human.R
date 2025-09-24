# Load necessary libraries
library(fastTopics)
library(Matrix)
library(tools)
library(Seurat)
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
  data <- readRDS("/gstore/data/astroMetaAnalysis/data/HumanSeurat.rds")
  counts <- GetAssayData(data,slot = "counts",assay = "RNA")
  
  # Load gene names if available
  # genes <- read.delim("genes.tsv", header = FALSE, stringsAsFactors = FALSE)$V1
  # rownames(counts) <- genes
  
  # --- Recommended Pre-processing ---
  # Filter out lowly expressed genes or cells if needed.
  #remove genes expressed in fewer than 10 cells.
  counts <- counts[rowSums(counts > 0) >= 50, ]
  
}, error = function(e) {
  stop("Error loading data files. Please check the file paths and format.\n", e$message)
})

# --- Fit the Models ---
# The results will be saved in a list.
fit_results <- list()

# 1. Poisson NMF Fit
message("Starting Poisson NMF fit...")
fit_poisson <- fit_poisson_nmf(X=counts,
                               k = k, numiter = 200, verbose = "progressbar")
fit_results[["poisson_nmf"]] <- fit_poisson

# 2. Multinomial Topic Model Fit
message("Starting Multinomial Topic Model fit...")
fit_multinom <- fit_topic_model(counts, k = k, numiter = 200, verbose = "progressbar")
fit_results[["multinom_topic_model"]] <- fit_multinom

# --- Save the Results ---
# Define output filename based on the number of topics (k)
results.dir <- "/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling/results"
output_filename <- file.path(results.dir,paste0("fasttopics_k", k, "_results.rds"))
message(paste("Saving results to", output_filename))

tryCatch({
  saveRDS(fit_results, file = output_filename)
}, error = function(e) {
  stop("Error saving the results. Please check file permissions or disk space.\n", e$message)
})

message("Script finished successfully.")