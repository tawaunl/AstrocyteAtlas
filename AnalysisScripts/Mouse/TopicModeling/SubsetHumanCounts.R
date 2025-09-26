counts <- readRDS("/gstore/data/astroMetaAnalysis/data/HumanCounts_subsetTopics.rds")
counts <- GetAssayData()
names(cell_clusters) <- paste0("Cell", 1:1000)

# --- Step 2: Sample a fixed percentage of cells from each cluster ---
# Define the percentage (as a fraction) to sample from each cluster.
fraction_to_sample <- 0.25 # Sample 25% of cells

# Split the cell names by their cluster.
cells_by_cluster <- split(names(cell_clusters), f = cell_clusters)

# Use lapply to sample the defined fraction from each cluster.
set.seed(824) # For reproducible sampling
sampled_cells_list_frac <- lapply(cells_by_cluster, function(cell_names) {
  num_current_cluster <- length(cell_names)
  # Calculate the number of cells to sample by multiplying by the fraction.
  num_to_sample <- floor(num_current_cluster * fraction_to_sample)
  
  sample(cell_names, size = num_to_sample, replace = FALSE)
})

# Combine the lists of sampled cells into a single vector.
sampled_cells_frac <- unlist(sampled_cells_list_frac)

# --- Step 3: Subset your data ---
# Subset the data based on the new vector of sampled cells.
sampled_counts_frac <- counts[, sampled_cells_frac]
sampled_clusters_frac <- cell_clusters[sampled_cells_frac]
saveRDS(sampled_counts_frac,"/gstore/data/astroMetaAnalysis/data/HumanCounts_subsetTopics.rds")
# Check the new dimensions and cluster distribution.
dim(sampled_counts_frac)
table(sampled_clusters_frac)