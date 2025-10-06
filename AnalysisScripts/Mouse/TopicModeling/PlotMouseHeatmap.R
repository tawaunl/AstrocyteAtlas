# Load necessary libraries
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(CountClust)
# Load Data-------- ---

# Main expression matrix (e.g., Log-Fold Change)
data <- readRDS("/gstore/data/astroMetaAnalysis/data/DS_50k.rds")
bulk <- AggregateExpression(data,assays = "RNA",group.by =c("Sample","finalClusters"),
                            return.seurat = T)
counts <- GetAssayData(bulk,slot = "counts",assay = "RNA")

out.dir <- "/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling"
k <- 14
topic_range <- 1:k

rda_fname <- file.path(out.dir,"results",paste0("GoM_k",k,"_results.rda"))
load(rda_fname)
usage <- read.csv(file.path(out.dir,"MouseResults",paste0("usage_k",k,".csv")))
theta <- read.csv( file.path(out.dir,"MouseResults",paste0("theta_k",k,".csv")))

fullScores <- ExtractTopFeatures(theta = Topic_clus$theta,
                                 top_features = dim(Topic_clus$theta)[1],
                                 method ="poisson",options = "min",shared = T)
colnames(usage)[1] <- "cell_id"
colnames(theta)[1] <- "gene"
genes <- theta$gene
num_genes <- length(genes)
num_topics <- nrow(fullScores$indices)

scores_matrix <- matrix(
  data = NA,
  nrow = num_genes,
  ncol = num_topics,
  dimnames = list(
    genes, # Set row names to your gene names
    paste0("Topic_", 1:num_topics) # Set column names
  )
)

for (i in 1:num_topics) {
  cluster_indices <- fullScores$indices[i, ]
  cluster_scores  <- fullScores$scores[i, ]
  scores_matrix[cluster_indices, i] <- cluster_scores
}

scores_matrix <- scores_matrix[-which(startsWith(genes,"Rps") |startsWith(genes,"Rpl")),]

plot_score_dict <- list()
for (topic in topic_range) {
  topic_scores = data.frame(gene=rownames(scores_matrix),
                            scores=scores_matrix[,topic])
  topic_dict = list()
  topic_scores <- topic_scores[order(topic_scores$scores,decreasing = T),]
  topic_dict[topic_scores$gene[1:20]] <- topic_scores$scores[1:20]
  plot_score_dict[[topic]] = topic_dict
}

names(plot_score_dict) <- paste0("lda_", 1:14)


# Convert the list of top genes into a single, clean data frame
# The `bind_rows` function with `.id = "topic_group"` is perfect for this
gene_groups_df <- enframe(unlist(plot_score_dict), name = "name", value = "score") %>%
  separate(name, into = c("topic_group", "gene"), sep = "\\.")
# Now gene_groups_df looks like this:
#   topic_group gene     score
# 1 lda_1       Gene454  0.040
# 2 lda_1       Gene23   0.015

top_genes_df <- gene_groups_df %>%
  group_by(topic_group) %>%
  slice_max(order_by = score, n = 15) %>%
  ungroup()
gene_groups_df <- top_genes_df
# 1. Find the dominant topic for each cell from the 'usage' data frame
rownames(usage) <- usage$cell_id
usage <- usage[,-1]
dominant_topic <- colnames(usage)[apply(usage, 1, which.max)]
names(dominant_topic) <- rownames(usage)

levels(data$finalClusters) <- c("0","1","2","3","4","0")

cells <- rownames(usage) # Get cell IDs from your usage matrix
cell_metadata <- data.frame(
  cell_id = cells,
  clusters = data$finalClusters
)

cell_metadata$pseudobulk_group <-cell_metadata$clusters


# 3. (Recommended) Normalize the pseudobulk matrix and scale

log_normalized_pseudobulk <- log1p(sweep(counts, 2, 1e4 / colSums(counts), `*`))
common_genes <- intersect(gene_groups_df$gene, rownames(log_normalized_pseudobulk))
log_normalized_pseudobulk <- log_normalized_pseudobulk[common_genes, ]
expression_matrix_scaled <- t(scale(t(log_normalized_pseudobulk)))

# 1. Create the column annotation data frame for the pseudobulks
# The rownames must match the column names of our final matrix
pseudobulk_metadata <- data.frame(
  group = colnames(expression_matrix_scaled)
) %>%
  separate(group, into = c("cell","clusters"), sep = "_", remove = F)

rownames(pseudobulk_metadata) <- pseudobulk_metadata$group

pseudobulk_metadata$clusters <- factor(pseudobulk_metadata$clusters)

levels(pseudobulk_metadata$clusters) = c("0","1","2","4","0")

# 3. Filter and order your gene annotation to match the final matrix
gene_groups_df_filtered <- gene_groups_df %>% 
  filter(gene %in% rownames(expression_matrix_scaled)) %>%
  arrange(topic_group)

expression_matrix_final <- expression_matrix_scaled[gene_groups_df_filtered$gene, ]

color_list <- list(
  clusters = c("0" = "dodgerblue", "1" = "goldenrod1","2"="red", "4" = "darkgreen") # Add cluster colors
)

# Update the HeatmapAnnotation object to include 'cluster'
col_annotation <- HeatmapAnnotation(
  clusters = pseudobulk_metadata[, c("clusters")], # Add cluster here
  col = color_list
)


ht <- Heatmap(
  expression_matrix_final,
  name = "Scaled Exp",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = T,
  cluster_columns = T,
  top_annotation = col_annotation,
  row_split = gene_groups_df_filtered$topic_group,
  column_split = pseudobulk_metadata$cluster,
  # --- AESTHETICS ---
  show_row_names = FALSE,
  show_column_names = F,
  show_row_dend = F,
  show_column_dend = F,
  column_names_rot = 45,
  row_title = NULL
)
ht
scores_matrix <- scores_matrix[rownames(expression_matrix_final), ]

# This object is the annotation heatmap for the topic weights
ht_weights <- Heatmap(
  t(scale(t(scores_matrix))),
  name = "Weight", # Legend title for this heatmap
  col = colorRamp2(c(-1,0,1), c("lightblue","white", "firebrick4")),
  cluster_columns = F,
  row_split = gene_groups_df_filtered$topic_group,
  # Aesthetics for this specific heatmap
  show_row_names = FALSE,
  width = unit(4, 'cm'), # Adjust the width of this heatmap
  column_title = "Topic Weights" # Title that appears above this heatmap
)

draw(ht_weights + ht)
