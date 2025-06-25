# Load necessary libraries
library(Seurat)
library(edgeR)
library(dplyr)
library(scater)

dir <- "~/Documents/AstrocytePaper"
data <- readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))

# Subset the Seurat object to include only Cortex and Hippocampus cells
sub <- subset(
  data,
  subset = BrainRegion %in% c("Cortex", "Hippocampus")& ClusterNames %in% c("DAA1")
)
# Extract metadata from Seurat object
meta.data <- sub@meta.data

# Get raw counts from the Seurat object
counts <- GetAssayData(sub, slot = "counts")

# Split metadata by `study` to compute per-study gene expression
gene_expression_filter <- meta.data %>%
  split(.$StudyName) %>%  # Group cells by study
  lapply(function(meta_subset) {
    # Subset the counts for cells in this study
    subset_counts <- counts[, rownames(meta_subset)]  # Extract counts for these cells
    # Calculate the proportion of cells with non-zero expression per gene
    gene_proportions <- rowSums(subset_counts > 0) / ncol(subset_counts)
    # Return a logical vector indicating genes expressed in at least % of cells
    gene_proportions >= 0.025 }) %>%
  as.data.frame()

# Identify genes expressed in at least % of cells in **all studies**
genes_expressed_in_all_studies <- rownames(counts)[rowSums(gene_expression_filter) == ncol(gene_expression_filter)]
# Filter counts to only include these genes before proceeding
se <- SingleCellExperiment(assays = list(counts = counts),
                           colData = sub@meta.data)
summed <- aggregateAcrossCells(se, 
                               id=colData(se)[,c("Sample","BrainRegion")])

summed.filt <- summed[,summed$ncells >= 10]
summed.filt <- summed.filt[genes_expressed_in_all_studies,]
current <- summed.filt
y <- DGEList(counts(current), samples=colData(current))
keep <- filterByExpr(y, group=current$BrainRegion)
y <- y[keep,]
y <- calcNormFactors(y)

design <- model.matrix(~0+factor(BrainRegion), y$samples)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
colnames(design) <- c("Cortex","Hippocampus")
my.contrasts <- makeContrasts(CortexvsHipp=Cortex-Hippocampus, levels = design)

res <- glmQLFTest(fit, coef=ncol(design),contrast = my.contrasts)

res <- topTags(res,n=60000,p.value = 2)
res <- res$table

res <- res %>% mutate(sig = ifelse((FDR <= 0.05 ) & abs(logFC) >= 1, "yes", "no"))

# pPlot Volcano ------------
source("~/Documents/scHelpers.R")

labels <- rownames(res)[which(res$sig=="yes")]
edgeRVolcano(res = res,plot.title = "Cortex vs. Hippocampus",labels = labels)
