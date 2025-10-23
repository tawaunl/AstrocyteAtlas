library(Seurat)
library(CountClust)
library(ggplot2)
library(MatrixExtra)
library(grid)
library(gridExtra)

# Load data----------------
k_values <- c(2,3, 10,14,16,18)
out.dir <- "/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling"

out <- readRDS(file.path(out.dir,"HumanResults","compGoM_results.rds"))
names(out) <- paste0("K_", k_values)
bic.plot <- sapply(names(out), function(x) out[[x]]$BIC)
bic.plot <- data.frame(bic.plot)
bic.plot$K <- rownames(bic.plot)
bic.plot$K <- factor(bic.plot$K, levels=names(out))
png(file.path(out.dir,"HumanResults","BIC_plot.png"), width=800, height=600)
ggplot(bic.plot, aes(x=K, y=bic.plot)) +
  geom_point() + geom_line() +
  xlab("Number of Topics (K)") + ylab("BIC") +
  theme_bw(base_size=18)
dev.off()
# Plot results for K=14----------------
k <- 10
data <- readRDS("/gstore/data/astroMetaAnalysis/data/HumanSubDataforTopics.rds")

usage <- read.csv(file.path(out.dir,"HumanResults",paste0("usage_k",k,".csv")))
theta <- read.csv( file.path(out.dir,"HumanResults",paste0("theta_k",k,".csv")))
scores <- read.csv(file.path(out.dir,"HumanResults",paste0("score_min_k",k,".csv")))

colnames(usage)[1] <- "cell_id"
colnames(theta)[1] <- "gene"
colnames(scores)[1] <- "topic"
n_genes = 25
n_features = 200
prolif_topic = 'lda_18'
prolif_cutoff = 0.08
topic_range <- 1:k

# join Usage to Seurat object
for (topic in topic_range) {
  data[[paste0("lda_",topic)]] <- usage[[paste0("lda_",topic)]]
  
}

# Prepare list for scores
topic_score_dict = list()
plot_score_dict = list()

for (topic in topic_range) {
  topic_scores = scores[topic,]
  
  topic_dict = list()
  for (i in 1:n_features) {
    gene_index <- as.numeric(topic_scores[paste0("indices.",i)])
    gene <- theta$gene[gene_index]
    score = as.numeric(topic_scores[paste0("scores.",i)])
    topic_dict[gene] = score
  }
  topic_score_dict[[topic]] = topic_dict
  
  plot_series = lapply(topic_dict,sort,decreasing=TRUE)
  plot_score_dict[[topic]] = plot_series
}

#Plot
topic_plots <- list()
x=1
for(topic in topic_range){
  df <- data.frame(gene=names(topic_score_dict[[topic]]),score=unlist(topic_score_dict[[topic]]))
  topic_plots[[x]] <- ggplot(df[1:20,],aes(x=reorder(gene,score),y=log(score+1)))+ geom_bar(stat="identity") +
    coord_flip() +
    scale_y_sqrt() + 
    xlab("Gene") +ylab("Topic Weight")+ ggtitle(paste0("Topic ",x)) +
    theme_light()+ theme(plot.title = element_text(hjust = 0.5)) 
  x=x+1
}


pdf(file.path(out.dir,"HumanResults","TopGenesinTopics.pdf"),width = 11,height=11)
do.call("grid.arrange",list(grobs=topic_plots,top=textGrob("Gene scores for each topic")))
dev.off()

topic_plots_umap <- list()
x=1
for(topic in topic_range){
  
  topic_plots_umap[[x]] <- FeaturePlot(data,paste0("lda_",topic),pt.size = 1.5) +
    ggtitle(paste0("Topic ",x)) +
    theme_void() + theme(plot.title = element_text(hjust = 0.5,face = "bold"), legend.position = "none") 
  x=x+1
}

png(file.path(file.path(out.dir,"HumanResults","UMAPscoringTopics.png")),
    width = 2000, height=2000)
do.call("grid.arrange",list(grobs=topic_plots_umap,
                            top=textGrob(expression(bold(underline("UMAP embedding for each topic"))))))
dev.off()


# GSEA on Topic loadings---------
library(clusterProfiler)
#Updated Pathways -------------------
rda_fname <- file.path(out.dir,"results",paste0("GoM_Human_k",k,"_results.rda"))
load(rda_fname)

fullScores <- ExtractTopFeatures(theta = Topic_clus$theta,
                                 top_features = dim(Topic_clus$theta)[1],
                                 method ="bernoulli",options = "min",shared = F)

genes <- theta$gene
num_genes <- length(genes)
num_topics <- nrow(fullScores$indices)

scores_matrix <- matrix(
  data = 0,
  nrow = num_genes,
  ncol = num_topics,
  dimnames = list(
    genes, # Set row names to your gene names
    paste0("Topic_", 1:num_topics) # Set column names
  )
)
for (i in 1:num_topics) {
  cluster_indices <- fullScores$indices[i, ]
  cluster_indices <- cluster_indices[!is.na(cluster_indices)]
  cluster_scores  <- fullScores$scores[i, ]
  cluster_scores <- cluster_scores[!is.na(cluster_scores)]
  scores_matrix[cluster_indices, i] <- cluster_scores
}
geneList <- list()
for(topic in 1:k){
  lda_data <- data.frame(scores=scores_matrix[,topic]) %>% arrange(desc(scores))
  
  genes <- lda_data$scores
  names(genes) <- rownames(lda_data)
  
  gene.df <- bitr( names(genes), fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  x <- match(gene.df$SYMBOL,names(genes))
  gsea <- genes[x]
  names(gsea) <- gene.df$ENTREZID
  geneList[[paste0("Topic_",topic)]] <- gsea[unique(names(gsea))]
}
## Reactome Enrichment ------------------
out.dir <- "/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling/HumanResults"

ck <- lapply(geneList, function(x){
  ReactomePA::enrichPathway(names(x[x>0]), organism = "human",
                            pvalueCutoff = .1, maxGSSize = 1000,readable = FALSE)
})

saveRDS(ck, file.path(out.dir,"HumanReactomeEnrichment_results.rds"))
ck <- readRDS(file.path(out.dir,"HumanReactomeEnrichment_results.rds"))

gsea_results <- lapply(names(ck), function(cluster){
  ck[[cluster]]@result
})
names(gsea_results) <- names(ck)
# Assuming gsea_results is a list of GSEA results per cluster
gsea_results <- gsea_results[topic_range]
source("~/scHelpers.R")

pdf(file.path(out.dir,"ReactomeEnrichmentonTopics_Human.pdf"),width=10,height = 10)
DotPlotCompare(
  gsea_list = gsea_results,
  n = 3,
  size_col = "RichFactor",
  color_col = "pvalue",
  size_cutoff = NULL,      # Filter to pathways with NES >= 1.5
  color_cutoff = 0.1,
  direction = "positive" # Filter to pathways with p.adjust <= 0.05
)
dev.off()

##GSEA Reactome ------------------
ck <- lapply(geneList, function(x){
  ReactomePA::gsePathway(x[x>0], organism = "human",seed = 824,
                            pvalueCutoff = .1, maxGSSize = 1000)
})

saveRDS(ck, file.path(out.dir,"HumanReactomeGSEA_results.rds"))
ck <- readRDS(file.path(out.dir,"HumanReactomeGSEA_results.rds"))

gsea_results <- lapply(names(ck), function(cluster){
  ck[[cluster]]@result
})
names(gsea_results) <- names(ck)
# Assuming gsea_results is a list of GSEA results per cluster
gsea_results <- gsea_results[topic_range]
source("~/scHelpers.R")

pdf(file.path(out.dir,"ReactomeGSEAonTopics_Human.pdf"),width=10,height = 10)
DotPlotCompare(
  gsea_list = gsea_results,
  n = 3,
  size_col = "NES",
  color_col = "pvalue",
  size_cutoff = NULL,      # Filter to pathways with NES >= 1.5
  color_cutoff = 0.1,
  direction = "positive" # Filter to pathways with p.adjust <= 0.05
)
dev.off()

## GSEA GO:BP----------
ck <- lapply(geneList, function(x){
  gseGO(x[x>0], OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "BP",
                            pvalueCutoff = .1, maxGSSize = 1000)
})

saveRDS(ck, file.path(out.dir,"HumanGO_BP_GSEA_results.rds"))
ck <- readRDS(file.path(out.dir,"HumanGO_BP_GSEA_results.rds"))

gsea_results <- lapply(names(ck), function(cluster){
  ck[[cluster]]@result
})
names(gsea_results) <- names(ck)
# Assuming gsea_results is a list of GSEA results per cluster
gsea_results <- gsea_results[topic_range]
source("~/scHelpers.R")

pdf(file.path(out.dir,"GO_BP_GSEAonTopics_Human.pdf"),width=8,height = 8)
DotPlotCompare(
  gsea_list = gsea_results,
  n = 3,
  size_col = "NES",
  color_col = "pvalue",
  size_cutoff = NULL,      # Filter to pathways with NES >= 1.5
  color_cutoff = 0.01,
  direction = "positive" # Filter to pathways with p.adjust <= 0.05
)
dev.off()

## Enrichment GO:BP----------
ck <- lapply(geneList, function(x){
  enrichGO(names(x[x>0]), OrgDb = org.Hs.eg.db::org.Hs.eg.db,ont = "BP",
        pvalueCutoff = .1, maxGSSize = 1000,readable = TRUE)
})

saveRDS(ck, file.path(out.dir,"HumanGO_BP_Enrichment_results.rds"))
ck <- readRDS(file.path(out.dir,"HumanGO_BP_Enrichment_results.rds"))

gsea_results <- lapply(names(ck), function(cluster){
  ck[[cluster]]@result
})
names(gsea_results) <- names(ck)
# Assuming gsea_results is a list of Enrichment results per cluster
gsea_results <- gsea_results[topic_range]
source("~/scHelpers.R")

pdf(file.path(out.dir,"GO_BP_EnrichmentonTopics_Human.pdf"),width=8,height = 8)
DotPlotCompare(
  gsea_list = gsea_results,
  n = 3,
  size_col = "Count",
  color_col = "pvalue",
  size_cutoff = NULL,      # Filter to pathways with NES >= 1.5
  color_cutoff = 0.01,
  direction = "positive" # Filter to pathways with p.adjust <= 0.05
)
dev.off()
