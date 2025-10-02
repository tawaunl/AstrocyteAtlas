library(Seurat)
library(CountClust)
library(ggplot2)
library(MatrixExtra)
library(grid)
library(gridExtra)

# Load data----------------
k_values <- c(2,3, 10,14,16,18)
out.dir <- "/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling"

data <- readRDS("/gstore/data/astroMetaAnalysis/data/DS_50k.rds")
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
k <- 18
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
topic_range <- 1:18 

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
  topic_plots[[x]] <- ggplot(df[1:20,],aes(x=reorder(gene,score),y=score)) + geom_bar(stat="identity") +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0),labels = function(x) format(x, scientific = TRUE),breaks = range(df$score)) + 
    xlab("Gene") + ggtitle(paste0("Topic",x)) +
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


# GSEA on top 50 Topic loadings
library(clusterProfiler)
# First get top 50 loadings from each topic
geneList <- list()
for(topic in 2:15){
  lda_data <- theta[,c(1,topic)]
  lda_data <- lda_data[order(lda_data[,2],decreasing = TRUE),]
  genes <- lda_data[,2]
  names(genes) <-lda_data[,1]
  
  gene.df <- bitr( names(genes), fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org.Mm.eg.db::org.Mm.eg.db)
  x <- match(gene.df$SYMBOL,names(genes))
  gsea <- genes[x]
  names(gsea) <- gene.df$ENTREZID
  geneList[[colnames(lda_data)[2]]] <- gsea
}

ck <- compareCluster(geneList,
                    fun = "gseGO",OrgDb = org.Mm.eg.db::org.Mm.eg.db ,ont="BP")

pdf(file.path(fig.dir,"GSEAonTopics.pdf"),width=10,height = 14)
dotplot(ck, by="NES") + ggtitle("GO:BP Pathways Topic Loadings") + 
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=20)) +
  xlab("LDA Topic")
dev.off()
