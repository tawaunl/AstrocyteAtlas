# Write Tables and Extract BIC 
library(Seurat)
library(fastTopics)
library(ggplot2)
library(MatrixExtra)
library(grid)
library(gridExtra)

data.dir <- "/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Rscripts/TopicModeling/14topics_tol0.1"
fig.dir <- "/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Rscripts/PaperFigures/SupplementalFigure_TopicModeling"
load("/gstore/data/astroMetaAnalysis/neurodegeneration_meta-analysis/Rscripts/TopicModeling/14topics_tol0.1/FitGoM_50K_k14_tol0.1.rda")
data <- readRDS("/gstore/data/astroMetaAnalysis/data/DS_50k.rds")
usage <- read.csv(file.path(data.dir,"usage_50k.csv"))
theta <- read.csv(file.path(data.dir,"theta_50k.csv"))
scores <- read.csv(file.path(data.dir,"score_min_50k.csv"))

colnames(usage)[1] <- "cell_id"
colnames(theta)[1] <- "gene"
colnames(scores)[1] <- "topic"
n_genes = 25
n_features = 200
prolif_topic = 'lda_14'
prolif_cutoff = 0.08
topic_range <- 1:14
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

pdf(file.path(fig.dir,"TopGenesinTopics.pdf"),width = 11,height=11)
do.call("grid.arrange",list(grobs=topic_plots,top=textGrob("Gene scores for each topic")))
dev.off()
# write gene loadings to CSV -------------
write.csv(theta,file = file.path(fig.dir,"TopicGeneLoadings.csv"))
write.csv(usage,file = file.path(fig.dir,"TopicCellScores.csv"))

topic_plots_umap <- list()
x=1
for(topic in topic_range){

    topic_plots_umap[[x]] <- FeaturePlot(data,paste0("lda_",topic),pt.size = 1.5) +
    ggtitle(paste0("Topic ",x)) +
    theme_void() + theme(plot.title = element_text(hjust = 0.5,face = "bold"), legend.position = "none") 
    x=x+1
}

png(file.path(fig.dir,"UMAPscoringTopics.png"),
    width = 1000, height=1000)
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
