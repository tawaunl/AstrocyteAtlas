# Reviewewr Response

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(wesanderson)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 
library(scales)
library(ggrepel)
library(scran.chan)
library(scuttle)
library(edgeR)
source("~/Documents/scHelpers.R")
library(openxlsx)

# Load Data ====================================================================
dir <- "~/Documents/AstrocytePaper"
data <- readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))

myoc.markers <- read.xlsx("~/Downloads/MyocMarkers.xlsx")
myoc.markers <- myoc.markers[order(myoc.markers$avg_log2FC,decreasing = T),]
mm <- myoc.markers$gene[1:30]
all <- table(data$Sample)
all <- data.frame(all)
colnames(all) <- c("Sample","Total")

# Score cells with Myoc Markers -------------
score_intersect <- intersect(mm,rownames(data))
scores_by_cell <- colMeans(
  GetAssayData(data,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
data[["MyocScore"]] <- scores_by_cell

data$ClusterNames <- factor(data$finalClusters)
levels(data$ClusterNames) <-c("Homeostatic","DAA1","DAA2","Synapse")

scvi_coords <- get_scvi_coords(data,data$finalClusters)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))
text <- data$ClusterNames
text_x <- vapply(split(scvi_coords$UMAP1, text), median, FUN.VALUE=0)
text_y <- vapply(split(scvi_coords$UMAP2, text), median, FUN.VALUE=0)
## Plot F -----------------------

ggplot(scvi_coords,aes(x=UMAP1, y=UMAP2, colour=MyocScore)) +
  geom_point() +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) +
  scale_colour_gradientn(limits= c(0,1), oob=squish,
                         colours =  RColorBrewer::brewer.pal(5,"Purples"))
DAA1 <- subset(data, subset = ClusterNames=="DAA1")

DAA1 <- Split_Layers(DAA1,split.by = "StudyID")


DAA1 <- runSeurat(DAA1)

DAA1 <- IntegrateLayers(object = DAA1, method = CCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose=FALSE,k.weight=40)

# re-join layers after integration
DAA1[["RNA"]] <- JoinLayers(DAA1[["RNA"]])

DAA1 <- FindNeighbors(DAA1, reduction = "integrated.cca", dims = 1:30)
DAA1 <- FindClusters(DAA1, resolution = c(0.1,0.15,0.20,.25,0.3,0.35,0.4,0.45,0.5))
DAA1 <- RunUMAP(DAA1,reduction = "integrated.cca",dims=1:30)
DimPlot(DAA1,group.by = "RNA_snn_res.0.2",label=T,label.size = 10.5,pt.size = 1)
DimPlot(DAA1,group.by = "StudyID")
DAA1 <- Add_Cell_QC_Metrics(object = DAA1, species = "mouse")

DotPlot(DAA1,features = "Myoc",group.by ="RNA_snn_res.0.25" )
DotPlot(DAA1,features = "Crym",group.by ="RNA_snn_res.0.25" )
DotPlot(DAA1,features = "Meg3",group.by ="RNA_snn_res.0.25" )

coldata <- as.data.frame(DAA1@meta.data)
study_cluster_tmp_by_study <- coldata %>%
  group_by(StudyName, RNA_snn_res.0.15) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

df_normalized <- study_cluster_tmp_by_study %>%
  group_by(RNA_snn_res.0.15) %>%
  mutate(TotalCellsInCluster = sum(freq),
         NormalizedCells = freq / TotalCellsInCluster)
#### Plot study Composition -----------

ggplot(df_normalized, aes(fill=StudyName, y=NormalizedCells, x=RNA_snn_res.0.15)) +
  geom_bar(position="stack", stat="identity", alpha=0.6) + theme_void()+ theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14)) +
  xlab("Clusters") + ylab("Proportion of cells within a study") +
  theme(axis.title.y =  element_text(angle = 90))



QC_Plots_Combined_Vln(seurat_object = DAA1, pt.size = 0.1,group.by ="RNA_snn_res.0.15",
                      colors_use = scales::hue_pal()(7))


counts <- round(GetAssayData(DAA1, slot="counts", assay="RNA"))   
counts <- initializeSparseMatrix(counts)
lognorm <- logNormCounts.chan(counts,batch = DAA1$StudyID)
markers <- scoreMarkers.chan(lognorm,
                             groups=DAA1$RNA_snn_res.0.25, batch=DAA1$StudyID,lfc=0)
saveRDS(DAA1,"~/Documents/AstrocytePaper/DAA1subcluster.rds")
# Get top 5 markers from each cluster 
topmarkers <- c()
for (cluster in 1:length(markers$statistics)) {
  clustermarkers <- data.frame(markers$statistics[[cluster]])
  clustermarkers <- clustermarkers[order(clustermarkers$logFC,decreasing = TRUE),]
  top <- rownames(clustermarkers)[1:5]
  topmarkers <- c(topmarkers,top)
}

topmarkers[23:25] <- c("Nrg3","Opcml","Ptprd")

topmarkers <- c(topmarkers,"Thbs4")

plot <- DotPlot(DAA1,features = c(unique(topmarkers)), assay = "RNA",
                cols = "RdYlBu",group.by = "RNA_snn_res.0.25",scale = TRUE)


plot <- DotPlot(DAA1,features = c("Myoc","Thbs4","Ttr","Crym","Vim"), assay = "RNA",
                cols = "RdYlBu",group.by = "RNA_snn_res.0.25",scale = TRUE)
plot + RotatedAxis() +coord_flip() + xlab('Gene') +  ylab('Cluster') +
  theme(axis.text = element_text(size = 14),axis.title = element_text(size=18,face='bold'))  

sce <- DAA1 %>% 
  Seurat::as.SingleCellExperiment() 

library(scran)
library(purrr)
library(dplyr)
library(AnnotationDbi)
features <- markers$statistics %>% 
  map(~ .x |>
        as.data.frame() |>
        arrange(desc(logFC))|>
        dplyr::slice(1:6) |>
        rownames())|> 
  unlist2()

features[26:28] <- c("Nrg3","Opcml","Ptprd")

sce$subclusters <- factor(sce$subclusters,levels = c("Myoc+","Reactive1",
                                                     "Reactive2","Interferon-Responsive",
                                                     "Meg3+"))

rowData(sce)$Marker <- features[match(rownames(sce), features)] |>
  names() %>% 
  factor(levels = levels(sce$subclusters))


sce %>% 
  scDotPlot::scDotPlot(features = unique(features),
                       group = "subclusters",
                       #block = "Sample",
                       scale = TRUE,
                       cluster = FALSE,
                       groupAnno = "subclusters",
                       featureAnno = "Marker",
                       featureLegends = FALSE,
                       annoHeight = 0.025,
                       annoColors = list("subclusters" = scales::hue_pal()(7),
                                                            "Marker" = scales::hue_pal()(7)),
                       annoWidth = 0.1,fontSize = 20,
                       dotColors=c( "blue","#FFFFBF", "red")) 



## Renaming DAA1 subclusters ---------
DAA1$subclusters <- DAA1$RNA_snn_res.0.25
DAA1$subclusters <- gsub("\\<0\\>", "Myoc+", DAA1$subclusters)
DAA1$subclusters <- gsub("\\<1\\>", "Reactive1", DAA1$subclusters)

DAA1$subclusters <- gsub("\\<2\\>", "Reactive2", DAA1$subclusters)
#DAA1$subclusters <- gsub("\\<5\\>", "Reactive1", DAA1$subclusters)
#DAA1$subclusters <- gsub("\\<3\\>", "Reactive1", DAA1$subclusters)

DAA1$subclusters <- gsub("\\<3\\>", "Interferon-Responsive", DAA1$subclusters)
DAA1$subclusters <- gsub("\\<4\\>", "Meg3+", DAA1$subclusters)



BRumaps <- lapply(unique(DAA1$BrainRegion), function(x){
  
  plot_umap_library(sobj=DAA1,meta="BrainRegion",sample=x)
})
names(BRumaps) <- unique(DAA1$BrainRegion)

plot_grid(BRumaps$Hippocampus,BRumaps$Cortex,BRumaps$CNS,BRumaps$`Corpus Callosum`,
          ncol = 2)


abundances<- table(DAA1$RNA_snn_res.0.25,factor(DAA1$Sample))
df_long <- GetAbundanceTable(abundances,DAA1@meta.data,samplecol = "Sample",
                             subset.cols = c("Sample","Disease","BrainRegion","StudyName"),
                             calc.col = "RNA_snn_res.0.25")
level_order <- c("AD_control","AD","MS_control","MS")
df_long$Disease <- factor(df_long$Disease,levels = level_order) 

fill.codes <- c("white","darkgreen","white", "purple")
color.codes <- c("darkgreen","darkgreen","purple","purple")
zone <- levels(factor(df_long$Disease))

filter <- c("Zhao 5XFAD","Choi 5XFAD","Shi P301S")

df_long <- df_long[! df_long$StudyName %in% filter,]

ggplot(df_long, aes(x=Disease, y=Percent,fill=Disease,colour = Disease)) + theme_classic() +
  geom_jitter(color="darkgrey",width = 0.4,size=3,alpha=0.6)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  facet_wrap(~Cluster,scales = "free",ncol=3) + xlab("Disease Label") +
  theme(axis.text.y = element_text(size=20)) +
  ylab("Cellularity Proportion") +  
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(values=setNames(fill.codes, zone))+
  scale_color_manual(values=setNames(color.codes, zone))

# LPS scoring ----------
lps.de <- read.csv("~/Documents/AstrocytePaper/Figure3/LPSvsSalineDE.csv",row.names = 1)
lps.genes <- rownames(lps.de[lps.de$diffexpressed=="LPS",])
score_intersect <- intersect(lps.genes,rownames(DAA1))
scores_by_cell <- colMeans(
  GetAssayData(DAA1,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
DAA1[["LPSCurrentAstroScore"]] <- scores_by_cell


scvi_coords <- get_scvi_coords(DAA1,DAA1$subclusters)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))
text <- DAA1$subclusters
text_x <- vapply(split(scvi_coords$UMAP1, text), median, FUN.VALUE=0)
text_y <- vapply(split(scvi_coords$UMAP2, text), median, FUN.VALUE=0)

ggplot(scvi_coords%>%
         arrange(LPSCurrentAstroScore),aes(x=UMAP1, y=UMAP2, colour=LPSCurrentAstroScore)) +
  geom_point() +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) +
  scale_colour_gradientn(limits= c(0,1), oob=squish,
                         colours =  RColorBrewer::brewer.pal(5,"Purples"))
dev.off()


#AD Scoring -----------------------
AD.DE <- readRDS(file.path(dir, "res_dsl_AD_full.011824.rds"))
MS.DE <- readRDS(file.path(dir, "res_dsl_MS_full.011824.rds"))
AD.DE_plot <- generatePlotTable(AD.DE)
AD.DE_plot$diffexpressed <- "unchanged"

AD.DE_plot$diffexpressed[AD.DE_plot$dl_mu >= 0.25 &
                           AD.DE_plot$FDR<=0.05 &
                           AD.DE_plot$sig=="yes" ] <- "up"

AD.DE_plot$diffexpressed[AD.DE_plot$dl_mu <= -0.25 &
                           AD.DE_plot$FDR<=0.05 &
                           AD.DE_plot$sig=="yes" ] <- "down"

AD.genes <-AD.DE_plot$symbol[AD.DE_plot$diffexpressed=="up"]
score_intersect <- intersect(AD.genes,rownames(DAA1))
scores_by_cell <- colMeans(
  GetAssayData(DAA1,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
DAA1[["ADCurrentAstroScore"]] <- scores_by_cell


scvi_coords <- get_scvi_coords(DAA1,DAA1$subclusters)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))
text <- DAA1$subclusters
text_x <- vapply(split(scvi_coords$UMAP1, text), median, FUN.VALUE=0)
text_y <- vapply(split(scvi_coords$UMAP2, text), median, FUN.VALUE=0)

ggplot(scvi_coords%>%
         arrange(ADCurrentAstroScore),aes(x=UMAP1, y=UMAP2, colour=ADCurrentAstroScore)) +
  geom_point() +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) +
  scale_colour_gradientn(limits= c(0,.5), oob=squish,
                         colours =  RColorBrewer::brewer.pal(9,"Purples"))

# MS Scoring --------------
MS.DE_plot <- generatePlotTable(MS.DE)
MS.DE_plot$diffexpressed <- "unchanged"

MS.DE_plot$diffexpressed[MS.DE_plot$dl_mu >= 0.25 &
                           MS.DE_plot$FDR<=0.05 &
                           MS.DE_plot$sig=="yes" ] <- "up"

MS.DE_plot$diffexpressed[MS.DE_plot$dl_mu <= -0.25 &
                           MS.DE_plot$FDR<=0.05 &
                           MS.DE_plot$sig=="yes" ] <- "down"

MS.genes <-MS.DE_plot$symbol[MS.DE_plot$diffexpressed=="up"]
score_intersect <- intersect(MS.genes,rownames(DAA1))
scores_by_cell <- colMeans(
  GetAssayData(DAA1,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
DAA1[["MSCurrentAstroScore"]] <- scores_by_cell


scvi_coords <- get_scvi_coords(DAA1,DAA1$subclusters)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))
text <- DAA1$subclusters
text_x <- vapply(split(scvi_coords$UMAP1, text), median, FUN.VALUE=0)
text_y <- vapply(split(scvi_coords$UMAP2, text), median, FUN.VALUE=0)

ggplot(scvi_coords%>%
         arrange(MSCurrentAstroScore),aes(x=UMAP1, y=UMAP2, colour=MSCurrentAstroScore)) +
  geom_point() +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) +
  scale_colour_gradientn(limits= c(0,1), oob=squish,
                         colours =  RColorBrewer::brewer.pal(9,"Purples"))

VlnPlot(DAA1,features= c("LPSCurrentAstroScore","ADCurrentAstroScore",
                         "MSCurrentAstroScore"),
        group.by = "subclusters")
# Abundance plot -----------------------

subDA1 <- subset(DAA1,subset= RNA_snn_res.0.25 %in% c(0,3),invert=T)
abundances<- table(DAA1$subclusters,factor(DAA1$Sample))

df_long <- GetAbundanceTable(abundances,coldata = DAA1@meta.data,samplecol = "Sample",
                             subset.cols = c("Sample","Disease","BrainRegion","StudyName"),
                             calc.col = "subclusters")
df_long$Disease <- factor(df_long$Disease,levels = level_order) 

ggplot(df_long, aes(x=Disease, y=Percent,fill=Disease,colour = Disease))+ theme_classic() +
  geom_jitter(color="darkgrey",width = 0.4,size=3,alpha=0.6)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  facet_wrap(~Cluster,scales = "free",ncol=3) + xlab("Disease Label") +
  theme(axis.text.y = element_text(size=20)) +
  ylab("Cellularity Proportion") +  
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(values=setNames(fill.codes, zone))+
  scale_color_manual(values=setNames(color.codes, zone))

df_transformed <- GetAbundanceTable(abundances,DAA1@meta.data,samplecol = "Sample",
                                    subset.cols = c("Sample","Disease","BrainRegion","StudyName"),
                                    calc.col = "subclusters",return_clr = T)

sig.clusters <- list()
for(cluster in 1:length(levels(factor(DAA1$subclusters)))){
  clusterName <- levels(factor(DAA1$subclusters))[cluster]
  dataCluster <- df_transformed[which(df_transformed$Cluster==clusterName),]
  res <- kruskal.test(Percent ~ Disease, data = dataCluster)
  # check if there is any significance
  if(res$p.value < 0.8/length(levels(factor(DAA1$subclusters)))){
    sig.clusters[[clusterName]] <- res
    # do a Welch's t-test to find which specific groups are different
    res.AD <- wilcox.test(dataCluster$Percent[which(dataCluster$Disease=="AD_control")],
                          dataCluster$Percent[which(dataCluster$Disease=="AD")])
    res.MS <- wilcox.test(dataCluster$Percent[which(dataCluster$Disease=="MS_control")],
                          dataCluster$Percent[which(dataCluster$Disease=="MS")])
    if(res.MS$p.value < 0.9){
      sig.clusters[[clusterName]][["MS"]] <- res.MS
    }
    if(res.AD$p.value < 0.9){
      sig.clusters[[clusterName]][["AD"]] <- res.AD
    }
  }
}


VlnPlot(DAA1,"Thbs4",group.by = "RNA_snn_res.0.25")

DimPlot(DAA1,group.by = "subclusters",pt.size=1.5,label=T,label.size = 8,repel = T)


#Percent of cells expressing Myoc---------
# Extract the expression data
expression_data <- FetchData(data, vars = "Myoc")

# Categorize cells as "Expressing" or "Non-Expressing"
expression_data$ExpressionCategory <- ifelse(expression_data$Myoc > 0, "Expressing", "Non-Expressing")

# Count the number of cells in each category
count_table <- table(expression_data$ExpressionCategory)
proportion_table <- prop.table(count_table)

# Convert to a data frame for plotting
plot_data <- as.data.frame(proportion_table)
colnames(plot_data) <- c("Category", "Proportion")

# Add percentages to the plot data
plot_data$Percent <- plot_data$Proportion * 100

# Create a bar plot
ggplot(plot_data, aes(x = Category, y = Percent, fill = Category)) + 
  geom_bar(stat = "identity", width = 0.6)  + 
  labs(title = "Proportion of DAA1 Cells Expressing 'Myoc' > 1", 
       x = "Category", 
       y = "Percent of Cells") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), 
            position = position_stack(vjust = 0.5))


#Percent of cells expressing Thbs4---------
# Extract the expression data
# Method 1: Using AverageExpression with slot parameter
library(ggplot2)
library(dplyr)
library(broom)

# First, prepare the data with statistics
gene_name <- "Thbs4"
gene_data <- GetAssayData(DAA1, slot = "data")[gene_name, ]
cluster_ids <- DAA1$RNA_snn_res.0.25

# Create data frame for analysis
expr_data <- data.frame(
  Expression = as.numeric(gene_data),
  Cluster = cluster_ids,
  Cell_ID = names(gene_data)
)

#%>%
#  filter(Expression > 0)  # Remove zeros for analysis

# Calculate summary statistics per cluster
cluster_stats <- expr_data %>%
  group_by(Cluster) %>%
  summarise(
    Mean_Expression = mean(Expression),
    Median_Expression = median(Expression),
    SE = sd(Expression) / sqrt(n()),
    N_expressing = n(),
    .groups = 'drop'
  )

# Get total cells per cluster for percentage calculation
total_cells_per_cluster <- table(DAA1$RNA_snn_res.0.25)
cluster_stats$Total_cells <- total_cells_per_cluster[cluster_stats$Cluster]
cluster_stats$Percent_expressing <- (cluster_stats$N_expressing / cluster_stats$Total_cells) * 100

print(cluster_stats)

# Method 1: Pairwise t-tests between all clusters
pairwise_tests <- pairwise.t.test(expr_data$Expression, 
                                  expr_data$Cluster,
                                  p.adjust.method = "bonferroni")

print(pairwise_tests)

# Method 2: ANOVA with post-hoc tests
library(multcomp)

# Fit ANOVA model
anova_model <- aov(Expression ~ Cluster, data = expr_data)
summary(anova_model)

# Tukey's HSD for multiple comparisons
tukey_results <- TukeyHSD(anova_model)
print(tukey_results)

# Convert to data frame for plotting
tukey_df <- as.data.frame(tukey_results$Cluster)
tukey_df$comparison <- rownames(tukey_df)
tukey_df$significant <- tukey_df$`p adj` < 0.05

# Function to convert p-values to significance stars
get_significance_stars <- function(p_val) {
  case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01 ~ "**", 
    p_val < 0.05 ~ "*",
    p_val < 0.1 ~ ".",
    TRUE ~ ""
  )
}

# Add significance stars to cluster stats
cluster_stats$stars <- ""  # Initialize

# For demonstration, let's compare each cluster to cluster 0 (reference)
if("0" %in% cluster_stats$Cluster) {
  ref_cluster <- "0"
  
  for(i in 1:nrow(cluster_stats)) {
    if(cluster_stats$Cluster[i] != ref_cluster) {
      # Perform t-test against reference cluster
      cluster_data <- expr_data$Expression[expr_data$Cluster == cluster_stats$Cluster[i]]
      ref_data <- expr_data$Expression[expr_data$Cluster == ref_cluster]
      
      if(length(cluster_data) > 0 & length(ref_data) > 0) {
        t_test <- t.test(cluster_data, ref_data)
        cluster_stats$stars[i] <- get_significance_stars(t_test$p.value)
      }
    }
  }
  
}


# Plot 1: Bar plot with error bars and significance stars
p1 <- ggplot(cluster_stats, aes(x = Cluster, y = Mean_Expression)) +
  geom_bar(stat = "identity", fill = "lightblue", color = "black", alpha = 0.7) +
  geom_errorbar(aes(ymin = Mean_Expression - SE, ymax = Mean_Expression + SE),
                width = 0.2, color = "black") +
  geom_text(aes(label = stars, y = Mean_Expression + SE + 0.1),
            size = 6, color = "red") +
  theme_minimal() +
  labs(title = paste("Average Expression of", gene_name, "(excluding zeros)"),
       subtitle = "Error bars = SE, *p<0.05, **p<0.01, ***p<0.001 vs Cluster 0",
       x = "Cluster",
       y = "Mean Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

# Plot 2: Combined expression and percentage plot
max_expr <- max(cluster_stats$Mean_Expression + cluster_stats$SE)
scale_factor <- max_expr / 100  # Scale percentage to fit on same axis

p2 <- ggplot(cluster_stats, aes(x = Cluster)) +
  geom_bar(aes(y = Mean_Expression), stat = "identity", 
           fill = "lightblue", color = "black", alpha = 0.7) +
  geom_errorbar(aes(ymin = Mean_Expression - SE, ymax = Mean_Expression + SE),
                width = 0.2, color = "black") +
  geom_point(aes(y = Percent_expressing * scale_factor), 
             color = "red", size = 3) +
  geom_line(aes(y = Percent_expressing * scale_factor, group = 1), 
            color = "red", linetype = "dashed") +
  geom_text(aes(label = stars, y = Mean_Expression + SE + max_expr * 0.05),
            size = 6, color = "red") +
  scale_y_continuous(
    name = "Mean Expression",
    sec.axis = sec_axis(~ . / scale_factor, name = "% Expressing Cells")
  ) +
  theme_minimal() +
  labs(title = paste("Expression Analysis of", gene_name),
       subtitle = "Blue bars = Mean expression, Red line = % expressing cells",
       x = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y.right = element_text(color = "red"),
        axis.text.y.right = element_text(color = "red"))

print(p2)

library(ggpubr)

# Create violin plot with statistical comparisons
p3 <- ggviolin(expr_data, x = "Cluster", y = "Expression", 
               fill = "Cluster", palette = "jco",
               add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(method = "anova", label.y = max(expr_data$Expression) * 1.1) +
  stat_compare_means(comparisons = list(c("0", "1"), c("0", "2"), c("1", "2")),
                     method = "t.test", label = "p.signif") +
  labs(title = paste("Distribution of", gene_name, "Expression by Cluster"),
       subtitle = "Only cells with expression > 0 shown",
       x = "Cluster",
       y = "Expression Level") +
  theme(legend.position = "none")

print(p3)
