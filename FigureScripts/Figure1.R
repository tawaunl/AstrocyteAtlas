# Script to generate all plots in Figure 1

# ======================== Load Libraries ======================================
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
library(cowplot)
library(ggrepel)
source("~/Documents/scHelpers.R")
# Load Data ====================================================================
dir <- "~/Documents/AstrocytePaper"
data<- readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))

# B. AD and MS Volcano Plots ===================================================
AD.DE <- readRDS(file.path(dir, "res_dsl_AD_full.011824.rds"))
MS.DE <- readRDS(file.path(dir, "res_dsl_MS_full.011824.rds"))

##  AD Volcano ================================================================
AD.DE_plot <- generatePlotTable(AD.DE)
mycolors <- c("blue", "firebrick3","grey")
names(mycolors) <- c("down", "up", "unchanged")
AD.DE_plot$diffexpressed <- "unchanged"

AD.DE_plot$diffexpressed[AD.DE_plot$dl_mu >= 0.25 &
                           AD.DE_plot$FDR<=0.05 &
                           AD.DE_plot$sig=="yes" ] <- "up"

AD.DE_plot$diffexpressed[AD.DE_plot$dl_mu <= -0.25 &
                             AD.DE_plot$FDR<=0.05 &
                             AD.DE_plot$sig=="yes" ] <- "down"

# Define specific genes to force label
specific_genes <- c("C4b", "Vim", "Serpina3n", "Gfap")
extrafont::font_import() # Uncomment and run this once to import fonts
extrafont::loadfonts(device = "pdf") # Ensure fonts are loaded

# Example thresholds
pos_logFC_thresh <- 4.5  # Positive threshold
neg_logFC_thresh <- -1 # Negative threshold
pvalue_thresh <- 1e-17    # P-value threshold

# Identify genes based on criteria for labeling
AD.DE_plot$label <- with(AD.DE_plot, (dl_mu > pos_logFC_thresh | (dl_mu < neg_logFC_thresh & PValue < 0.05) | PValue < pvalue_thresh) | symbol %in% specific_genes)

AD.DE_plot$label[which(AD.DE_plot$label==TRUE &AD.DE_plot$diffexpressed=="unchanged")] <- FALSE

p1 <- ggplot(data=AD.DE_plot, aes(x=dl_mu, y=-log10(PValue), col=diffexpressed, label=symbol)) + 
  ggrastr::geom_point_rast(alpha=0.6,size=3.5) + 
  theme_classic() +
  geom_text_repel(data = subset(AD.DE_plot,label),size=7, fontface = "italic",
                  aes(label = symbol),show.legend = FALSE,box.padding = 0.5,
                  force=3, max.overlaps = 15)+ 
  scale_colour_manual(values = mycolors)+
  theme(legend.position = "none",
        axis.title = element_text(size=22),
        plot.title = element_text(size=22,face='bold',hjust = 0.5)) +
  xlab("Meta-LogFoldChange") + ylab(bquote(bold(-log[10](Meta-PValue)))) +
  ylim(0,28) +xlim(-5,8)


p1

cairo_pdf(file.path(dir,"Figure1","A_ADVolcanoPlot.pdf"),
    width = 11.5,height=10)
print(p1)
dev.off()
## MS Volcano =================================================================
MS.DE_plot <- generatePlotTable(MS.DE)

MS.DE_plot <- MS.DE_plot %>% mutate(sig = ifelse((FDR<=0.05 | FDR <= 0.05) & n_up | n_down >= n_tested/2 & abs(dl_mu) >= 0.5, "yes", "no"))


MS.DE_plot$diffexpressed <- "unchanged"

MS.DE_plot$diffexpressed[MS.DE_plot$dl_mu >= 0.25 &
                             MS.DE_plot$FDR<=0.05 &
                             MS.DE_plot$sig=="yes" ] <- "up"

MS.DE_plot$diffexpressed[MS.DE_plot$dl_mu <= -0.25 &
                           MS.DE_plot$FDR<=0.05 &
                           MS.DE_plot$sig=="yes" ] <- "down"

DiseaseDE <- list(AD=AD.DE_plot,MS=MS.DE_plot)


# Example thresholds
pos_logFC_thresh <- 1  # Positive threshold
neg_logFC_thresh <- -0.3 # Negative threshold
pvalue_thresh <- 1e-7    # P-value threshold

# Identify genes based on criteria for labeling
MS.DE_plot$label <- with(MS.DE_plot, (dl_mu > pos_logFC_thresh | (dl_mu < neg_logFC_thresh & PValue < 0.05) | PValue < pvalue_thresh) | symbol %in% specific_genes)

MS.DE_plot$label[which(MS.DE_plot$label==TRUE &MS.DE_plot$diffexpressed=="unchanged")] <- FALSE

p2 <- ggplot(data=MS.DE_plot, aes(x=dl_mu, y=-log10(PValue), col=diffexpressed, label=symbol)) + 
  ggrastr::geom_point_rast(alpha=0.6,size=3.5) + 
  theme_classic() +
  geom_text_repel(data = subset(MS.DE_plot,label),size=7, fontface = "italic",
                  aes(label = symbol),show.legend = FALSE,box.padding = 0.5,
                  force=3, max.overlaps = 15)+ 
  scale_colour_manual(values = mycolors)+
  theme(legend.position = "none",
        axis.title = element_text(size=22),
        plot.title = element_text(size=22,face='bold',hjust = 0.5)) +
  xlab("Meta-LogFoldChange") + ylab(bquote(bold(-log[10](Meta-PValue)))) +
  ylim(0,11) +xlim(-5,8)

p2


cairo_pdf(file.path(dir,"Figure1","A_MSVolcanoPlots.pdf"),
    width = 11.5,height=10)
print(p2)
dev.off()
# C. Pathway Analysis for AD and MS DE =========================================

library(clusterProfiler)

AD.DE_plot <- AD.DE_plot[order(AD.DE_plot$dl_mu,decreasing = TRUE),]
AD.genes <- AD.DE_plot$dl_mu
names(AD.genes) <- AD.DE_plot$symbol

gene.df <- bitr( names(AD.genes), fromType = "SYMBOL",
                 toType = c("ENSEMBL", "ENTREZID"),
                 OrgDb = org.Mm.eg.db::org.Mm.eg.db)
x <- match(gene.df$SYMBOL,names(AD.genes))
AD.gsea <- AD.genes[x]
names(AD.gsea) <- gene.df$ENTREZID


MS.DE_plot <- MS.DE_plot[order(MS.DE_plot$dl_mu,decreasing = TRUE),]
MS.genes <- MS.DE_plot$dl_mu
names(MS.genes) <- MS.DE_plot$symbol

gene.df <- bitr( names(MS.genes), fromType = "SYMBOL",
                 toType = c("ENSEMBL", "ENTREZID"),
                 OrgDb = org.Mm.eg.db::org.Mm.eg.db)
x <- match(gene.df$SYMBOL,names(MS.genes))
MS.gsea <- MS.genes[x]
names(MS.gsea) <- gene.df$ENTREZID
## Individual Pathways============

### GO:BP =====================
ego3 <- gseGO(geneList     = AD.gsea,
              OrgDb        = org.Mm.eg.db::org.Mm.eg.db,
              ont          = "BP",
              minGSSize    = 50,
              maxGSSize    = 500,
              eps          = 0,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

AD.paths <- dotplot(ego3,x="GeneRatio",color="pvalue")+ ggtitle("AD GO:BP Pathways") + 
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=20))

ego4 <- gseGO(geneList     = MS.gsea,
              OrgDb        = org.Mm.eg.db::org.Mm.eg.db,
              ont          = "BP",
              minGSSize    = 50,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              eps          = 0,
              verbose      = FALSE)

DiseasePaths <- list(AD_GSEAGO=ego3@result,MS_GSEAGO =ego4@result)
saveRDS(DiseasePaths,file.path(dir,"Figure1","DiseasePathways.rds"))

library(cluster)
## Grouped Pathways =================================
DiseasePaths <- readRDS(file.path(dir,"Figure1","DiseasePathways.rds"))
ck <- compareCluster(list("AD Models"=AD.gsea,"MS Models"=MS.gsea),
                     fun = "gseGO",OrgDb = org.Mm.eg.db::org.Mm.eg.db ,ont="BP")

ck1 <- filter(ck,NES>1)
ck2 <- filter(ck,NES<-1)

go_up <- dotplot(ck1, by="NES",showCategory = 5) + ggtitle("GO:BP Pathways Upregulated") + 
  theme(axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=20)) +xlab("Disease")
go_down <- dotplot(ck2, by="NES",showCategory = 5) + ggtitle("GO:BP Pathways Downregulated") + 
  theme(axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=18)) +
  theme(plot.title = element_text(hjust = 0.5,size=20))+xlab(" ")

png(file.path(dir,"Figure1","B_DotPlotPathways_ADvsMS.png"),
    width = 3500,height=2000,res = 300)
plot_grid(go_up,go_down)
dev.off()


cairo_pdf(file.path(dir,"Figure1","B_DotPlotPathways_ADvsMS.pdf"),
    width = 15,height=7)
plot_grid(go_up,go_down)
dev.off()
