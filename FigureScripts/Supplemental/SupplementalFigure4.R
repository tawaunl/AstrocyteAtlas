#D.  4-ways ------
## Load libraries -----------------
library(ggplot2)
library(clusterProfiler)
library(rstatix)
library(ggrepel)
library(UpSetR)
library(dplyr)
library(tidyverse)
source("~/Documents/scHelpers.R")
dir <- "~/Documents/AstrocytePaper/Human"
human.AD <- readRDS(file.path(dir,"res_dsl.AD.cain_fixed.072924.rds"))

human.MS <- readRDS(file.path(dir,"res_dsl.MS.080124.rds"))
human.PD <- readRDS(file.path(dir,"res_dsl.PD.072924.rds"))

human.AD <- generatePlotTable(human.AD)
human.MS <- generatePlotTable(human.MS)
human.PD <- generatePlotTable(human.PD)

feat <- readRDS(file.path(dir,"humanFeatures.rds"))
human.AD$diffexpressed <- ifelse(human.AD$FDR<=0.05 & human.AD$dl_mu>=0.5 & human.AD$n_up>=3, "up",
                                 ifelse(human.AD$FDR<=0.05 & human.AD$dl_mu<= -0.5 & human.AD$n_down>=3, "down","unchanged"))
human.MS$diffexpressed <- ifelse(human.MS$FDR<=0.05 & human.MS$dl_mu>=0.5 & human.MS$n_up>=3, "up",
                                 ifelse(human.MS$FDR<=0.05 & human.MS$dl_mu<= -0.5 & human.MS$n_down>=3, "down","unchanged"))
human.PD$diffexpressed <- ifelse(human.PD$FDR<=0.05 & human.PD$dl_mu>=0.5 & human.PD$n_up>=2, "up",
                                 ifelse(human.PD$FDR<=0.05 & human.PD$dl_mu<= -0.5 & human.PD$n_down>=2, "down","unchanged"))
protein_coding <- feat %>% dplyr::filter(type == "protein_coding")
human.AD <- human.AD %>% dplyr::filter(symbol %in% protein_coding$symbol)
human.MS <- human.MS %>% dplyr::filter(symbol %in% protein_coding$symbol)
human.PD <- human.PD %>% dplyr::filter(symbol %in% protein_coding$symbol)
plot_4way <- function(res1, res2, disease1, disease2, lfcCutoff=0.5,PvalCutoff=0.05) {
  require(ggrepel)
  require(tidyr)
  
  res1 <- as.data.frame(res1) %>% dplyr::mutate(diffexpressed=ifelse(dl_mu>=lfcCutoff & PValue<=PvalCutoff , "up", 
                                                                     ifelse(dl_mu<=-lfcCutoff & PValue<=0.05 , "down", "unchanged"))) 
  
  res2 <- as.data.frame(res2) %>% dplyr::mutate(diffexpressed=ifelse(dl_mu>=lfcCutoff & PValue<=PvalCutoff , "up", 
                                                                     ifelse(dl_mu<=-lfcCutoff & PValue<=PvalCutoff , "down", "unchanged")))                                    
  res_merged <- merge(res1,res2,by="symbol")
  
  res_merged <- res_merged %>% dplyr::mutate(sig=
                                               ifelse(diffexpressed.x %in% c("up", "down") & diffexpressed.y %in% c("up", "down"), "Significant in both", 
                                                      ifelse(diffexpressed.x %in% c("up", "down") & diffexpressed.y == "unchanged", paste("Changed only in ", disease1, sep=""), 
                                                             ifelse(diffexpressed.y %in% c("up", "down") & diffexpressed.x == "unchanged", paste("Changed only in ", disease2, sep=""), "Unchanged in both"))))
  
  res_merged <- res_merged  %>% dplyr::select(dl_mu.x,dl_mu.y,diffexpressed.x, diffexpressed.y,sig,symbol) %>% unique %>% drop_na()
  
  p <- ggplot(res_merged, aes(x=dl_mu.x, y = dl_mu.y, color=sig)) +
    ggrastr::geom_point_rast(size=3,alpha=0.6) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 30),
          legend.text = element_text(size=24)) +
    scale_color_manual(values=c("goldenrod1",  "red", "mediumblue", "gray67"), name = "") + geom_hline(yintercept =0) + 
    geom_vline(xintercept=0) + 
    geom_text_repel(data = res_merged[res_merged$sig == "Significant in both",],size=9, aes(x=dl_mu.x,y=dl_mu.y,label=symbol),show.legend = FALSE,box.padding = 0.8, max.overlaps = 20, fontface = "italic") +
    xlab(paste("LogFC in ", disease1, sep="")) + ylab(paste("LogFC in ", disease2, sep="")) + guides(fill=guide_legend(title=NULL))
  
  list_out <- list()
  list_out[[1]] <- res_merged
  list_out[[2]] <- p
  return(list_out)
}

# Plot D ----------
p1 <- plot_4way(res1 = human.AD,res2 = human.MS,"AD","MS")

p2 <- plot_4way(res1 = human.AD,res2 = human.PD,"AD","PD")
p3 <- plot_4way(res1 = human.PD,res2 = human.MS,"PD","MS")

cairo_pdf("~/Documents/AstrocytePaper/Supplemental/SuppFig4/Disease4ways.pdf",
          width=45,height = 10)
cowplot::plot_grid(p1[[2]],p2[[2]],p3[[2]],ncol = 3)
dev.off()



# E. Upset -----------
##DE ----------------
AD <- human.AD %>% filter(diffexpressed=="up" )
MS <- human.MS %>% filter(diffexpressed=="up" ) 
PD <- human.PD %>% filter(diffexpressed=="up" ) 

plotData <- list(AD=AD$symbol,MS=MS$symbol,PD=PD$symbol)

plotDE <- upset(fromList(plotData), order.by = "freq",
            matrix.color = "black",main.bar.color = "black",sets.bar.color = "black",
            mainbar.y.label = "Gene Intersections",sets.x.label = "# Enriched Genes",
            text.scale = 2,point.size = 5,line.size = 1.5)
plot

## Pathways -------------

DiseasePaths <- readRDS(file.path(dir,"HumanDiseasePathwaysGO.rds"))

AD <- DiseasePaths$AD_GSEAGO@result %>% filter(NES > 0  & pvalue < 0.05)
MS <- DiseasePaths$MS_GSEAGO@result %>% filter(NES > 0  & pvalue < 0.05) 
PD <- DiseasePaths$PD_GSEAGO@result %>% filter(NES > 0  & pvalue < 0.05) 

plotData <- list(AD=AD$ID,MS=MS$ID,PD=PD$ID)

plotPath <- upset(fromList(plotData), order.by = "freq",
                matrix.color = "black",main.bar.color = "black",sets.bar.color = "black",
                mainbar.y.label = "Pathway Intersections",sets.x.label = "# Enriched Pathways",
                text.scale = 2,point.size = 5,line.size = 1.5)

pdf("~/Documents/AstrocytePaper/Supplemental/SuppFig4/DEupsetplot.pdf")
print(plotDE)
dev.off()

pdf("~/Documents/AstrocytePaper/Supplemental/SuppFig4/PaathwaysUpset.pdf")
print(plotPath)
dev.off()

# F. Cluster 2 markers
