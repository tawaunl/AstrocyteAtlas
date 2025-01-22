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

# Plot D ----------
p1 <- plot_4way(res1 = human.AD,res2 = human.MS,"AD","MS")

p2 <- plot_4way(res1 = human.AD,res2 = human.PD,"AD","PD")
p3 <- plot_4way(res1 = human.PD,res2 = human.MS,"PD","MS")

cairo_pdf("~/Documents/AstrocytePaper/Supplemental/SuppFig4/Disease4ways.pdf",
          width=40,height = 20)
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
