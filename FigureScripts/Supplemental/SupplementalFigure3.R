#A. LPS Volcao plot -------------------------------------
library(ggplot2)
library(dplyr)
library(scales)
dir <- "~/Documents/AstrocytePaper/Supplemental/SuppFig3"
lps.de <- read.csv("~/Documents/AstrocytePaper/Figure3/LPSvsSalineDE.csv")
mycolors <- c(hue_pal()(2), "grey")
names(mycolors) <- c("Saline", "LPS", "unchanged")
colnames(lps.de)[1] <- "symbol"
# Define thresholds for labeling outliers
positive_logFC_thresh <- 2  # Threshold for upregulated genes
negative_logFC_thresh <- -1.5 # Threshold for downregulated genes
pvalue_thresh <- 0.05         # Common p-value threshold for both groups
tolabel <- c("C4b","Serpina3n","Gfap")
# Create a logical vector for labeling outlier points
lps.de$label <- with(lps.de, 
                     ((logFC > positive_logFC_thresh & PValue < pvalue_thresh) | 
                        (logFC < negative_logFC_thresh & PValue < pvalue_thresh) |
                        symbol %in% tolabel))
lfc_max <- max(abs(lps.de$logFC)) * 1.1  # Adding 10% for padding
x_limits <- c(-lfc_max, lfc_max)

cairo_pdf(file.path(dir,"LPSvsSaline_VolcanoPlot.pdf" ),
          width = 16,height=12,)
ggplot(data=lps.de, aes(x=logFC, y=-log10(FDR), col=diffexpressed, label=symbol)) + 
  ggrastr::geom_point_rast(alpha=0.6,size=3.5) + 
  theme_classic()+ geom_text_repel(size=6,data = subset(lps.de, label),
                                   show.legend = FALSE) +
  scale_colour_manual(values = mycolors) + ggtitle("LPS vs Saline")+
  theme(legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.key.size = unit(1 ,'cm'),
        axis.title = element_text(size=24),
        plot.title = element_text(size=33,hjust = 0.5)) +
  xlab("Avg. LogFoldChange") + ylab(bquote(-log[10](PValue))) +
  scale_x_continuous(limits = x_limits) 

dev.off()

#B UMAP scoring --------------------
## Bulk LPS ==================================
lps.genes <- c("CD44", "GALNTL2", "SULF2", "CRISPLD2", "DCN", "THBS2", "LRG1",
               "GGTA1", "AMIGO2", "FBL5", "ICAM1", "NFASC", "SORBS1", "CHI3L1",
               "CD14", "SRGN", "S1PR3", "TLR2", "OASL2", "MPA2L", "SAA3.ZC3HAV1",
               "TRIM30A", "LCN2", "TGTP1", "IRGM1", "IRGM1", "TNFAIP2", "PTX3",
               "TAPBP", "TAPBPL", "B2M", "H2-D1", "H2-K1", "H2-Q6", "H2-T10", 
               "H2-T23", "PMSB8", "PSMB9", "PROCR","STEAP4", "CP", "KCTD1",
               "SLC22A4", "SLC39A14", "SLC1A5", "SLC10A6","CXCL1", "CXCL10",
               "CXCL2", "SPP1", "IL1R1", "IL13RA", "OSMR","IFI202B", "IFI44", 
               "IFITM3", "IGTP", "IIGP1", "GBP2", "GBP3","BCL3", "CDKN1A",
               "PARP14", "XAF1", "NEK6","C1RB", "C1RA", "C3", "SERPING1", "C4B",
               "C1S","GFAP", "VIM", "TAGLN2", "SYNPO","CD109", "A2M", "SERPINA3N", 
               "TIMP1","FKBP5", "HSP1", "HSPB6","SPHKAP", "MAP3K6", "RHOJ","GAP43",
               "OLFM1", "SEMA4C","GPX3", "GSR",	"ACSL5", "HPGD", "ASPG", "CELA1",
               "USP18", "CPNE8", "LY6E", "UGT1A1", "ENDOU", "PLIN4", "TSPO", "TGM1",
               "NT5E", "ANGPT1", "TGM2", "S100A10", "TSPAN4", "SLC43A3", "IER3")
lps.genes <- str_to_title(lps.genes)

score_intersect <- intersect(lps.genes,rownames(data))
scores_by_cell <- colMeans(
  GetAssayData(data,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
data[["LPSAstroScore"]] <- scores_by_cell

data$ClusterNames <- factor(data$finalClusters)
levels(data$ClusterNames) <-c("Homeostatic","DAA1","DAA2","Synapse")

scvi_coords <- get_scvi_coords(data,data$finalClusters)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))
text <- data$ClusterNames
text_x <- vapply(split(scvi_coords$UMAP1, text), median, FUN.VALUE=0)
text_y <- vapply(split(scvi_coords$UMAP2, text), median, FUN.VALUE=0)

png(file.path(dir,"Figure3","C_LPSAstroScoringUMAP.png"),
    width = 1500,height=1500,res=300)
ggplot(scvi_coords%>%
         arrange(LPSAstroScore),aes(x=UMAP1, y=UMAP2, colour=LPSAstroScore)) +
  geom_point() +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) +
  scale_colour_gradientn(limits= c(0,1), colours = RColorBrewer::brewer.pal(5,"Purples"),
                         oob=squish)
dev.off()
##Mcao -------------------
mcao.genes <- c("CD24A", "CD44", "JUB", "ChL4", "COL6A1", "IGFBP3", "COL6A2", "COL12A1",
                "ICAM1", "LGALS1", "PVR", "THBS1", "THBS2", "TGFBI", "VCAN",
                "CHI3L1", "DPYSL3", "ECM1", "ADAM12", "ADAMTS4", "CELA1", "ADAMTS5",
                "GALNTL2", "LGALS3", "NAV2", "FGL2", "NETO2", "OLFML3","AKR1B8",
                "ASPG", "BCAT1", "ODC1", "ESD", "PDE3B", "ASNS", "CYP1B1", "B3GNT5",
                "GCNT2", "GGTA1", "OASL2", "CTPS", "UCK2","BCL3", "CDBPD", "KLF5",
                "KLF6", "TGIF1", "AHR", "FOSL1", "HMGA", "FOSL2", "SBNO2", "NEK6",
                "S100A6", "ZWINT", "CDT1", "CCND1", "CDK6", "CDKN1A", "GADD45A",
                "HMGA2", "EMP1","BDKRB2", "GADD45B", "SPHK1", "SOCS3", "CAMK2D",
                "RHOJ", "SPHKAP", "ODZ2", "SHISA6","CD14", "C3", "NUPR1",
                "PROCR", "THBD", "TLR4", "MPA2L", "PLP2", "PTX3", "LCN2", 
                "CLCF1", "CCL2", "CXCL1", "CXCL10", "CXCL2", "IL6", "LIF", "SPP1",
                "ARPC1B", "ACTN1", "FSCN1", "FLNC", "FLNA", "SYNPO", "MSN", "TAGLN2",
                "STEAP4", "CP", "STEAP1", "SLC10A6", "SLC39A14", "SLC5A3","AKAP12",
                "GFAP", "LMNA", "NES", "TUBA1A", "TUBB6", "VIM", "AGPAT9", "CAV1",
                "LASS6", "PLA2G4A", "CH25H", "HPGD", "PTGS2","CD109", "A2M",
                "SERPINA3N", "SERPINE1", "SERPING1", "TIMP1","ANXA1", "ANXA2",
                "ANXA7", "ANXA3", "S100A10", "S100A11","GPX1", "HMOX1", "SRXN1",
                "TXNRD1","GCH1", "MCL1", "LITAF", "EDA2R", "SULF1","IFI202B",
                "IFITM3", "IFI203", "IIGP1", "GBP3","PAPPA", "PRSS23", "USP18",
                "LONRF1","GAP43", "VGF", "CACNG5", "SYT4","NHP2", "NOP58", 
                "EIF1A", "PCBP3","IL13RA1", "IL6RA", "OSMR", "S1PR3", "MET", 
                "TNFRSF12A","BDNF", "CTGF", "GDF15","NT5E", "RNF125", "RNF19B",
                "OCIAD2", "MRPS6", "MTHFD2","CPNE8", "EII2", "FAM129B", "GRB10",
                "HSPB1", "HSPB6", "TMEM74", "KLHDC8A", "LRRC59", "LRG1", "TSPAN4",
                "MEST", "TM4SF1", "LRRFIP1", "STX11", "TGM1", "SPATA13", "IER3", 
                "SLC44A3", "SLC7A1", "AHNAK", "AHNAK2", "H19", "BTG3", "CHAC1")

mcao.genes <- str_to_title(mcao.genes)

score_intersect <- intersect(mcao.genes,rownames(data))
scores_by_cell <- colMeans(
  GetAssayData(data,assay = "RNA",
               layer = "data")[score_intersect,],na.rm = TRUE)
data[["MCAOAstroScore"]] <- scores_by_cell

data$ClusterNames <- factor(data$finalClusters)
levels(data$ClusterNames) <-c("Homeostatic","DAA1","DAA2","Synapse")

scvi_coords <- get_scvi_coords(data,data$finalClusters)
colnames(scvi_coords) <- make.unique(colnames(scvi_coords))

png(file.path(dir,"Figure3","C_MCAOAstroScoringUMAP.png"),
    width = 1500,height=1500,res=300)
ggplot(scvi_coords %>%
         arrange(MCAOAstroScore),aes(x=UMAP1, y=UMAP2, colour=MCAOAstroScore)) +
  geom_point() +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", legend.title =element_text(size=14),
        legend.key.size = unit(1.5, 'cm')) +
  scale_colour_gradientn(limits= c(0,1), colours =  RColorBrewer::brewer.pal(5,"Purples"),
                         oob=squish)
dev.off()


# C. Bar Plots of AStro types ------

igfbp5.dir <- file.path("~","Documents/InsituAnalysis/Igfbp5_results")
tg_animals_igf <- c("20_01","20_02","21_01","21_02")
p_cutoffs <- c(0,5,10,15,25,50) # plaque size cutoffs
breaks <- c(0,10,20,30,40,50,Inf) # distance bins
load(file.path(igfbp5.dir,"data.files.Rdata"))

astro_data <- lapply(astro_data, function(x){
  x$normch3 <- x$Subcellular..Channel.3..Num.spots.estimated / x$Nucleus..Area
  x$normch4 <- x$Subcellular..Channel.4..Num.spots.estimated / x$Nucleus..Area
  x$normastroType <- "Other Astros"
  x$normastroType[x$normch3 >= 1 &
                    x$normch4 <=2.5] <- "DAA1"
  x$normastroType[x$normch3 < 2 &
                    x$normch4 >= 2] <- "DAA2"
  x
})

data <- bind_rows(astro_data, .id = "column_label")
astro_counts <- data %>%
  group_by(column_label,Genotype, normastroType) %>%
  tally(name = "Count")

total_astrocytes_per_genotype <- astro_counts %>%
  group_by(Genotype) %>%
  summarize(total_astrocytes = sum(Count), .groups = 'drop')

# Join the total astrocytes per genotype to the original dataframe
df_with_totals <- astro_counts %>%
  left_join(total_astrocytes_per_genotype, by = "Genotype")

# Calculate the proportion of each astrocyte subtype to the total within each genotype
df_ratios <- df_with_totals %>%
  mutate(proportion_to_total_astrocytes = Count / total_astrocytes)

## Plot C ---------
color.codes <- c("goldenrod1","red","black")
zone <- levels(factor(astro_counts$normastroType))
cairo_pdf("~/Documents/AstrocytePaper/Supplemental/SuppFig3/AstrocyteCount.pdf",
          height = 10,width=10)
ggplot(astro_counts,aes(x=Genotype,y=Count)) + facet_wrap(~normastroType,scales = "free") +
  geom_boxplot(aes(fill=normastroType),width=0.5) + geom_point(aes(color=column_label)) +
  scale_fill_manual(values=setNames(color.codes, zone)) +theme_pubclean() +
  ggtitle("Normalized Astrocyte Count")+ 
  theme(legend.position="right",
        plot.title = element_text(hjust = 0.5,size=20),
        axis.title.x = element_blank())
dev.off()
cairo_pdf("~/Documents/AstrocytePaper/Supplemental/SuppFig3/AstrocyteCount_proportion.pdf",
          width = 5,height = 8)
ggplot(df_ratios,aes(fill = normastroType,y=proportion_to_total_astrocytes,x=Genotype))+
  geom_bar(position="stack", stat="identity",width=0.5) +
  scale_fill_manual(values=setNames(color.codes, zone)) +theme_pubclean() +
  ggtitle("Igfbp5")+theme(legend.position="right") 
dev.off()
## Gene expression ----------
astro <- bind_rows(astro_data, .id = "column_label")
df_long <- astro %>%
  pivot_longer(cols = c(normch3, normch4),
               names_to = "gene",
               values_to = "expression")
plot <- df_long %>%
  group_by(column_label,Genotype, gene) %>%
  summarize(mean_expression = mean(expression), .groups = 'drop')


cairo_pdf("~/Documents/AstrocytePaper/Supplemental/SuppFig3/GeneExpression.pdf")
ggplot(plot, aes(x =Genotype, y = mean_expression, fill = Genotype)) +
  geom_boxplot(width=0.5) + geom_point(size=2)+
facet_wrap(~ gene, scales = "free") +
  theme_minimal() +
  labs(title = "Boxplot of Gene Expressions with Averages",
       x = "Genotype",
       y = "Expression")
dev.off()

# D. S100a6 -------

s100.dir <- file.path("~","Documents/InsituAnalysis/S100a6_results")
tg_animals_s100 <- c("31_02", "31_03", "32_01", "32_02" ,"33_02")
load(file.path(s100.dir,"data.files.Rdata"))

astro_data <- lapply(astro_data, function(x){
  x$normch3 <- x$Subcellular..Channel.3..Num.spots.estimated / x$Nucleus..Area 
  x$normch4 <- x$Subcellular..Channel.4..Num.spots.estimated / x$Nucleus..Area
  x$normastroType <- "DN"
  x$normastroType[x$normch3 >= 1 &
                    x$normch4 <=2] <- "DAA1"
  x$normastroType[x$normch3 < 0.5 &
                    x$normch4 >= 2 ] <- "DAA2"
  x
})

data <- bind_rows(astro_data, .id = "column_label")
astro_counts <- data %>%
  group_by(column_label,Genotype, normastroType) %>%
  tally(name = "Count")
ggplot(astro_counts,aes(x=Genotype,y=Count)) + facet_wrap(~normastroType,scales = "free") +
  geom_boxplot(aes(fill=normastroType)) + geom_point(aes(color=column_label)) +
  scale_fill_manual(values=setNames(color.codes, zone)) +theme_pubclean() +
  ggtitle("Normalized Astrocyte Count")
