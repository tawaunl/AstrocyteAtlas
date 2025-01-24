dir <- "~/Documents/AstrocytePaper"
sav.dir <- file.path(dir,"Supplemental","SuppFig1")
#A.  4-way plot for Mouse AD and MS --------------------------------------------
library(gg4way)
library(ggplot2)
library(dplyr)
source("~/Documents/scHelpers.R")
AD.DE <- readRDS(file.path(dir, "res_dsl_AD_full.011824.rds"))
MS.DE <- readRDS(file.path(dir, "res_dsl_MS_full.011824.rds"))

AD.DE_plot <- generatePlotTable(AD.DE)
MS.DE_plot <- generatePlotTable(MS.DE)

x <- list("ADvsControl"=AD.DE_plot,"MSvsControl"=MS.DE_plot)

## Plot A ------

cairo_pdf(file.path(sav.dir,"Mouse_ADvsMS_4wayplot.pdf"),
    width=14,height = 12)
plot_4way(AD.DE_plot,MS.DE_plot,"AD models","MS models",0.25)
dev.off()

# B. Pathway Analysis -------------------

AD_plot <- gg4way(DGEdata = x,
                  x = "ADvsControl", FDR = "PValue",
                  y = "MSvsControl", sep = "vs",
                  logFC = "dl_mu",logFCcutoff = 0.25, FDRcutoff = 0.05,
                  label=T,
                  textSize =8,
                  colorVector = c("grey80", "goldenrod1", "red", "mediumblue")) +
  ggplot2::theme(legend.title = ggplot2::element_text(size = 14),
                 legend.text = ggplot2::element_text(size = 14),
                 axis.title = element_text(size=14)) +ggrastr::rasterise(geom_point(size=2))

library(clusterProfiler)
ff <- AD_plot |> getTotals()

data1 <- AD_plot$data
shared_up <- data1 %>% filter(`ADvsControl Direction` =="Up" &
                                `MSvsControl Direction`=="Up" &
                                Significant=="Significant in Both")
AD_up <- dat1a %>% filter(`ADvsControl Direction`=="Up" &
                            Significant == "Significant in ADvsControl")
MS_up <- data1 %>% filter(`MSvsControl Direction`=="Up" &
                            Significant == "Significant in MSvsControl")
pathData <- list("AD Distinct"=AD_up,"Shared in Both"=shared_up,"MS Distinct"=MS_up)

enrichPaths <- lapply(pathData,function(x){
  enrichGO(x$symbol,OrgDb = org.Mm.eg.db::org.Mm.eg.db,keyType = "SYMBOL",ont="BP",readable = T)
})

enrichPlots <- list()
for(set in 1:length(enrichPaths)){
  name <- names(enrichPaths)[set]
  enrichPlots[[name]] <- dotplot(enrichPaths[[set]]) +
    ggtitle(paste("GO:BP",name,"Pathways")) +
    theme(axis.title.x = element_text(size = 18,face="bold"),
          axis.text.x = element_text(angle = 90,size=14),
          axis.text.y = element_text(size=14)) +
    theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
}
topPaths <- lapply(enrichPaths,function(x){
  data<- x@result
  data <- data[order(data$p.adjust),]
})
names(topPaths) <- names(enrichPaths)
df <- bind_rows(topPaths, .id = "Cluster")

df <- df %>% 
  separate(GeneRatio,c("count","total"),sep="/") %>% 
  mutate(GeneRatio = as.double(count)/as.double(total))

## Function to filter first N GO term under specified p.adjust cutoff.
filterCategory <- function(df, N, cutoff){
  df_topn <- NULL
  for(i in levels(factor(df$Cluster))){
    df_slice <- df %>% filter(Cluster == i) %>% 
      arrange(p.adjust) %>% 
      filter(p.adjust< cutoff) %>% 
      head(N)
    df_topn <- bind_rows(df_topn, df_slice)
  }
  return(df_topn)
}
df_top10 <- filterCategory(df,5,0.05)      # Filter the results to show only top 10 GO terms less than p.adjust 0.05 in all cluster.
descriptions <- df_top10$ID

full <- df %>% filter(ID %in% descriptions)

## Order the GO term in a style similar to dotplot.

full1 <- filterCategory(full,10,0.05)
full1$Name <- str_wrap(full1$Description, width = 30) # wrap axis labels
full1 <- full1 %>% filter(Cluster != "Shared in Both") # not enough shared to genes to be meaningful

full1$Cluster <- factor(full1$Cluster,levels = c("MS Distinct",
                                                 "AD Distinct"))
full1$Name <- factor(full1$Name,
                     levels = c(full1$Name[3],full1$Name[1],full1$Name[2],
                                full1$Name[4],full1$Name[5],full1$Name[8],
                                full1$Name[9],full1$Name[10],full1$Name[6],
                                full1$Name[7]))
## Plot B -------------
pdf(file.path(sav.dir,"DistinctDEPathways_AD_MS.pdf"),width=7,height=7)
ggplot(full1, aes(x=Cluster, y=Name)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_colour_gradient(limits=c(min(full1$p.adjust), max(full1$p.adjust)), low="red",high="blue") +
  theme_bw() + ggtitle("MS/AD Distinct Pathways") +
  theme(axis.title = element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90,size=14),
        axis.text.y = element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))+ 
  scale_size(range = c(3,11))
dev.off()



write.csv(data,file.path,"AD_MS_sharedPathways")

#C.  Mouse Umaps -------------------------------------
## by study --------------------------------
data<- readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))

library(gridExtra)
studies <- unique(data$StudyName)
list_of_plots <- list()
for (i in 1:length(studies)) {
  p <- plot_umap_library(data,"StudyName",studies[i],color = "red")
  list_of_plots[[i]] <- p
}
### Plot C -----------
cairo_pdf(file.path(sav.dir ,"Mouse_UMAPS_byStudy.pdf"),
    width = 21, height = 15)
do.call("grid.arrange", c(list_of_plots, ncol=5))
dev.off()

## by Region --------------------------------
region <- unique(data$BrainRegion)
list_of_plots <- list()
for (i in 1:length(region)) {
  p <- plot_umap_library(data,"BrainRegion",region[i],color = "darkgreen")
  list_of_plots[[i]] <- p
}

### Plot D -------------------
cairo_pdf(file.path(sav.dir ,"Mouse_UMAPS_byRegion.pdf"),
    width = 13, height = 10)
do.call("grid.arrange", c(list_of_plots, ncol=3))
dev.off()


