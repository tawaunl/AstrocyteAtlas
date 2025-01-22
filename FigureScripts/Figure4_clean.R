#Fig 4
pkgs <-c("ggplot2","ggpubr","viridis","gridExtra",
         "dplyr","tidyverse","Seurat","SingleCellExperiment")
pacman::p_load(pkgs,character.only = TRUE)
dir <- "~/Documents/AstrocytePaper"
data <- readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))

# A. DotPlot ===================
data_3xtg <- subset(data,subset= StudyName=="Lee TauPS2APP")

data_3xtg <- subset(data_3xtg,subset= Model=="Control" | Model=="TauPS2APP")

## Plot A -----------------------------------
pdf(file.path(dir,"Figure4","Dotplot.pdf"))
DotPlot(data_3xtg,features=c("Apoe","S100a6","Igfbp5"),col.min = -1,col.max = 1,
        cols = "RdYlBu",group.by = "finalClusters",scale = TRUE) + RotatedAxis() +coord_flip() + xlab('Gene') +  ylab('Cluster') +
  theme(axis.text = element_text(size = 14),axis.title = element_text(size=18,face='bold')) 
dev.off()

# C. Scatter Plot of expression ----------------------
igfbp5.dir <- file.path("~","Documents/InsituAnalysis/Igfbp5_results")
tg_animals_igf <- c("20_01","20_02","21_01","21_02")
p_cutoffs <- c(0,5,10,15,25,50) # plaque size cutoffs
breaks <- c(0,10,20,30,40,50,Inf) # distance bins
load(file.path(igfbp5.dir,"data.files.Rdata"))

astro_data <- lapply(astro_data, function(x){
  x$normch3 <- x$Subcellular..Channel.3..Num.spots.estimated / x$Nucleus..Area
  x$normch4 <- x$Subcellular..Channel.4..Num.spots.estimated / x$Nucleus..Area
  x$normastroType <- "DN"
  x$normastroType[x$normch3 >= 1 &
                    x$normch4 <=2.5] <- "DAA1"
  x$normastroType[x$normch3 < 1 &
                    x$normch4 >= 2] <- "DAA2"
  x
})

data <- bind_rows(astro_data, .id = "column_label")

astro_counts <- data %>%
  group_by(column_label,Genotype, normastroType) %>%
  tally(name = "Count")


ggplot(astro_counts,aes(x=Genotype,y=Count)) + facet_wrap(~normastroType,scales = "free") +
  geom_boxplot() + geom_point()




### Plot C ---------------------------------
png(file.path(dir,"Figure4","Igfbp5_ScatterSubtype.png"),
    width =2400,height = 2000,res=300)

ggscatter(data, x = "normch3",
          y = "normch4",color="normastroType",
          add = "reg.line", conf.int = FALSE,
          palette =c("goldenrod1","red","black"),
          cor.coef = F, cor.method = "pearson",
          xlab = "Estimated Igfbp5 spots",
          ylab = "Estimated Apoe spots",size = 4,alpha=0.7,cor.coef.size = 10)
dev.off()

# D. Boxplot of Distance ----------------

binnedData <- lapply(astro_data[tg_animals_igf],function(x){
  x$dis2nearestPlaque <- x[paste0("dis2nearestPlaque_",0)]
  test <- x %>%
    group_by(dis2nearestPlaque, normastroType) %>%
    summarise(Count = n())
  colnames(test)[1] <- "dis2nearestPlaque"
  aggregate(test[,-c(1,2)], by = list(distance=cut(unlist(as.vector(test$dis2nearestPlaque)),
                                                   breaks,include.lowest=T),
                                      Subtype = test$normastroType ), FUN = sum )
})


x <- list_rbind(x=binnedData, 
                names_to = "id") %>% arrange(c(distance)) %>% arrange(c(Subtype)) 

idv <- vector()
for(type in levels(factor(x$Subtype))){
  for(dis in levels(factor(x$distance))){
    idv = c(idv,(x$Count[x$distance== dis & x$Subtype == type]/sum(x$Count[x$distance== dis])))
  }
}

df <- data.frame(Distance = (x$distance),
                 Subtype = x$Subtype,
                 Sample=factor(x$id),
                 ratios = idv)

color.codes <- c("goldenrod1","red","black")
zone <- levels(factor(df$Subtype))
### Plot D -----------------------------------------------
pdf(file.path(dir,"Figure4","Igfbp5_DistancePlot.pdf"),
    width=10,height = 6)
ggplot(df, aes(x=Distance,y=ratios,fill=Subtype)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(group=Subtype),size=2,alpha=1, color="darkgrey",
             position = position_dodge(width = .78)) +
  scale_fill_manual(values=setNames(color.codes, zone)) +
  theme_pubclean()+
  ylab("Ratio of Astrocyte type") +
  xlab(paste("Distance to nearest plaque "))+
  theme(axis.title = element_text(size=24, face = "bold"),
        axis.text = element_text(size=14),
        legend.position = "right")
dev.off()

# E Stacked Bar plot of distance--------
color.codes <- c("goldenrod1","red","black")
zone <- levels(factor(df$Subtype))
png(file.path(dir,"Figure4","Igfbp5_DistancePlot_Stacked.png"),
    res=300,width=1500,height=2000)
ggplot(df,aes(fill = Subtype,y=ratios,x=Distance))+
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=setNames(color.codes, zone)) +theme_pubclean() +
  ggtitle("Igfbp5")+theme(legend.position="right") 
dev.off()
#Statistics ---------------
library(rstatix)
df %>%
  group_by(Subtype, Distance) %>%
  identify_outliers(ratios) # there are outliers 

df %>%
  group_by(Subtype, Distance) %>%
  shapiro_test(ratios) # normailty is present at most levels aside from 2


ggqqplot(df, "ratios", ggtheme = theme_bw()) +
  facet_grid(Distance ~ Subtype, labeller = "label_both") 
# plot shows points fall approximately along the reference line, we can assume normality.

## Test with outliers -----------------

res.aov <- anova_test(
  data = df, dv = ratios, wid = Sample,
  within = c(Subtype, Distance)
)
get_anova_table(res.aov)  
# There is a statistically significant two-way interactions between Subtype and Distance
# F(6,18) = 4.54 p<0.000001

### Post hoc ======================
# effect of Subtype at each distance 
one.way <- df %>%
  group_by(Distance) %>%
  anova_test(dv = ratios, wid = Sample, within = Subtype) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

pwc <- df %>%
  group_by(Distance) %>%
  pairwise_t_test(
    ratios ~ Subtype, paired = TRUE,
    p.adjust.method = "bonf"
  )
pwc

pwc <- pwc %>% add_xy_position(x = "Distance")
ggboxplot(
  df, x = "Distance", y = "ratios",
  color = "Subtype", palette = "jco") + 
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

write.csv(df , "~/Documents/AstrocytePaper/2024-11-26 lucast3/Figure4D_DistanceDataforStatistics.csv")
