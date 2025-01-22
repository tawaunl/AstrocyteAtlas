# C. Bar Plots of AStro types ------

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


## Plot C ---------
color.codes <- c("goldenrod1","red","black")
zone <- levels(factor(astro_counts$normastroType))
cairo_pdf("~/Documents/AstrocytePaper/Supplemental/SuppFig3/AstrocyteCount.pdf",
          height = 10,width=15)
ggplot(astro_counts,aes(x=Genotype,y=Count)) + facet_wrap(~normastroType,scales = "free") +
  geom_boxplot(aes(fill=normastroType)) + geom_point(aes(color=column_label)) +
  scale_fill_manual(values=setNames(color.codes, zone)) +theme_pubclean() +
  ggtitle("Normalized Astrocyte Count")+ 
  theme(legend.position="right",
        plot.title = element_text(hjust = 0.5,size=20),
        axis.title.x = element_blank())
dev.off()
