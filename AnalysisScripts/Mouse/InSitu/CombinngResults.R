# Combining results
# Load packages ===========================================
pkgs <-c("ggplot2","ggpubr","viridis","gridExtra",
         "dplyr","tidyverse")
pacman::p_load(pkgs,character.only = TRUE)

igfbp5.dir <- file.path("~","Documents/InsituAnalysis/Igfbp5_results")
s100.dir <- file.path("~","Documents/InsituAnalysis/S100a6_results")
tg_animals_igf <- c("20_01","20_02","21_01","21_02")
tg_animals_s100 <- c("31_02", "31_03", "32_01", "32_02" ,"33_02")
p_cutoffs <- c(0,5,10,15,25,50) # plaque size cutoffs
breaks <- c(0,10,20,50,Inf) # distance bins

# load Igfbp5 data
load(file.path(igfbp5.dir,"data.files.Rdata"))
astro_data <- lapply(astro_data, function(x){
  x$normch3 <- x$Subcellular..Channel.3..Num.spots.estimated / x$Nucleus..Perimeter 
  x$normch4 <- x$Subcellular..Channel.4..Num.spots.estimated / x$Nucleus..Perimeter
  x$normastroType <- "DN"
  x$normastroType[x$normch3 >= 1 &
                    x$normch4 <=1.5] <- "DAA1"
  x$normastroType[x$normch3 < 0.5 &
                    x$normch4 >= 2 ] <- "DAA2"
  x
})

igf <- astro_data

# load S100a6 data 

load(file.path(s100.dir,"data.files.Rdata"))
astro_data <- lapply(astro_data, function(x){
  x$normch3 <- x$Subcellular..Channel.3..Num.spots.estimated / x$Nucleus..Perimeter 
  x$normch4 <- x$Subcellular..Channel.4..Num.spots.estimated / x$Nucleus..Perimeter
  x$normastroType <- "DN"
  x$normastroType[x$normch3 >= 1 &
                    x$normch4 <=2.5] <- "DAA1"
  x$normastroType[x$normch3 < 0.5 &
                    x$normch4 >= 1.5 ] <- "DAA2"
  x
})

s100 <- astro_data


# Run Analysis ============================================

StandardLabelPlots <- list()

for (cutoff in p_cutoffs) {
  binnedData_igf <- lapply(igf[tg_animals_igf],function(x){
    x$dis2nearestPlaque <- x[paste0("dis2nearestPlaque_",cutoff)]
    test <- x %>%
      group_by(dis2nearestPlaque, normastroType) %>%
      summarise(Count = n())
    colnames(test)[1] <- "dis2nearestPlaque"
    aggregate(test[,-c(1,2)], by = list(distance=cut(unlist(as.vector(test$dis2nearestPlaque)),
                                                     breaks, include.lowest=T),
                                        Subtype = test$normastroType ), FUN = sum )
  })
  binnedData_s100 <- lapply(s100[tg_animals_s100],function(x){
    x$dis2nearestPlaque <- x[paste0("dis2nearestPlaque_",cutoff)]
    test <- x %>%
      group_by(dis2nearestPlaque, normastroType) %>%
      summarise(Count = n())
    colnames(test)[1] <- "dis2nearestPlaque"
    aggregate(test[,-c(1,2)], by = list(distance=cut(unlist(as.vector(test$dis2nearestPlaque)),
                                                     breaks, include.lowest=T),
                                        Subtype = test$normastroType ), FUN = sum )
  })
  
  
  x <- list_rbind(x=binnedData_igf, 
                  names_to = "id") %>% arrange(c(distance)) %>% arrange(c(Subtype)) 
  y <- list_rbind(x=binnedData_s100, 
                  names_to = "id") %>% arrange(c(distance)) %>% arrange(c(Subtype))
  
  summary_igf <- c((x$Count[x$distance=="[0,10]" & x$Subtype =="DAA1"]/sum(x$Count[x$Subtype=="DAA1"])),
           (x$Count[x$distance=="(10,20]" & x$Subtype =="DAA1"]/sum(x$Count[x$Subtype=="DAA1"])),
           (x$Count[x$distance=="(20,50]" & x$Subtype =="DAA1"]/sum(x$Count[x$Subtype=="DAA1"])),
           (x$Count[x$distance=="(50,Inf]" & x$Subtype =="DAA1"]/sum(x$Count[x$Subtype=="DAA1"])),
           (x$Count[x$distance=="[0,10]" & x$Subtype =="DAA2"]/sum(x$Count[x$Subtype=="DAA2"])),
           (x$Count[x$distance=="(10,20]" & x$Subtype =="DAA2"]/sum(x$Count[x$Subtype=="DAA2"])),
           (x$Count[x$distance=="(20,50]" & x$Subtype =="DAA2"]/sum(x$Count[x$Subtype=="DAA2"])),
           (x$Count[x$distance=="(50,Inf]" & x$Subtype =="DAA2"]/sum(x$Count[x$Subtype=="DAA2"])),
           (x$Count[x$distance=="[0,10]" & x$Subtype =="DN"]/sum(x$Count[x$Subtype=="DN"])),
           (x$Count[x$distance=="(10,20]" & x$Subtype =="DN"]/sum(x$Count[x$Subtype=="DN"])),
           (x$Count[x$distance=="(20,50]" & x$Subtype =="DN"]/sum(x$Count[x$Subtype=="DN"])),
           (x$Count[x$distance=="(50,Inf]" & x$Subtype =="DN"]/sum(x$Count[x$Subtype=="DN"])))
  summary_s100 <- c((y$Count[y$distance=="[0,10]" & y$Subtype =="DAA1"]/sum(y$Count[y$Subtype=="DAA1"])),
                   (y$Count[y$distance=="(10,20]" & y$Subtype =="DAA1"]/sum(y$Count[y$Subtype=="DAA1"])),
                   (y$Count[y$distance=="(20,50]" & y$Subtype =="DAA1"]/sum(y$Count[y$Subtype=="DAA1"])),
                   (y$Count[y$distance=="(50,Inf]" & y$Subtype =="DAA1"]/sum(y$Count[y$Subtype=="DAA1"])),
                   (y$Count[y$distance=="[0,10]" & y$Subtype =="DAA2"]/sum(y$Count[y$Subtype=="DAA2"])),
                   (y$Count[y$distance=="(10,20]" & y$Subtype =="DAA2"]/sum(y$Count[y$Subtype=="DAA2"])),
                   (y$Count[y$distance=="(20,50]" & y$Subtype =="DAA2"]/sum(y$Count[y$Subtype=="DAA2"])),
                   (y$Count[y$distance=="(50,Inf]" & y$Subtype =="DAA2"]/sum(y$Count[y$Subtype=="DAA2"])),
                   (y$Count[y$distance=="[0,10]" & y$Subtype =="DN"]/sum(y$Count[y$Subtype=="DN"])),
                   (y$Count[y$distance=="(10,20]" & y$Subtype =="DN"]/sum(y$Count[y$Subtype=="DN"])),
                   (y$Count[y$distance=="(20,50]" & y$Subtype =="DN"]/sum(y$Count[y$Subtype=="DN"])),
                   (y$Count[y$distance=="(50,Inf]" & y$Subtype =="DN"]/sum(y$Count[y$Subtype=="DN"])))
  
  df <- data.frame(Distance = c((x$distance),y$distance),
                   Subtype = c(x$Subtype,y$Subtype),
                   Sample=c(factor(x$id),factor(y$id)),
                   ratios = c(summary_igf,summary_s100))
  
  color.codes <- c("goldenrod1","red","black")
  zone <- levels(factor(df$Subtype))
  
  p <-  ggplot(df, aes(x=Distance,y=ratios,fill=Subtype)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(color=Sample, group=Subtype),size=3,
               position = position_dodge2(width =.75)) +
    scale_fill_manual(values=setNames(color.codes, zone)) +
    theme_pubclean()+
    ylab("Ratio of Astrocyte type") +
    xlab(paste("Distance to nearest plaque center >",cutoff,"um"))+
    theme(axis.title = element_text(size=24, face = "bold"),
          axis.text = element_text(size=14),
          legend.position = "right")
  
  StandardLabelPlots[[as.character(cutoff)]] <- p
}
StandardLabelPlots$`0`
allStandard <- patchwork::wrap_plots(StandardLabelPlots,ncol = 3)
