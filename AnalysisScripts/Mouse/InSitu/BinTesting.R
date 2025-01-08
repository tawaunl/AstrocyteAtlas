# Bin Testing 
# Load packages ===========================================
pkgs <-c("ggplot2","ggpubr","viridis","gridExtra",
         "dplyr","tidyverse")
pacman::p_load(pkgs,character.only = TRUE)

igfbp5.dir <- file.path("~","Documents/InsituAnalysis/Igfbp5_results")
s100.dir <- file.path("~","Documents/InsituAnalysis/S100a6_results")
tg_animals_igf <- c("20_01","20_02","21_01","21_02")
tg_animals_s100 <- c("31_02", "31_03", "32_01", "32_02" ,"33_02")
p_cutoffs <- c(0,5,10,15,25,50) # plaque size cutoffs
breaks <- list(orig=c(0,10,20,50,Inf),
               By5 = c(0,5,10,15,20,25,Inf),
               By10 = c(0,10,20,30,40,50,Inf)) # distance bins

# load Igfbp5 data =================
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

BinPlots <- lapply(breaks,function(y){
  binnedData <- lapply(astro_data[tg_animals_igf],function(x){
    x$dis2nearestPlaque <- x[paste0("dis2nearestPlaque_",0)]
    test <- x %>%
      group_by(dis2nearestPlaque, normastroType) %>%
      summarise(Count = n())
    colnames(test)[1] <- "dis2nearestPlaque"
    aggregate(test[,-c(1,2)], by = list(distance=cut(unlist(as.vector(test$dis2nearestPlaque)),
                                                     y, include.lowest=T),
                                        Subtype = test$normastroType ), FUN = sum )
  })
  
  
  x <- list_rbind(x=binnedData, 
                  names_to = "id") %>% arrange(c(distance)) %>% arrange(c(Subtype)) 
  
  idv <- vector()
  for(type in levels(factor(x$Subtype))){
    for(dis in levels(factor(x$distance))){
      idv = c(idv,(x$Count[x$distance== dis & x$Subtype == type]/sum(x$Count[x$Subtype==type])))
    }
  }
  
  df <- data.frame(Distance = (x$distance),
                   Subtype = x$Subtype,
                   Sample=factor(x$id),
                   ratios = idv)
  
  color.codes <- c("goldenrod1","red","black")
  zone <- levels(factor(df$Subtype))
  
  ggplot(df, aes(x=Distance,y=ratios,fill=Subtype)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(color=Sample, group=Subtype),size=3,
               position = position_dodge2(width =.75)) +
    scale_fill_manual(values=setNames(color.codes, zone)) +
    theme_pubclean()+
    ylab("Ratio of Astrocyte type") +
    xlab(paste("Distance to nearest plaque "))+
    theme(axis.title = element_text(size=24, face = "bold"),
          axis.text = element_text(size=14),
          legend.position = "right")
  
  
})

patchwork::wrap_plots(BinPlots,ncol=1)


# S100 ======================
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

BinPlots <- lapply(breaks,function(y){
  binnedData <- lapply(astro_data[tg_animals_s100],function(x){
    x$dis2nearestPlaque <- x[paste0("dis2nearestPlaque_",0)]
    test <- x %>%
      group_by(dis2nearestPlaque, normastroType) %>%
      summarise(Count = n())
    colnames(test)[1] <- "dis2nearestPlaque"
    aggregate(test[,-c(1,2)], by = list(distance=cut(unlist(as.vector(test$dis2nearestPlaque)),
                                                     y, include.lowest=T),
                                        Subtype = test$normastroType ), FUN = sum )
  })
  
  
  x <- list_rbind(x=binnedData, 
                  names_to = "id") %>% arrange(c(distance)) %>% arrange(c(Subtype)) 
  
  idv <- vector()
  for(type in levels(factor(x$Subtype))){
    for(dis in levels(factor(x$distance))){
      idv = c(idv,(x$Count[x$distance== dis & x$Subtype == type]/sum(x$Count[x$Subtype==type])))
    }
  }
  
  df <- data.frame(Distance = (x$distance),
                   Subtype = x$Subtype,
                   Sample=factor(x$id),
                   ratios = idv)
  
  color.codes <- c("goldenrod1","red","black")
  zone <- levels(factor(df$Subtype))
  
  ggplot(df, aes(x=Distance,y=ratios,fill=Subtype)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(color=Sample, group=Subtype),size=3,
               position = position_dodge2(width =.75)) +
    scale_fill_manual(values=setNames(color.codes, zone)) +
    theme_pubclean()+
    ylab("Ratio of Astrocyte type") +
    xlab(paste("Distance to nearest plaque "))+
    theme(axis.title = element_text(size=24, face = "bold"),
          axis.text = element_text(size=14),
          legend.position = "right")
  
  
})

patchwork::wrap_plots(BinPlots,ncol=1)
