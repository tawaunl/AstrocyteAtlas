library(tidyverse)
library(caret)
library(terra)
library(dplyr)
library(fields)

#########################################
# ### (A) IMPORT FILES ----------------------------------------------------


#############  NOTE: (1) files MUST be in name format 20_20240117_071551_01.txt, where 20 is the sample number and 01 is the replicate within sample
#############        (2) change file path for file location for variable "det_files" and "det_data"
#############        (3) creates names for datasets from from input det_files and stores in  det_names
#############        (4) loads data file into dataframe  :: det_temp[i]
dir <- file.path("~","Documents/InsituAnalysis/Igfbp5_results")
det_files <- list.files(file.path(dir,"detection_results"))                                                  # Extract file names from this file path
det_files 
det_names<-c()                                                                           # storage for names of datasets
det_data<-list()                                                                         # storage for datasets


for(i in 1:length(det_files)) {                                                          
  
  # (1) create dataset name from filename/datafile, store in variable det_temp
  if(nchar(det_files[i])>65){
    det_names<-append(det_names, paste0(substring(det_files[i], 1, 2), sep = "_", substring(det_files[i], 67, 68))) # creates+stores names for dataset size:length(det_data)
  }else{
    det_names<-append(det_names, paste0(substring(det_files[i], 1, 2), sep = "_", substring(det_files[i], 60, 61))) # creates+stores names for dataset size:length(det_data)
  }
  
  # (2.2) try read table for .txt files 
  
  det_data[[i]]<- as_tibble(read.table(paste0(dir,"/detection_results/",  
                                              det_files[i]), sep = '\t', header = TRUE))                                        # stores each .txt data set 
  
} ## LOOP END

names(det_data) <- det_names

## Read Plaque data =-------------------
plaque_data <- list()
plaque_files <- list.files(file.path(dir,"annotation_results_plaques")) 
plaque_names <- c()
for(i in 1:length(plaque_files)) {                                                          
  # (1) create dataset name from filename/datafile, store in variable det_temp
  if(nchar(plaque_files[i])>65){
    plaque_names<-append(plaque_names, paste0(substring(plaque_files[i], 1, 2), sep = "_", substring(plaque_files[i], 67, 68))) # creates+stores names for dataset size:length(det_data)
  }else{
    plaque_names<-append(plaque_names, paste0(substring(plaque_files[i], 1, 2), sep = "_", substring(plaque_files[i], 60, 61))) # creates+stores names for dataset size:length(det_data)
  }
  
  plaque_data[[i]]<- as_tibble(read.table(paste0(dir,"/annotation_results_plaques/",  
                                                 plaque_files[i]), sep = '\t', header = TRUE))                                        # stores each .txt data set 
  
} 
names(plaque_data) <- plaque_names

# ### (B) GROUP CELLS AND DETECTIONS --------------------------------------


##############(1)Add values to 'Names' to group cells and associated detections
##############(2) adds letter A to Cell1, Detection1, Detection2, letter B to Cell2, Detection1, Detection2 (and so on)
##############(3) this is unique within a single det_data within the list but not between det_data1 and det_data2 for example (i.e. use filename or another var to distinguish before query for cells/detections by Name) )

n=1

for (n in 1:length(det_data)){                                                   #loop through each dataset in list
  
  cell_inds<-0
  cell_diffs<-0
  cell_groups<-0 
  det_temp<-data.frame()
  
  det_temp<-det_data[[n]]                                                          #pass current data frame to temp var.
  
  cell_inds<-which(det_temp$Object.type == 'Cell')
  cell_diffs <- diff(c(cell_inds, length(det_temp$Object.type)+1))
  cell_groups <- rep(seq_along(cell_inds),cell_diffs)
  det_data[[n]]$Name<- cell_groups
  
}


# ### (C) ANALYSIS OF ASTROCYTE DISTANCE TO PLAQUES--------------------------


###############(1) Analyze FITC signal based on selected parameter from QuPath data across all samples
###############(2) subset cells,astrocytes, and plaque coordinates
###############(3) calculate distance metrics between cells and plaques
###############(4) Save distance to nearest plaque for plotting                       ##


### #  2. Nucleus..FITC.mean


astro_cutoff<-c(7000)
index<-0
j=1
i=1


df_cells2<-list()
df_astrocytes2<-list()


dim_cells2<-list()
dim_ast2<-list()
dim_ratio2<-list()


cutoff.tracker<-c()
sample.tracker<-c()

df_table2<-data.frame()
cell_maps <-list()

distance_plots_ch3 <-list()
distance_plots_ch4 <-list()
astro_data <- list()
cells_data <- list()
orig_plaques <- list()
p_cutoffs <- c(0,5,10,15,25,50) # plaque size cutoffs
tg_animals <- c("20_01","20_02","21_01","21_02")
# runs for i= number of input datasets
for (dataset in 1:length(det_data)){     # runs for i= number of input datasets
  
  cells <-  filter(det_data[[dataset]], Object.type  == "Cell")  
  
  plaques <- filter(det_data[[dataset]], Classification  == "plaques")  
  
  if("Nucleus..FITC.mean" %in% colnames(cells)){
    astros <-  filter(cells, Nucleus..FITC.mean > astro_cutoff[j])
  } else{
    astros <-  filter(cells, Nucleus..FITC..C2..mean > astro_cutoff[j])
  }
  
  cells <- cells[-which(cells$Object.ID %in% astros$Object.ID),]
  
  cell_coords <- data.frame(x=cells$Centroid.X.µm,
                            y=cells$Centroid.Y.µm,
                            size=as.vector(predict((preProcess(as.data.frame(cells$Cell..Area), method=c("range"))), as.data.frame(cells$Cell..Area))),
                            color="red",
                            cell="cell")
  
  colnames(cell_coords)[3]="size"
  
  astro_coords <- data.frame(x = astros$Centroid.X.µm,
                             y = astros$Centroid.Y.µm,
                             size = predict(
                               (preProcess(as.data.frame(astros$Cell..Area), method=c("range"))),
                               as.data.frame(astros$Cell..Area)),
                             color ="yellow",
                             cell = "astro",
                             row.names=astros$Object.ID)
  colnames(astro_coords)[3]="size"
  
  if(det_names[dataset] %in% tg_animals){
    orig_plaques[[names(det_data)[dataset]]] <- plaques
    objs <- plaque_data[[names(det_data)[dataset]]]
    
    sizes = predict(
      (preProcess(as.data.frame(objs$Area.µm.2), method=c("range"))),
      as.data.frame(objs$Area.µm.2))
    
    plaque_coords <- data.frame(x=objs$Centroid.X.µm,
                                y=objs$Centroid.Y.µm,
                                size=predict(
                                  (preProcess(as.data.frame(objs$Area.µm.2), method=c("range"))),
                                  as.data.frame(objs$Area.µm.2)),
                                color="blue",cell="plaque",
                                row.names = objs$Object.ID)
    colnames(plaque_coords)[3] <- "size"
    
    # plot map 
    
    plots <- list()
    for(cutoff in p_cutoffs){
      stations <- rbind(cell_coords,
                        astro_coords,
                        plaque_coords[which(objs$Area.µm.2 >= cutoff),])
      plots[[as.character(cutoff)]] <- ggplot(stations,aes(x=x,y=y))+ geom_point(aes(size=size,colour=color),show.legend = FALSE) +
        scale_radius(range = c(0, 1)) + 
        scale_colour_identity() + 
        ggtitle(paste(det_names[dataset],"cutoff >",cutoff,"um"))+theme_classic()+
        theme(plot.title = element_text(size = 20, face="bold",hjust = 0.5))
    }
    
    cell_maps[[det_names[dataset]]] <- patchwork::wrap_plots(plots,ncol=3)
    
    # calculate all distances between plaques and astros
    for(cutoff in p_cutoffs){
      p_sub <- plaque_coords[which(objs$Area.µm.2 >= cutoff),]
      dis <- rdist(astro_coords[,1:2],p_sub[,1:2]) # Euclidean distance
      dis_cells <- rdist(cell_coords[,1:2],p_sub[,1:2])
      # Add min distances back to astrocyte coords
      astro_coords[paste0("dis2nearestPlaque_",cutoff)] <- apply(dis,1,min)
      cell_coords[paste0("dis2nearestPlaque_",cutoff)] <- apply(dis_cells,1,min)
    }
    
    # add in metadata for spot counts
    astro_coords <- cbind(astro_coords,astros)
    cell_coords <- cbind(cell_coords,cells)
    
    astro_data[[det_names[dataset]]] <- astro_coords
    cells_data[[det_names[dataset]]] <- cell_coords
    
  } else{
    stations <- rbind(cell_coords,astro_coords)
    cell_maps[[det_names[dataset]]] <- ggplot(stations,aes(x=x,y=y))+ geom_point(aes(size=size,colour=color),show.legend = FALSE) +
      scale_radius(range = c(0, 1)) + 
      scale_colour_identity() + ggtitle(det_names[dataset])+theme_classic()+
      theme(plot.title = element_text(size = 20, face="bold",hjust = 0.5))
    
    cells_data[[det_names[dataset]]] <-cbind(cell_coords,cells)
    astro_data[[det_names[dataset]]] <-cbind(astro_coords,astros)
  }
}

# (D) Correlation between marker expression --------------------
# Are the cells that are high for DAA1 also higher for DAA2 ?
## Scatter plot DAA1 and DAA2 spots ==============
library(ggpubr)
library(gridExtra)

scatters <- list()

scatters <- lapply(astro_data, function(x){
  ggscatter(x, x = "Subcellular..Channel.3..Num.spots.estimated",
            y = "Subcellular..Channel.4..Num.spots.estimated", 
            color = "Nucleus..Area",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Estimated Igfbp5 spots",
            ylab = "Estimated Apoe spots")
  
})

do.call("grid.arrange", c(scatters, list(ncol=3)))

## Corrected scatters =========================
astro_data <- lapply(astro_data, function(x){
  
  x$normch3 <- x$Subcellular..Channel.3..Num.spots.estimated / x$Nucleus..Area
  x$normch4 <- x$Subcellular..Channel.4..Num.spots.estimated / x$Nucleus..Area
  x 
})

norm_scatters <- lapply(astro_data, function(x){
  ggscatter(x, x = "normch3",
            y = "normch4", 
            color = "Nucleus..Area",
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Estimated Igfbp5 spots",
            ylab = "Estimated Apoe spots")
  
})

do.call("grid.arrange", c(norm_scatters, list(ncol=3)))

#(E) TauPS2APP vs NonTG expression -------------------
library(viridis)
for(dataset in 1:length(astro_data)){
  if(names(astro_data[dataset]) %in% tg_animals){
    astro_data[[dataset]]$Genotype <- "TauPS2APP"
  }
  else{
    astro_data[[dataset]]$Genotype <- "NonTG"
  }
}

expressionMeans <- data.frame(matrix(ncol = 4))

for (dataset in 1:length(astro_data)) {
  name <- names(astro_data[dataset])
  genotype <- astro_data[[dataset]]$Genotype[1]
  meanch3 <- mean(astro_data[[dataset]]$normch3)
  meanch4 <- mean(astro_data[[dataset]]$normch4)
  expressionMeans[dataset,] <- c(name,genotype,meanch3,meanch4)
}
colnames(expressionMeans) <- c("Sample","Genotype","Igfbp5","Apoe")

expressionMeans <- pivot_longer(expressionMeans,
                                cols=c(Igfbp5,Apoe),
                                names_to = "ChannelName")


expressionMeans$value <- as.numeric(expressionMeans$value)
ggplot(expressionMeans, aes(x=Genotype, y=value,fill=Genotype)) + theme_classic() +
  geom_boxplot(outlier.shape=NA) + scale_fill_viridis(discrete = TRUE, alpha=1) +
  geom_jitter(width = 0.4,size=3,aes(color=Sample))  +
  facet_wrap(~factor(ChannelName),scales = "free", ncol=2) + xlab("Disease Label") +
  theme(axis.text.x = element_text(angle=75,vjust = 0.5),axis.text.y = element_text(size=20)) +
  ylab("Mean Expression") +  
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),legend.text = element_text(size=14),
        legend.title = element_text(size=16)) 


# (F) Expression Relation to pathology ---------------------------
breaks <- c(0,10,20,50,Inf)

RelationPlots <- list()
for (cutoff in p_cutoffs) {
  plots <- list()
  binnedData <- lapply(astro_data[tg_animals],function(x){
    ag <- aggregate(x,by=list(distance=cut(unlist(as.vector(x[paste0("dis2nearestPlaque_",cutoff)])),
                                           breaks, include.lowest=T),
                              Cy3spots=x$Subcellular..Channel.3..Num.spots.estimated,
                              Cy5spots=x$Subcellular..Channel.4..Num.spots.estimated),
                    FUN=mean)
    y <- aggregate(ag[,-1], by=list(binnedDistance=ag$distance),FUN=mean)
    y
  })
  
  x <- list_rbind(x=binnedData, 
                  names_to = "id") %>% arrange(c(distance))
  
  p1 <- ggplot(x,aes(x=binnedDistance,y=Cy3spots)) +
    geom_boxplot(aes(fill=binnedDistance),outlier.shape = NA) +
    geom_point(aes(color=id),position=position_dodge(width=0.75),size=3) +
    scale_fill_viridis(discrete = TRUE, alpha=1,name="Distance",
                       breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                       labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) + 
    scale_x_discrete(breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                     labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) +
    xlab(paste("Distance to nearest plaque center >",cutoff,"um")) + 
    ylab("Avg. DAA1(Igfbp5) spots in astrocyte")+
    ggtitle("Average Igfbp5 expression in astrocytes")+
    theme_pubclean() +
    theme(plot.title = element_text(hjust=0.5,face="bold",size=20),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = "right") 
  
  p2 <- ggplot(x,aes(x=binnedDistance,y=Cy5spots)) +
    geom_boxplot(aes(fill=binnedDistance),outlier.shape = NA) +
    geom_point(aes(color=id),position=position_dodge(width=0.75),size=3) +
    scale_fill_viridis(discrete = TRUE, alpha=1,name="Distance",
                       breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                       labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) + 
    scale_x_discrete(breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                     labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) +
    xlab(paste("Distance to nearest plaque center >",cutoff,"um")) +
    ylab("Avg. DAA2 (Apoe) spots in astrocyte") +
    ggtitle("Average Apoe expression in astrocytes") +
    theme_pubclean() +
    theme(plot.title = element_text(hjust=0.5,face="bold",size=20),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = "right") 
  RelationPlots[[as.character(cutoff)]] <- patchwork::wrap_plots(p1,p2,ncol=2)
  
}

##  Region specificity -------------
RegionPlots <- list()
for (cutoff in p_cutoffs) {
  binnedData_region <- lapply(astro_data[tg_animals],function(x){
    ag <- aggregate(x,by=list(distance=cut(unlist(as.vector(x[paste0("dis2nearestPlaque_",cutoff)])),
                                           breaks, include.lowest=T),
                              Cy3spots=x$Subcellular..Channel.3..Num.spots.estimated,
                              Cy5spots=x$Subcellular..Channel.4..Num.spots.estimated,
                              Region=x$Parent),
                    FUN=mean)
    
    aggregate(ag[,-c(1,4)], by=list(Distance=ag$distance,
                                    Region=ag$Region),FUN=mean)
    
  })
  
  region_data <- list_rbind(x=binnedData_region, 
                            names_to = "id") %>% arrange(c(Distance)) %>% arrange(c(Region)) 
  
  p1 <- ggplot(region_data,aes(x=Distance,y=Cy3spots)) +
    geom_boxplot(aes(fill=Distance),outlier.shape = NA) +
    geom_point(aes(color=id),position=position_dodge(width=0.75),size=3) +
    facet_wrap(~Region,ncol=3) +
    scale_fill_viridis(discrete = TRUE, alpha=1,name="Distance",
                       breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                       labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) +
    xlab(paste("Distance to nearest plaque center >",cutoff,"um")) + 
    ylab("Avg. DAA1(Igfbp5) spots in astrocyte")+
    ggtitle("Average Igfbp5 expression in astrocytes")+
    theme_pubclean() +
    theme(plot.title = element_text(hjust=0.5,face="bold",size=20),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "right",
          panel.spacing = unit(0,"lines"))
  
  p2 <- ggplot(region_data,aes(x=Distance,y=Cy3spots)) +
    geom_boxplot(aes(fill=Distance),outlier.shape = NA) +
    geom_point(aes(color=id),position=position_dodge(width=0.75),size=3) +
    facet_wrap(~Region,ncol=3) +
    scale_fill_viridis(discrete = TRUE, alpha=1,name="Distance",
                       breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                       labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) +
    xlab(paste("Distance to nearest plaque center >",cutoff,"um")) + 
    ylab("Avg. DAA1(Igfbp5) spots in astrocyte")+
    ggtitle("Average Igfbp5 expression in astrocytes")+
    theme_pubclean() +
    theme(plot.title = element_text(hjust=0.5,face="bold",size=20),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "right",
          panel.spacing = unit(0,"lines"))
  RegionPlots[[as.character(cutoff)]] <- patchwork::wrap_plots(p1,p2,ncol=2)
  
}




# (G) Label Astros ============
## Standard label ========================

# we will iterate over plaque size with standard cutoffs
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

data <- bind_rows(astro_data, .id = "column_label")
png(file.path(dir,"Igfbp5_ScatterSubtype.png"),
    width =800,height = 500)

ggscatter(data, x = "normch3",
          y = "normch4",color="normastroType",
          add = "reg.line", conf.int = TRUE,
          palette =c("goldenrod1","red","black"),
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Estimated Igfbp5 spots",
          ylab = "Estimated Apoe spots")
dev.off()

StandardLabelPlots <- list()

sub <- data[data$column_label %in% tg_animals_igf,]
cells <- data.frame(sub[,c(1:7,14,15,103,106,116,117,119,176)])

cells_21 <- cells[cells$column_label =="21_01",]

ggscatter(cells_21, x = "normch3",
          y = "normch4",color="normastroType",
          add = "reg.line", conf.int = TRUE,
          palette =c("goldenrod1","red","black"),
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Estimated Igfbp5 spots",
          ylab = "Estimated Apoe spots")

for (cutoff in p_cutoffs) {
  binnedData <- lapply(astro_data[tg_animals],function(x){
    x$dis2nearestPlaque <- x[paste0("dis2nearestPlaque_",cutoff)]
    test <- x %>%
      group_by(dis2nearestPlaque, normastroType) %>%
      summarise(Count = n())
    colnames(test)[1] <- "dis2nearestPlaque"
    aggregate(test[,-c(1,2)], by = list(distance=cut(unlist(as.vector(test$dis2nearestPlaque)),
                                                     breaks, include.lowest=T),
                                        Subtype = test$normastroType ), FUN = sum )
  })
  
  
  x <- list_rbind(x=binnedData, 
                  names_to = "id") %>% arrange(c(distance)) %>% arrange(c(Subtype)) 
  
  idv <- c((x$Count[x$distance=="[0,10]" & x$Subtype =="DAA1"]/sum(x$Count[x$Subtype=="DAA1"])),
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
  
  df <- data.frame(Distance = (x$distance),
                   Subtype = x$Subtype,
                   Sample=factor(x$id),
                   ratios = idv)
  
  color.codes <- c("goldenrod1","red","black")
  zone <- levels(factor(df$Subtype))
  
  p <-  ggplot(df, aes(x=Distance,y=ratios,fill=Subtype)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(color=Sample, group=Distance),size=3,
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

png(file.path(dir,"Igfbp5_DistancePlotSubtype.png"),
    widt=1000,height = 500)

StandardLabelPlots$`0`

dev.off()

save(list=c("cell_maps","astro_data","cells_data","allStandard",
            "plaque_data","orig_plaques","RegionPlots","RelationPlots"),
     file=file.path(dir,"data.files.Rdata"))



## iterate over cutoffs for cell calls  =======================
calcMeans <- function(cutoff){
  binnedData <- lapply(astro_data[tg_animals],function(x){
    x$dis2nearestPlaque <- x[paste0("dis2nearestPlaque_",cutoff)]
    test <- x %>%
      group_by(dis2nearestPlaque, normastroType) %>%
      summarise(Count = n())
    colnames(test)[1] <- "dis2nearestPlaque"
    aggregate(test[,-c(1,2)], by = list(distance=cut(unlist(as.vector(test$dis2nearestPlaque)),
                                                     breaks, include.lowest=T),
                                        Subtype = test$normastroType ), FUN = sum )
  })
  
  
  x <- list_rbind(x=binnedData, 
                  names_to = "id") %>% arrange(c(distance)) %>% arrange(c(Subtype))
  
  
  
  idv <- c((x$Count[x$distance=="[0,10]" & x$Subtype =="DAA1"]/sum(x$Count[x$Subtype=="DAA1"])),
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
  
  
  df <- data.frame(Distance = (x$distance),
                   Subtype = x$Subtype,
                   Sample=factor(x$id),
                   ratios = idv)
  return(df)
}
xcutoffs <- c(1,1.25,1.5,1.75,2)
ycutoffs <- c(1.5,1.75,2,2.25,2.5)
iterativePlots <- list()
for(cutoff in p_cutoffs){
  for(igf in xcutoffs){
    for (apoe in ycutoffs) {
      astro_data <- lapply(astro_data, function(x){
        x$normch3 <- x$Subcellular..Channel.3..Num.spots.estimated / x$Nucleus..Perimeter
        x$normch4 <- x$Subcellular..Channel.4..Num.spots.estimated / x$Nucleus..Perimeter
        x$normastroType <- "DN"
        x$normastroType[x$normch3 >= igf &
                          x$normch4 <= 1.5] <- "DAA1"
        x$normastroType[x$normch3 < 0.5 &
                          x$normch4 >= apoe ] <- "DAA2"
        x
      })
      data <- bind_rows(astro_data, .id = "column_label")
      
      scat <- ggscatter(data, x = "normch3",
                        y = "normch4",color="normastroType",
                        add = "reg.line", conf.int = TRUE,
                        palette =c("goldenrod1","red","black"),
                        cor.coef = TRUE, cor.method = "pearson",
                        xlab = "Estimated Igfbp5 spots",
                        ylab = "Estimated Apoe spots")
      
      plotdata <- calcMeans(cutoff=cutoff)
      p <-  ggplot(df, aes(x=Distance,y=ratios,fill=Subtype)) +
        geom_boxplot(outlier.shape = NA) + 
        geom_point(aes(color=Sample, group=Distance),size=3,
                   position = position_dodge2(width =.75)) +
        scale_fill_manual(values=setNames(color.codes, zone)) +
        theme_pubclean()+
        ylab("Ratio of Astrocyte type") +
        xlab(paste("Distance to nearest plaque center >",cutoff,"um"))+
        theme(axis.title = element_text(size=24, face = "bold"),
              axis.text = element_text(size=14),
              legend.position = "right")
      plot_row <- plot_grid(scat, p)
      
      title <- ggdraw() + 
        draw_label(
          paste("Apoe cutoff=",apoe,"   ",
                "Igfbp5 cutoff =",igf), fontface = 'bold', x = 0,hjust = -1.5) +
        theme(plot.margin = margin(0, 0, 0, 7))
      
      iterativePlots[[paste0("PlaqueeSize_",cutoff)]] <- c(iterativePlots[[paste0("PlaqueeSize_",cutoff)]], list(plot_grid( title, plot_row, ncol = 1,
                                                                                                                            rel_heights = c(0.1, 1))))
    }
  }
}

# Final Save ===================================
save(list=c("cell_maps","astro_data","cells_data","allStandard",
            "plaque_data","orig_plaques","RegionPlots","RelationPlots"),
     file=file.path(dir,"data.files.Rdata"))


# plot iteratively into PDF
for(cutoff in p_cutoffs){
  pdf(paste0(dir,"/DistancePlots_PlaqueCutoff_",cutoff,".pdf"),width = 14 )
  for(plot in 1:length(iterativePlots[[paste0("PlaqueeSize_",cutoff)]])){
    print(iterativePlots[[paste0("PlaqueeSize_",cutoff)]][[plot]]) }
  dev.off()
}





# Misc Plots =========================================================
### Density plot across cells =========================================

density_data <- lapply(astro_data[tg_animals], function(x){
  df <- data.frame(Name      = x$Name,
                   Region    = x$Parent,
                   DAA1Spots = x$Subcellular..Channel.3..Num.single.spots,
                   DAA2Spots = x$Subcellular..Channel.4..Num.single.spots,
                   PlaqueDis = x$dis2nearestPlaque)
  df1 <- df %>% 
    pivot_longer(
      cols      = DAA1Spots:DAA2Spots,
      names_to  = "channel",
      values_to = "estimatedSpots"
    )
  df1
})
density_plots <- lapply(density_data,function(x){
  
  ggplot(data=x, aes(x=estimatedSpots, group=channel, fill=channel)) +
    geom_density(adjust=1.5, alpha=.4)
  
})

library(gridExtra)
library(ggpubr)
title=text_grob("Igfbp5 Results", size = 20, face = "bold") 

png(file.path(dir,"DensityPlots_AllSamples.png"),width = 800,height=500)
do.call("grid.arrange", c(density_plots, list(ncol=2,top=title)))
dev.off()
### Region ==========
density_plots_region <- lapply(density_plots, function(x){
  x + facet_wrap(~Region)
})

png(file.path(dir,"DensityPlots_Region_AllSamples.png"),width = 800,height=500)
do.call("grid.arrange", c(density_plots_region, list(nrow=2,top=title)))
dev.off()
### Proportions ==========

pp <- lapply(astro_data[tg_animals], function(x){
  df <- data.frame(Name      = x$Name,
                   Region    = x$Parent,
                   DAA1Spots = x$Subcellular..Channel.3..Num.single.spots,
                   DAA2Spots = x$Subcellular..Channel.4..Num.single.spots,
                   PlaqueDis = x$dis2nearestPlaque,
                   SpotPro  = x$Subcellular..Channel.3..Num.single.spots/x$Subcellular..Channel.4..Num.single.spots)
  df
})

data <- map_df(.x = seq(1:length(pp)),
               .f = function(x) bind_rows(pp) %>%
                 mutate(id = row_number(),
                        df = x)) %>% arrange(id) 

ggplot(data=data, aes(x=SpotPro)) +
  geom_density(adjust=1.5, alpha=.4) + ggtitle("Aggregated Igfbp5 Proportions") + 
  theme(plot.title = element_text(size=20,face="bold",hjust=0.5))

### look at aggregated across animals =========================
density_data_long <- map_df(.x = seq(1:length(density_data)),
                            .f = function(x) bind_rows(density_data) %>%
                              mutate(id = row_number(),
                                     df = x)) %>% arrange(id) 

png(file.path(dir,"DensityPlots_Aggregated.png"),width = 800,height=500)
ggplot(data=density_data_long, aes(x=estimatedSpots, group=channel, fill=channel)) +
  geom_density(adjust=1.5, alpha=.4) + ggtitle("Aggregated Igfbp5 Data") + 
  theme(plot.title = element_text(size=20,face="bold",hjust=0.5))
dev.off()

png(file.path(dir,"DensityPlots_AggregatedRegion.png"),width = 800,height=500)
ggplot(data=density_data_long, aes(x=estimatedSpots, group=channel, fill=channel)) +
  geom_density(adjust=1.5, alpha=.4) + ggtitle("Aggregated Igfbp5 Data") + 
  theme(plot.title = element_text(size=20,face="bold",hjust=0.5)) +
  facet_wrap(~Region)
dev.off()



## collapse all astrocytes
library(dplyr)
data <- bind_rows(astro_data, .id = "column_label")
ggscatter(data, x = "Subcellular..Channel.3..Num.single.spots",
          y = "Subcellular..Channel.4..Num.single.spots", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Estimated Igfbp5 spots",
          ylab = "Estimated Apoe spots")

## Boxplots of dispersions  ========
ggplot(data.frame(data$Subcellular..Channel.3..Num.single.spots),aes(y=data.Subcellular..Channel.3..Num.single.spots)) +
  geom_boxplot() +ggtitle("Igfbp5") + ylab("# spots") + theme_bw() + 
  theme(plot.title = element_text(hjust=0.5,size=22,face = "bold"))
ggplot(data.frame(data$Subcellular..Channel.4..Num.single.spots),aes(y=data.Subcellular..Channel.4..Num.single.spots)) +
  geom_boxplot() +ggtitle("Apoe") + ylab("# spots") + theme_bw() + 
  theme(plot.title = element_text(hjust=0.5,size=22,face = "bold"))

#cutoffs determined based on third quantile
cutoffc3 <-10
cutoffc4 <-10

## Label Astro types =================
astro_data <- lapply(astro_data, function(x){
  x$astroType <- "DN"
  x$astroType[x$Subcellular..Channel.3..Num.single.spots >= cutoffc3 &
                x$Subcellular..Channel.4..Num.single.spots < cutoffc4 ] <- "DAA1"
  x$astroType[x$Subcellular..Channel.3..Num.single.spots < cutoffc3 &
                x$Subcellular..Channel.4..Num.single.spots >= cutoffc4 ] <- "DAA2"
  x$astroType[x$Subcellular..Channel.3..Num.single.spots >= cutoffc3 &
                x$Subcellular..Channel.4..Num.single.spots >= cutoffc4 ] <- "DP"
  x
})
data <- bind_rows(astro_data, .id = "column_label")
data$totalSpotsch3 <- data$Subcellular..Channel.3..Num.single.spots + data$Subcellular..Channel.3..Num.clusters
data$totalSpotsch4 <- data$Subcellular..Channel.4..Num.single.spots + data$Subcellular..Channel.4..Num.clusters

data$astroType <- "DN"
data$astroType[data$totalSpotsch3 >= cutoffc3 &
                 data$totalSpotsch4 < cutoffc4 ] <- "DAA1"
data$astroType[data$totalSpotsch3 < cutoffc3 &
                 data$totalSpotsch4 >= cutoffc4 ] <- "DAA2"
data$astroType[data$totalSpotsch3 >= cutoffc3 &
                 data$totalSpotsch4 >= cutoffc4 ] <- "DP"
x
data$normch3 <- data$Subcellular..Channel.3..Num.spots.estimated / data$Nucleus..Perimeter
data$normch4 <- data$Subcellular..Channel.4..Num.spots.estimated / data$Nucleus..Perimeter
data$shuffCh3 <- sample(data$normch3)
data$shuffCh4 <- sample(data$normch4)

## colored scatter ===============
ggscatter(data, x = "normch3",
          y = "shuffCh4",color="Nucleus..Area",
          add = "reg.line", conf.int = TRUE, palette = "Reds",
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Estimated Igfbp5 spots",
          ylab = "Estimated Apoe spots")+ scale_color_gradient(low="yellow",high="red")


ggplot(data[order(data$Nucleus..Area,decreasing = F),], aes(x=Subcellular..Channel.3..Num.spots.estimated,y=Subcellular..Channel.4..Num.spots.estimated))+
  geom_point(aes(color=Nucleus..Area)) + scale_color_gradient(low="yellow",high="red")

# count cell numbers across tissue

ex <- astro_data[[1]]
daa1 <- sum(ex$astroType=="DAA1")

cell_sums <- lapply(astro_data, function(x){
  data.frame(DP   = sum(x$astroType=="DP"),
             DN   = sum(x$astroType=="DP"),
             DAA1 = sum(x$astroType=="DAA1"),
             DAA2 = sum(x$astroType=="DAA2"))
})
sums_data <- bind_rows(cell_sums, .id = "column_label")

plotdata <- pivot_longer(sums_data,cols=starts_with("D"))

ggplot(plotdata,aes(x=name,y=value)) +geom_boxplot() +geom_jitter(aes(color=column_label))



## Ratio analysis ================
astro_data <- lapply(astro_data,function(x){
  x$astroRatio <- x$Subcellular..Channel.3..Num.single.spots/x$Subcellular..Channel.4..Num.single.spots
  x$normastroRatio <- rescale(x$Subcellular..Channel.3..Num.single.spots,to=c(1,2))/rescale(x$Subcellular..Channel.4..Num.single.spots,to=c(1,2))
  
  x
})



astro_data <- lapply(astro_data,function(x){
  x$ShuffledRatio <- sample(x$Subcellular..Channel.3..Num.single.spots)/sample(x$Subcellular..Channel.4..Num.single.spots)
  x
})

ratiobreaks <- c(0,.1,.25,.5,.75,1,Inf)
normratiobreaks <- c(0,.6,.75,1,1.5,2)
binnedData <- lapply(astro_data[tg_animals],function(x){
  x$binnedRatio <- cut(x$astroRatio,ratiobreaks, include.lowest=T)
  x$normbinnedRatio <- cut(x$normastroRatio,normratiobreaks, include.lowest=T)
  
  x$shuffledbinnedRatio <- cut(x$ShuffledRatio,ratiobreaks, include.lowest=T)
  x
  
})

x <- map_df(.x = seq(1:length(binnedData)),
            .f = function(x) bind_rows(binnedData[x]) %>%
              mutate(id = row_number(),
                     df = x)) %>% arrange(id) 


y <- pivot_longer(x,cols = c("binnedRatio","shuffledbinnedRatio"), names_to ="type")

ggplot(x,aes(x=binnedRatio))+
  geom_histogram(stat = "count")


ggplot(y) +
  geom_histogram(aes(x = value, fill = type),
                 alpha = 0.5,
                 position = "identity",
                 stat="count") +
  labs(title = "Overlaid Histograms",
       x = "Bins",
       y = "Count") +
  scale_fill_manual(values = c("binnedRatio" = "blue", "shuffledbinnedRatio" = "green"))

## Permutation test for spot relationships -------------------------

# Density of ch3 
ggplot(data.frame(rescale(x$Subcellular..Channel.2..Num.single.spots)), aes(x=rescale(x$Subcellular..Channel.2..Num.single.spots))) +
  geom_density()
ggplot(data.frame(rescale(x$Subcellular..Channel.2..Num.single.spots)), aes(x=rescale(x$Subcellular..Channel.2..Num.single.spots))) +
  geom_density()
ggplot(x, aes(x=Subcellular..Channel.4..Num.single.spots)) +
  geom_density()

library(surveillance)
ptest_spots <- permutationTest(as.numeric(x$Subcellular..Channel.3..Num.single.spots +21.8), (as.numeric(x$Subcellular..Channel.4..Num.single.spots)), nPermutation = 9999,
                               plot = TRUE, verbose = TRUE)

# try on normalized values
ptest_norm <- permutationTest(as.numeric(rescale(x$Subcellular..Channel.3..Num.single.spots)), (as.numeric(rescale(x$Subcellular..Channel.4..Num.single.spots))), nPermutation = 100000,
                              plot = TRUE, verbose = TRUE)
# Ratio random test ===================
normRatios <- (as.numeric(rescale(x$Subcellular..Channel.3..Num.single.spots,to = c(1,2)))) / (as.numeric(rescale(x$Subcellular..Channel.4..Num.single.spots,to = c(1,2))))


randRatios <- runif(n=length(normRatios),min = 0.5,max = 2)
ptest <- permutationTest(normRatios, randRatios, nPermutation = 1000,
                         plot = T, verbose = T)

perm <- 10000
n=1
diffs <- vector()
pvals <- vector()
while (n < perm) {
  randRatios <- runif(n=length(normRatios),min = 0.5,max = 2)
  ptest <- permutationTest(normRatios, randRatios, nPermutation = 100000,
                           plot = FALSE, verbose = FALSE)
  pvals <- c(pvals,ptest$pVal.permut)
  diffs <- c(diffs,ptest$diffObs)
  n=n+1
}


# Cell Count Distances


t.test(x=x$Count[1:4], y=x$Count[17:20],paired = FALSE)
t.test(x=x$Count[5:8], y=x$Count[21:24],paired = FALSE)
t.test(x=x$Count[9:12], y=x$Count[25:28],paired = FALSE)
t.test(x=x$Count[13:16], y=x$Count[29:32],paired = FALSE)

