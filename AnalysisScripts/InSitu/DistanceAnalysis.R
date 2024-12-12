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
dir <- file.path(getwd(),"Documents/InsituAnalysis/S100a6_results")
det_files <- list.files(file.path(dir,"detection_results"))                                                  # Extract file names from this file path
det_files 
det_names<-c()                                                                           # storage for names of datasets
det_data<-list()                                                                         # storage for datasets


for(i in 1:length(det_files)) {                                                          
  
  # (1) create dataset name from filename/datafile, store in variable det_temp
  
  det_names<-append(det_names, paste0(substring(det_files[i], 1, 2), sep = "_", substring(det_files[i], 20, 21))) # creates+stores names for dataset size:length(det_data)
  
  # (2.2) try read table for .txt files 
  
  det_data[[i]]<- as_tibble(read.table(file.path(dir,"detection_results",  
                                                 det_files[i]), sep = '\t', header = TRUE))                                        # stores each .txt data set 
  
} ## LOOP END



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


cutoff<-c(5000)
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
plaque_data <- list()
# runs for i= number of input datasets
for (dataset in 1:length(det_data)){     # runs for i= number of input datasets
  
  cells <-  filter(det_data[[dataset]], Object.type  == "Cell")  
  
  plaques <- filter(det_data[[dataset]], Classification  == "Plaques")  
  
  if("Nucleus..FITC.mean" %in% colnames(cells)){
    astros <-  filter(cells, Nucleus..FITC.mean > cutoff[j])
  } else{
    astros <-  filter(cells, Nucleus..FITC..C2..mean > cutoff[j])
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
  
  if(dim(plaques)[1]>1){
    plaque_coords <- data.frame(x=plaques$Centroid.X.µm,
                                y=plaques$Centroid.Y.µm,
                                size=0.1,color="blue",cell="plaque",
                                row.names = plaques$Object.ID)
    # plot map 
    stations <- rbind(cell_coords,astro_coords,plaque_coords)
    cell_maps[[det_names[dataset]]] <- ggplot(stations,aes(x=x,y=y))+ geom_point(aes(size=size,colour=color),show.legend = FALSE) +
      scale_radius(range = c(0, 1)) + 
      scale_colour_identity() + ggtitle(det_names[dataset])+theme_classic()+
      theme(plot.title = element_text(size = 20, face="bold",hjust = 0.5))
    # calculate all distances between plaques and astros
    dis <- rdist(astro_coords[,1:2],plaque_coords[,1:2]) # Euclidean distance
    dis_cells <- rdist(cell_coords[,1:2],plaque_coords[,1:2])
    
    # Add min distances back to astrocyte coords
    astro_coords$dis2nearestPlaque <- apply(dis,1,min)
    cell_coords$dis2nearestPlaque <- apply(dis_cells,1,min)
    
    # add in metadata for spot counts
    astro_coords <- cbind(astro_coords,astros)
    cell_coords <- cbind(cell_coords,cells)
    
    astro_data[[det_names[dataset]]] <- astro_coords
    cells_data[[det_names[dataset]]] <- cell_coords
    plaque_data[[det_names[dataset]]] <- plaques
    
    distance_plots_ch3[[det_names[dataset]]] <- ggplot(astro_coords,aes(x=dis2nearestPlaque,y=Subcellular..Channel.3..Num.single.spots)) +
      stat_summary_bin(bins = 100)
    distance_plots_ch4[[det_names[dataset]]] <- ggplot(astro_coords,aes(x=dis2nearestPlaque,y=Subcellular..Channel.4..Num.single.spots)) +
      stat_summary_bin(bins = 100)
    
  } else{
    stations <- rbind(cell_coords,astro_coords)
    cell_maps[[det_names[dataset]]] <- ggplot(stations,aes(x=x,y=y))+ geom_point(aes(size=size,colour=color),show.legend = FALSE) +
      scale_radius(range = c(0, 1)) + 
      scale_colour_identity() + ggtitle(det_names[dataset])+theme_classic()+
      theme(plot.title = element_text(size = 20, face="bold",hjust = 0.5))
    astro_coords <- cbind(astro_coords,astros)
    cell_coords <- cbind(cell_coords,cells)
    cells_data[[det_names[dataset]]] <- cell_coords
    astro_data[[det_names[dataset]]] <- astro_coords
    
  }
}




# Get spot avarages across samples -----------------------
breaks <- c(0,10,20,50,Inf)
tg_animals <- c("31_02", "31_03", "32_01", "32_02" ,"33_01","33_02")

binnedData <- lapply(astro_data[tg_animals],function(x){
  ag <- aggregate(x[,c(6,100,103,106)],by=list(distance=cut(x$dis2nearestPlaque,breaks, include.lowest=T),
                                               Cy3spots=x$Subcellular..Channel.3..Num.spots.estimated,
                                               Cy5spots=x$Subcellular..Channel.4..Num.spots.estimated),
                  FUN=mean)
  aggregate(ag[,-1], by=list(binnedDistance=ag$distance),FUN=mean)
  
})

x <- map_df(.x = seq(1:length(binnedData)),
            .f = function(x) bind_rows(binnedData[x]) %>%
              mutate(id = row_number(),
                     df = x)) %>%
  arrange(id) 

## Plot across samples -------------------------------------------
library(viridis)
png(file.path(dir,"DistancePlotApoe.png"),height = 500,width=500)
ggplot(x,aes(x=binnedDistance,y=Cy5spots)) +
  geom_boxplot(aes(fill=binnedDistance)) + geom_point(position=position_dodge(width=0.75))+
  scale_fill_viridis(discrete = TRUE, alpha=1,name="Distance",breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                     labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) + 
  scale_x_discrete(breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                   labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) +
  xlab("Distance to nearest plaque center") + ylab("Avg. DAA2 (Apoe) spots in astrocyte") +
  ggtitle("Average Apoe expression in astrocytes") +
  theme(plot.title = element_text(hjust=0.5,face="bold",size=20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

dev.off()

png(file.path(dir,"DistancePlotS100a6.png"),height = 500,width=500)
ggplot(x,aes(x=binnedDistance,y=Cy3spots)) +
  geom_boxplot(aes(fill=binnedDistance)) + geom_point(position=position_dodge(width=0.75))+
  scale_fill_viridis(discrete = TRUE, alpha=1,name="Distance",breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                     labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) + 
  scale_x_discrete(breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                   labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) +
  xlab("Distance to nearest plaque center") + ylab("Avg. DAA1(S100a6) spots in astrocyte")+
  ggtitle("Average S100a6 expression in astrocytes")+
  theme(plot.title = element_text(hjust=0.5,face="bold",size=20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))+ scale_color_manual(
          breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
          labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) 
dev.off()
#  Region specificity -------------

breaks <- c(0,10,20,50,Inf)

binnedData_region <- lapply(astro_data[tg_animals],function(x){
  ag <- aggregate(x[,c(6,100,103,106)],by=list(distance=cut(x$dis2nearestPlaque,breaks, include.lowest=T),
                                               Cy3spots=x$Subcellular..Channel.3..Num.spots.estimated,
                                               Cy5spots=x$Subcellular..Channel.4..Num.spots.estimated,
                                               Region=x$Parent),
                  FUN=mean)
  aggregate(ag[,-c(1,4)], by=list(Distance=ag$distance,
                                  Region=ag$Region),FUN=mean)
  
})

region_data <- map_df(.x = seq(1:length(binnedData_region)),
                      .f = function(x) bind_rows(binnedData_region[x]) %>%
                        mutate(id = row_number(),
                               df = x)) %>% arrange(id) 

png(file.path(dir,"DistancePlotS100a6_region.png"),height = 500,width=500)
ggplot(region_data,aes(x=Region,y=Cy3spots)) +
  geom_boxplot(aes(fill=Distance)) + 
  geom_point(aes(color=Distance),position=position_dodge(width=0.75),show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, alpha=1,breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                     labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) +
  xlab("Distance to nearest plaque center") + ylab("Avg. DAA1(S100a6) spots in astrocyte")+
  ggtitle("Average S100a6 expression in astrocytes")+
  theme(plot.title = element_text(hjust=0.5,face="bold",size=20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_color_manual(values =setNames(viridis(4),c("[0,10]","(10,20]","(20,50]","(50,Inf]")) )
dev.off()

png(file.path(dir,"DistancePlotApoe_region.png"),height = 500,width=500)
ggplot(region_data,aes(x=Region,y=Cy5spots)) +
  geom_boxplot(aes(fill=Distance)) + 
  geom_point(aes(color=Distance),position=position_dodge(width=0.75),show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, alpha=1,breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                     labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) +
  xlab("Distance to nearest plaque center") + ylab("Avg. DAA2(Apoe) spots in astrocyte")+
  ggtitle("Average Apoe expression in astrocytes")+
  theme(plot.title = element_text(hjust=0.5,face="bold",size=20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_color_manual(values =setNames(viridis(4),c("[0,10]","(10,20]","(20,50]","(50,Inf]")) )
dev.off()

save(list=c("cell_maps","distance_plots_ch3","distance_plots_ch4",
            "astro_data","cells_data","binnedData","binnedData_region","plaque_data"),
     file=file.path(dir,"data.files.Rdata"))

#save.image(file.path(dir,"data.files.Rdata")) 

load(file.path(dir,"data.files.Rdata"))

# Density plot across cells =========================================
library(ggplot2)

density_data <- lapply(astro_data, function(x){
  df <- data.frame(Name      = x$Name,
                   Region    = x$Parent,
                   DAA1Spots = x$Subcellular..Channel.3..Num.spots.estimated,
                   DAA2Spots = x$Subcellular..Channel.4..Num.spots.estimated,
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
title=text_grob("S100a6 Results", size = 20, face = "bold") 

png(file.path(dir,"DensityPlots_AllSamples.png"),width = 800,height=500)
do.call("grid.arrange", c(density_plots, list(ncol=2,top=title)))
dev.off()
## Region ==========
density_plots_region <- lapply(density_plots, function(x){
  x + facet_wrap(~Region)
})

png(file.path(dir,"DensityPlots_Region_AllSamples.png"),width = 1000,height=500)
do.call("grid.arrange", c(density_plots_region, list(nrow=2,top=title)))
dev.off()

## Proportions ==========

pp <- lapply(astro_data[tg_animals], function(x){
  df <- data.frame(Name      = x$Name,
                   Region    = x$Parent,
                   DAA1Spots = x$Subcellular..Channel.3..Num.spots.estimated,
                   DAA2Spots = x$Subcellular..Channel.4..Num.spots.estimated,
                   PlaqueDis = x$dis2nearestPlaque,
                   SpotPro  = x$Subcellular..Channel.3..Num.spots.estimated/x$Subcellular..Channel.4..Num.spots.estimated)
  df
})

data <- map_df(.x = seq(1:length(pp)),
               .f = function(x) bind_rows(pp) %>%
                 mutate(id = row_number(),
                        df = x)) %>% arrange(id) 

ggplot(data=data, aes(x=SpotPro)) +
  geom_density(adjust=1.5, alpha=.4) + ggtitle("Aggregated S100a6 Proportions") + 
  theme(plot.title = element_text(size=20,face="bold",hjust=0.5))

# look at aggregated across animals =========================
density_data_long <- map_df(.x = seq(1:length(density_data)),
                            .f = function(x) bind_rows(density_data) %>%
                              mutate(id = row_number(),
                                     df = x)) %>% arrange(id) 

png(file.path(dir,"DensityPlots_Aggregated.png"),width = 800,height=500)
ggplot(data=density_data_long, aes(x=estimatedSpots, group=channel, fill=channel)) +
  geom_density(adjust=1.5, alpha=.4) + ggtitle("Aggregated S100a6 Data") + 
  theme(plot.title = element_text(size=20,face="bold",hjust=0.5))
dev.off()

png(file.path(dir,"DensityPlots_AggregatedRegion.png"),width = 800,height=500)
ggplot(data=density_data_long, aes(x=estimatedSpots, group=channel, fill=channel)) +
  geom_density(adjust=1.5, alpha=.4) + ggtitle("Aggregated S100a6 Data") + 
  theme(plot.title = element_text(size=20,face="bold",hjust=0.5)) +
  facet_wrap(~Region)
dev.off()

# Scatter plot of cells with DAA1 and DAA2 spots ==============

scatters <- list()

scatters <- lapply(astro_data, function(x){
  ggscatter(x, x = "Subcellular..Channel.3..Num.spots.estimated",
            y = "Subcellular..Channel.4..Num.spots.estimated", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Estimated S1006 spots",
            ylab = "Estimated Apoe spots")
  
})

do.call("grid.arrange", c(scatters, list(ncol=3,top=title)))

## collapse all astrocytes
library(dplyr)
data <- bind_rows(astro_data, .id = "column_label")
ggscatter(data, x = "Subcellular..Channel.3..Num.spots.estimated",
          y = "Subcellular..Channel.4..Num.spots.estimated", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Estimated S1006 spots",
          ylab = "Estimated Apoe spots")
## Boxplots of dispersions  ========
ggplot(data.frame(data$Subcellular..Channel.3..Num.spots.estimated),aes(y=data.Subcellular..Channel.3..Num.spots.estimated)) +
  geom_boxplot() +ggtitle("S100a6") + ylab("# spots") + theme_bw() + 
  theme(plot.title = element_text(hjust=0.5,size=22,face = "bold"))
ggplot(data.frame(data$Subcellular..Channel.4..Num.spots.estimated),aes(y=data.Subcellular..Channel.4..Num.spots.estimated)) +
  geom_boxplot() +ggtitle("Apoe") + ylab("# spots") + theme_bw() + 
  theme(plot.title = element_text(hjust=0.5,size=22,face = "bold"))

#cutoffs determined based on third quantile
cutoffc3 <- 15.9
cutoffc4 <- 56.61

## Label Astro types =================
astro_data <- lapply(astro_data, function(x){
  x$astroType <- "DN"
  x$astroType[x$Subcellular..Channel.3..Num.spots.estimated > cutoffc3+1 &
                x$Subcellular..Channel.4..Num.spots.estimated < cutoffc4 ] <- "DAA1"
  x$astroType[x$Subcellular..Channel.3..Num.spots.estimated < cutoffc3 &
                x$Subcellular..Channel.4..Num.spots.estimated > cutoffc4+1 ] <- "DAA2"
  x$astroType[x$Subcellular..Channel.3..Num.spots.estimated > cutoffc3 &
                x$Subcellular..Channel.4..Num.spots.estimated > cutoffc4 ] <- "DP"
  x
})
data <- bind_rows(astro_data, .id = "column_label")



## colored scatter ===============
ggscatter(data, x = "Subcellular..Channel.3..Num.spots.estimated",
          y = "Subcellular..Channel.4..Num.spots.estimated", 
          color = "astroType",
          palette =c("goldenrod1","red","black","dodgerblue"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Estimated S1006 spots",
          ylab = "Estimated Apoe spots")

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


# Get spot avarages across samples -----------------------
breaks <- c(0,10,20,50,Inf)

binnedData <- lapply(astro_data,function(x){
  test <- x %>%
    group_by(dis2nearestPlaque, astroType) %>%
    summarise(Count = n())
  aggregate(test,by=list(distance=cut(test$dis2nearestPlaque,breaks, include.lowest=T),
                         Subtype=test$astroType),FUN=count)
  
  
})


x <- map_df(.x = seq(1:length(binnedData)),
            .f = function(x) bind_rows(binnedData[x]) %>%
              mutate(id = row_number(),
                     df = x)) %>%
  arrange(id) 

ggplot(x,aes(x=distance,y=Count)) +
  geom_boxplot(aes(fill=distance)) + geom_point(position=position_dodge(width=0.75))+
  scale_fill_viridis(discrete = TRUE, alpha=1,name="Distance",breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                     labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) + 
  scale_x_discrete(breaks=c("[0,10]","(10,20]","(20,50]","(50,Inf]"),
                   labels=c("0-10 μm","10-20 μm","20-50 μm","50+ μm")) + facet_wrap(~Subtype)+
  xlab("Distance to nearest plaque center") + ylab("# of Cells") +
  ggtitle("# of astrocytes in relation to pathology") +
  theme(plot.title = element_text(hjust=0.5,face="bold",size=20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Norm Ratio labeling ================
astro_data <- lapply(astro_data, function(x){
  x$normch3 <- x$Subcellular..Channel.3..Num.spots.estimated / x$Nucleus..Perimeter
  x$normch4 <- x$Subcellular..Channel.4..Num.spots.estimated / x$Nucleus..Perimeter
  
  x$normastroType <- "DN"
  x$normastroType[x$normch3 >= 1.2  ] <- "DAA1"
  x$normastroType[x$normch3 >= 1.2 &
                    x$normch4 <=1.5] <- "DAA1"
  x$normastroType[x$normch3 < 0.5 &
                    x$normch4 >= 1.5 ] <- "DAA2"
  x
})

data <- bind_rows(astro_data, .id = "column_label")
ggscatter(data, x = "normch3",
          y = "normch4",color="normastroType",
          add = "reg.line", conf.int = TRUE,
          palette =c("goldenrod1","red","black"),
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Estimated S100a6 spots",
          ylab = "Estimated Apoe spots")

binnedData <- lapply(astro_data[tg_animals],function(x){
  test <- x %>%
    group_by(dis2nearestPlaque, normastroType) %>%
    summarise(Count = n())
  aggregate(test[,-c(1,2)], by = list(distance = cut(test$dis2nearestPlaque,breaks, include.lowest=T),
                                      Subtype = test$normastroType ), FUN = sum )
  
  
})


x <- map_df(.x = seq(1:length(binnedData)),
            .f = function(x) bind_rows(binnedData[x]) %>%
              mutate(id = row_number(),
                     df = x)) %>%
  arrange(id)
means <- c(mean(x$Count[1:6]),mean(x$Count[7:12]),mean(x$Count[13:18]),
           mean(x$Count[19:24]),mean(x$Count[25:30]),mean(x$Count[31:36]),
           mean(x$Count[37:42]),mean(x$Count[43:48]),mean(x$Count[49:54]),
           mean(x$Count[55:60]),mean(x$Count[61:66]),mean(x$Count[67:72]))

sd <-  c(sd(x$Count[1:6]),sd(x$Count[7:12]),sd(x$Count[13:18]),
         sd(x$Count[19:24]),sd(x$Count[25:30]),sd(x$Count[31:36]),
         sd(x$Count[37:42]),sd(x$Count[43:48]),sd(x$Count[49:54]),
         sd(x$Count[55:60]),sd(x$Count[61:66]),sd(x$Count[67:72]))


df2 <- data.frame(Distance = rep(levels(factor(x$distance)),3),
                  Subtype = c(rep("DAA1",4),rep("DAA2",4),rep("DN",4)),
                  Count = means,
                  SD = sd)
pd <- df2[c(1:8),]
pd$Distance <- factor(pd$Distance,
                      levels =c("[0,10]","(10,20]", "(20,50]","(50,Inf]"))

ggplot(pd,aes(x=Distance,y=Count,group=Subtype,color=Subtype)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin=Count-SD, ymax=Count+SD), width=.2,
                position=position_dodge(0.05))

# Ratio of total cells ================

means <- c(mean(x$Count[1:6]/sum(x$Count[1:24])),mean(x$Count[7:12]/sum(x$Count[1:24])),
           mean(x$Count[13:18]/sum(x$Count[1:24])),
           mean(x$Count[19:24]/sum(x$Count[1:24])),mean(x$Count[25:30]/sum(x$Count[25:48])),
           mean(x$Count[31:36]/sum(x$Count[25:48])),
           mean(x$Count[37:42]/sum(x$Count[25:48])),mean(x$Count[43:48]/sum(x$Count[25:48])),
           mean(x$Count[49:54]/sum(x$Count[49:72])),
           mean(x$Count[55:60]/sum(x$Count[49:72])),mean(x$Count[61:66]/sum(x$Count[49:72]))
           ,mean(x$Count[67:72]/sum(x$Count[49:72])))

sd <- c(sd(x$Count[1:6]/sum(x$Count[1:24])),sd(x$Count[7:12]/sum(x$Count[1:24])),
        sd(x$Count[13:18]/sum(x$Count[1:24])),
        sd(x$Count[19:24]/sum(x$Count[1:24])),sd(x$Count[25:30]/sum(x$Count[25:48])),
        sd(x$Count[31:36]/sum(x$Count[25:48])),
        sd(x$Count[37:42]/sum(x$Count[25:48])),sd(x$Count[43:48]/sum(x$Count[25:48])),
        sd(x$Count[49:54]/sum(x$Count[49:72])),
        sd(x$Count[55:60]/sum(x$Count[49:72])),sd(x$Count[61:66]/sum(x$Count[49:72]))
        ,sd(x$Count[67:72]/sum(x$Count[49:72])))

rt <- data.frame(Distance = rep(levels(factor(x$distance)),3),
                 Subtype = c(rep("DAA1",4),rep("DAA2",4),rep("DN",4)),
                 Count = means,
                 SD = sd)
pd <- rt[c(1:8),]
pd$Distance <- factor(pd$Distance,
                      levels =c("[0,10]","(10,20]", "(20,50]","(50,Inf]"))
color.codes <- c("goldenrod1","red")
zone <- levels(factor(pd$Subtype))
ggplot(pd,aes(x=Distance,y=Count,group=Subtype,color=Subtype)) +
  geom_line() + geom_point(size=3) +
  geom_errorbar(aes(ymin=Count-SD, ymax=Count+SD), width=.5,
                position=position_dodge(0.05)) + theme_classic() +
  ylab("Ratio of Astrocyte type") + xlab("Distance from Pathology") +
  scale_colour_manual(values=setNames(color.codes, zone)) + 
  theme(axis.title = element_text(size=24, face = "bold"),
        axis.text = element_text(size=14))

t.test(x=x$Count[1:6], y=x$Count[25:30],paired = FALSE)
t.test(x=x$Count[7:12], y=x$Count[31:36],paired = FALSE)
t.test(x=x$Count[13:18], y=x$Count[37:42],paired = FALSE)
t.test(x=x$Count[19:24], y=x$Count[43:48],paired = FALSE)

sub <- subset(readRDS("~/Documents/AstrocytePaperAstrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"),
              subset= StudyName=="Lee TauPS2APP")

coldata <- as.data.frame(sub@meta.data)
study_cluster_tmp_by_study <- coldata %>%
  group_by(finalClusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 
sc <- rbind((study_cluster_tmp_by_study[1,]+study_cluster_tmp_by_study[4,]),
            study_cluster_tmp_by_study[2,],study_cluster_tmp_by_study[3,])
sc$dataType <- "sCell"
sc[1,1] <- "homeo_synapse"  
sc[2,1] <- "DAA1" 
sc[3,1] <- "DAA2" 
colnames(sc)[1]<- "astroType"
insitu <- bind_rows(astro_data[tg_animals], .id = "column_label")
insitu <- insitu %>%
  group_by(normastroType) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
insitu$dataType <- "in situ"
insitu[3,1] <- "homeo_synapse"  
colnames(insitu)[1]<- "astroType"

plot <- rbind(sc,insitu)

color.codes <- c("goldenrod1","red","dodgerblue")
zone <- levels(factor(plot$astroType))
ggplot(plot, aes(fill=astroType, y=freq,x=dataType)) +
  geom_bar(position="stack", stat="identity") + theme_void()+ theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14)) +
  xlab("Clusters") + ylab("Proportion of cells within a study") +
  theme(axis.title.y =  element_text(angle = 90))+
  scale_fill_manual(values=setNames(color.codes, zone))

# Control vs Transgenic expression -------------------------------
tg_animals <- c("31_02", "31_03", "32_01", "32_02" ,"33_01","33_02")

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
colnames(expressionMeans) <- c("Sample","Genotype","S100a6","Apoe")

expressionMeans$S100a6 <- as.numeric(expressionMeans$S100a6)
expressionMeans$Apoe <- as.numeric(expressionMeans$Apoe)
expressionMeans$S100a6_norm <- expressionMeans$S100a6
expressionMeans$Apoe_norm <- expressionMeans$Apoe

expressionMeans$S100a6_norm[which(expressionMeans$Genotype=="TauPS2APP")] <- expressionMeans$S100a6[which(expressionMeans$Genotype=="TauPS2APP")]/PLC$numPlaques
expressionMeans$Apoe_norm[which(expressionMeans$Genotype=="TauPS2APP")] <- expressionMeans$Apoe[which(expressionMeans$Genotype=="TauPS2APP")]/PLC$numPlaques
expAvgbyAnimal <- data.frame(Sample = unique(gsub("_.*","",expressionMeans$Sample)),
                             Genotype = c(rep("NonTG",3),rep("TauPS2APP",3)),
                             S100a6 = .colMeans(expressionMeans$S100a6, 2, length(expressionMeans$S100a6) / 2),
                             Apoe  = .colMeans(expressionMeans$Apoe, 2, length(expressionMeans$Apoe) / 2))

expressionMeans <- pivot_longer(expressionMeans,
                                cols=c(S100a6_norm,Apoe_norm),
                                names_to = "ChannelName")
expAvgbyAnimal <- pivot_longer(expAvgbyAnimal,
                                cols=c(S100a6,Apoe),
                                names_to = "ChannelName")
sub <- expressionMeans[-which(startsWith(expressionMeans$Sample,"33")),]


ggplot(expressionMeans, aes(x=Genotype, y=value,fill=Genotype)) + theme_classic() +
  geom_boxplot(outlier.shape=NA) + scale_fill_viridis(discrete = TRUE, alpha=1) +
  geom_jitter(width = 0.4,size=3,aes(color=Sample))  +
  facet_wrap(~factor(ChannelName),scales = "free", ncol=2) + xlab("Disease Label") +
  theme(axis.text.x = element_text(angle=75,vjust = 0.5),axis.text.y = element_text(size=20)) +
  ylab("Mean Expression") +  
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),legend.text = element_text(size=14),
        legend.title = element_text(size=16)) 



ggplot(expAvgbyAnimal, aes(x=Genotype, y=value,fill=Genotype)) + theme_classic() +
  geom_boxplot(outlier.shape=NA) + scale_fill_viridis(discrete = TRUE, alpha=1) +
  geom_jitter(width = 0.4,size=3,aes(color=Sample))  +
  facet_wrap(~factor(ChannelName),scales = "free", ncol=2) + xlab("Disease Label") +
  theme(axis.text.x = element_text(angle=75,vjust = 0.5),axis.text.y = element_text(size=20)) +
  ylab("Mean Expression") +  
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),legend.text = element_text(size=14),
        legend.title = element_text(size=16)) 

# Plaque load correaltion -------------------------------------------
numPlaques <- vector()
for(data in 1:length(plaque_data)){
  numPlaques <- c(numPlaques,dim(plaque_data[[data]])[1])
}
PLC <- data.frame(Sample = expressionMeans$Sample[which(expressionMeans$Genotype=="TauPS2APP" &
                                                       expressionMeans$ChannelName == "S100a6") ],
                  S100a6 = expressionMeans$value[which(expressionMeans$Genotype=="TauPS2APP" &
                                                    expressionMeans$ChannelName == "S100a6") ],
                  Apoe = expressionMeans$value[which(expressionMeans$Genotype=="TauPS2APP" &
                                                     expressionMeans$ChannelName == "Apoe") ],
                  numPlaques = numPlaques )

PLC$numPlaques[6] <- 6683
ggscatter(PLC, x = "numPlaques",
          y = "S100a6", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Estimated number of Plaques",
          ylab = "Estimated S100a6 norm Expression")
ggscatter(PLC, x = "numPlaques",
          y = "Apoe", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Estimated number of Plaques",
          ylab = "Estimated Apoe norm Expression")
ggplot(PLC,aes(x=Sample,y=numPlaques)) +
  geom_bar(stat = "count")

barplot(numPlaques ~ Sample,data=PLC[c()])

#plaque distance
data <- astro_data[tg_animals]
data <- bind_rows(data, .id = "column_label")

ggscatter(data, x = "normch3",
          y = "dis2nearestPlaque", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Estimated S1006 spots",
          ylab = "Distance to nearest plaque spots")
