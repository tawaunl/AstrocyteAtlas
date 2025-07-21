# Extract metadata
metadata <- data@meta.data %>%
  mutate(Cluster = data$ClusterNames)  # Add cluster (identity) information

# Count the number of cells per BrainRegion, Sample, and Cluster
cell_counts <- metadata %>%
  group_by(Cluster, BrainRegion, Sample) %>%
  summarize(CellCount = n(), .groups = "drop")

# Calculate the proportion of cells for each Sample within each Cluster and BrainRegion
proportion_data <- cell_counts %>%
  group_by(Cluster, Sample) %>%
  mutate(Proportion = CellCount / sum(CellCount))  # Normalize to get proportions per sample

# Create the bar plot
ggplot(proportion_data, aes(x = BrainRegion, y = Proportion, fill = BrainRegion)) +
  geom_bar(stat = "identity", position = "stack") +# Stacked bar chart for sample contributions
  facet_wrap(~ Cluster, scales = "free_y") +  # Format y-axis as percentages
  labs(title = "Percentage of Cells in Brain Regions by Sample and Cluster",
       x = "Brain Region",
       y = "Percentage of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text =element_text(size = 16),
        axis.title = element_text(size=18),
        plot.title = element_text(size=20,hjust = 0.5),
        legend.position = "none") +  # Rotate x-axis labels for readability
  scale_fill_brewer(palette = "Set2")  # Prettier color palette



#Agt expression ----------------
plot_data <- FetchData(data, vars = c("Agt", "BrainRegion")) %>%
  mutate(Cluster = data$ClusterNames)  # Add cluster IDs to data

# Calculate mean expression of "Agt" within each BrainRegion and Cluster
mean_expression <- plot_data %>%
  group_by(BrainRegion, Cluster) %>%
  summarize(MeanExpression = mean(Agt), .groups = "drop")  # Calculate mean expression per group

# Create the plot
ggplot(plot_data, aes(x = BrainRegion, y = Agt)) +
  # Add individual cell expression as points (jittered)
  geom_jitter(aes(color = BrainRegion), width =.1, alpha =0.9, size = 1) +
  # Add bar plot for mean expression
  geom_bar(data = mean_expression, 
           aes(x = BrainRegion, y = MeanExpression, fill = BrainRegion),
           stat = "identity", alpha = .7, width = 0.7) +
  # Adjustments for facets by cluster
  facet_wrap(~ Cluster, scales = "free_y") +
  # Add labels and themes
  labs(title = "Expression of 'Agt' by Brain Region and Cluster",
       x = "Brain Region",
       y = "Expression of 'Agt'") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text =element_text(size = 16),
        axis.title = element_text(size=18),
        plot.title = element_text(size=20,hjust = 0.5),
        legend.position = "none",
        strip.text = element_text(size=14)) +
  scale_fill_brewer(palette = "Set2") +  # Adjust bar colors
  scale_color_brewer(palette = "Set2")   # Adjust dot colors

#Ogt expression -----------------------
# Extract information: normalized expression of "Agt", BrainRegion, and clusters
plot_data <- FetchData(data, vars = c("Ogt", "BrainRegion")) %>%
  mutate(Cluster = data$ClusterNames)  # Add cluster IDs to data

# Calculate mean expression of "Agt" within each BrainRegion and Cluster
mean_expression <- plot_data %>%
  group_by(BrainRegion, Cluster) %>%
  summarize(MeanExpression = mean(Ogt), .groups = "drop")  # Calculate mean expression per group

# Create the plot
ggplot(plot_data, aes(x = BrainRegion, y = Ogt)) +
  # Add individual cell expression as points (jittered)
  geom_jitter(aes(color = BrainRegion), width =.1, alpha =0.9, size = 1) +
  # Add bar plot for mean expression
  geom_bar(data = mean_expression, 
           aes(x = BrainRegion, y = MeanExpression, fill = BrainRegion),
           stat = "identity", alpha = .7, width = 0.7) +
  # Adjustments for facets by cluster
  facet_wrap(~ Cluster, scales = "free_y") +
  # Add labels and themes
  labs(title = "Expression of 'Ogt' by Brain Region and Cluster",
       x = "Brain Region",
       y = "Expression of 'Ogt'") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text =element_text(size = 16),
        axis.title = element_text(size=18),
        plot.title = element_text(size=20,hjust = 0.5),
        legend.position = "none",
        strip.text = element_text(size=14)) +
  scale_fill_brewer(palette = "Set2") +  # Adjust bar colors
  scale_color_brewer(palette = "Set2")   # Adjust dot colors

# Proportion by BrainRegion and Cluster

#C. Study Bar Plot ==============
coldata <- as.data.frame(data@meta.data)
study_cluster_tmp_by_study <- coldata %>%
  group_by(Model, finalClusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

cluster1 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==1),]
cluster2 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==2),]
cluster3 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==3),]
cluster4 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==4),]

cluster1$scaled <- cluster1$freq/sum(cluster1$freq)
cluster2$scaled <- cluster2$freq/sum(cluster2$freq)
cluster3$scaled <- cluster3$freq/sum(cluster3$freq)
cluster4$scaled <- cluster4$freq/sum(cluster4$freq)


study_cluster_tmp_by_study_scaled <- rbind(cluster1,cluster2,cluster3,cluster4)
#### Plot C -----------

ggplot(study_cluster_tmp_by_study_scaled, aes(fill=Model, y=scaled, x=finalClusters)) +
  geom_bar(position="stack", stat="identity", alpha=0.9) + theme_void()+ theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14)) +
  xlab("Clusters") + ylab("Proportion of cells across all studies") +
  theme(axis.title.y =  element_text(angle = 90))

# Brain Region Stacked ------------------
coldata <- as.data.frame(data@meta.data)
study_cluster_tmp_by_study <- coldata %>%
  group_by(BrainRegion, finalClusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

cluster1 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==1),]
cluster2 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==2),]
cluster3 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==3),]
cluster4 <- study_cluster_tmp_by_study[which(study_cluster_tmp_by_study$finalClusters==4),]

cluster1$scaled <- cluster1$freq/sum(cluster1$freq)
cluster2$scaled <- cluster2$freq/sum(cluster2$freq)
cluster3$scaled <- cluster3$freq/sum(cluster3$freq)
cluster4$scaled <- cluster4$freq/sum(cluster4$freq)


study_cluster_tmp_by_study_scaled <- rbind(cluster1,cluster2,cluster3,cluster4)
#### Plot C -----------

ggplot(study_cluster_tmp_by_study_scaled, aes(fill=BrainRegion, y=scaled, x=finalClusters)) +
  geom_bar(position="stack", stat="identity", alpha=0.9) + theme_void()+ theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) + theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 14)) + theme(axis.text.y = element_text(size = 14)) +
  xlab("Clusters") + ylab("Proportion of cells within a Cluster") +
  theme(axis.title.y =  element_text(angle = 90))
library(tidyr)

# Cellularity Brain Region ----------
normalized_data <- metadata %>%
  group_by(Sample, BrainRegion, ClusterNames) %>% 
  summarise(CellCount = n(), .groups = "drop") %>%
  
  # (2) Normalize: Divide cell counts by the total number of cells in each sample
  group_by(Sample) %>% 
  mutate(TotalCellsPerSample = sum(CellCount)) %>% 
  mutate(Percentage = CellCount / TotalCellsPerSample *100) %>%
  ungroup()



# () Create the box plot
ggplot(normalized_data, aes(x = BrainRegion, y = Percentage, fill = BrainRegion)) +
  geom_boxplot(outlier.shape = NA, alpha = .5) + # Boxplot without outliers
  geom_jitter(width = 0.2, size =1 ) + # Jitter to show sample-specific values
  facet_wrap(~ClusterNames, scales = "free_y",ncol=2) + # Facet by cluster
  labs(
    title = "Cluster Distribution by Brain Region",
    x = "Brain Region",
    y = "Percentage of Cells in Regaion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # Rotate x-axis labels if needed
    strip.text = element_text(size =13 , face = "bold")# Format facet labels
  )

