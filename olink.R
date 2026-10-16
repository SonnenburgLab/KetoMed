library(vegan)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(nlme)
library(pheatmap)
library(reshape2)

# Read olink data
olink_data <- read.csv("keto_med_olink_deltas.csv")
olink_data_filtered <- read.csv("keto_med_olink_filtered.csv")
olink_data_filtered2 <- read.csv("keto_med_olink_filtered2.csv")

# Separate metadata and immune markers
metadata <- olink_data[, 1:6]  # Adjust the column indices based on meta data columns
immune_markers <- olink_data[, 7:84]  # Adjust to include all immune marker columns

metadata_filtered2 <- olink_data_filtered2[, 1:4]  # Adjust the column indices based on meta data columns
immune_markers_filtered2 <- olink_data_filtered2[, 5:82]  # Adjust to include all immune marker columns

# Scale the immune markers
scaled_markers_filtered2 <- scale(immune_markers_filtered2)
# Set PPID as row names
rownames(scaled_markers_filtered2) <- olink_data_filtered2$sample
# Add metadata to scaled filtered immune markers
immune_combined2 <- cbind(data.frame(metadata_filtered2, scaled_markers_filtered2))

###### Paper Figures #######
## Generate box plots of signficant different olink markers ##
custom_colors <- c("bl" = "#b7b7b7ff","keto_del" = "#c7885fff", "keto_own" = "#f5c09fff", "med_del" = "#669289ff", "med_own" = "#afd3caff")

## Generate box plots of signficant different olink markers ##
custom_colors2 <- c("bl" = "#b7b7b7ff","keto" = "#c7885fff", "med" = "#669289ff")

# Exclude all samples except baseline and keto_del and Keto_own
immune_blketo <- subset(immune_combined, diet.phase %in% c("bl", "keto_del", "keto_own"))

# Exclude all samples except baseline and keto_del and Keto_own med_del and med_own
immune_combined <- subset(immune_combined, diet.phase %in% c("bl", "keto_del", "keto_own", "med_del", "med_own"))

# test for signficant changes within keto diet immune markers #


# Signficant markers from paired t test in excel
# ttest paired	bl v keto_del	rank	adj p
#il18_inflam
#cst5	
#ccl25	
#fgf21_inflam	

#plot signficant markers
BKD_cst5 <- immune_combined2[, c("diet.phase","cst5")]
BKD_il18 <- immune_combined2[, c("diet.phase","il18_inflam")]
BKD_ccl25 <- immune_combined2[, c("diet.phase","ccl25")]
BKD_fgf21 <- immune_combined2[, c("diet.phase","fgf21_inflam")]


ggplot(BKD_fgf21, aes(x = diet.phase, y = fgf21_inflam, fill = diet.phase)) + 
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.7) + 
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.8, color = "black") + 
  scale_fill_manual(values = custom_colors2) + 
  labs( x = NULL, y = "fgf21_inflam") + 
  theme_classic(base_size = 12)

#plot signficant markers all diet phases
BKD_cst5 <- immune_combined[, c("diet.phase","cst5")]
BKD_ccl11 <- immune_combined[, c("diet.phase","ccl11")]
BKD_il12b <- immune_combined[, c("diet.phase","il12b")]


ggplot(BKD_cst5, aes(x = diet.phase, y = cst5, fill = diet.phase)) + 
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.7) + 
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8, color = "black") + 
  scale_fill_manual(values = custom_colors) + 
  labs( x = NULL, y = "cst5") + 
  theme_classic(base_size = 12)





#set colors
# Define the number of colors for negative and positive values
num_neg_colors <- 50  # Number of colors for negative values
num_pos_colors <- 50  # Number of colors for positive values
# Create the custom color palette
custom_colors <- c(colorRampPalette(c("blue", "white"))(num_neg_colors), 
                    "white", 
                    colorRampPalette(c("white", "red"))(num_pos_colors))


# Generate the heatmap
pheatmap(scaled_markers, 
         main = "olink changes",
         #color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE)

# Combine metadata with scaled immune markers
scaled_data <- cbind(metadata, scaled_markers)

#subset data
subset_data <- scaled_data[scaled_data$diet.phase == "med_del", ]
immune_data_for_heatmap <- subset_data[, -c(1:6)]

# Generate the heatmap
pheatmap(immune_data_for_heatmap, 
         main = "olink changes",
         color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE)

##### Heatmap of significant olink values ##########
# Read olink data of only complete subject dataset
olink_data_complete <- read.csv("km_olink_deltas_sigave.csv")

# Separate metadata and immune markers
metadata_complete <- olink_data_complete[, 1:6]  # Adjust the column indices based on meta data columns
immune_markers_complete <- olink_data_complete[, 7:23]  # Adjust to include all immune marker columns

# Scale the immune markers
#scaled_markers_complete <- scale(immune_markers_complete)
# Set PPID as row names
#rownames(scaled_markers_complete) <- olink_data_complete$sample

# Combine metadata with scaled immune markers
data_complete <- cbind(metadata_complete, immune_markers_complete)

#set colors
# Define the number of colors for negative and positive values
num_neg_colors <- 50  # Number of colors for negative values
num_pos_colors <- 50  # Number of colors for positive values
# Create the custom color palette
custom_colors <- c(colorRampPalette(c("blue", "white"))(num_neg_colors), 
                   "white", 
                   colorRampPalette(c("white", "red"))(num_pos_colors))

num_colors <- 100  # the same number for both negative and positive can be adjusted

# Create a custom color palette
custom_colors <- c(
  colorRampPalette(c("blue", "white"))(num_colors),  # for values from 0 to -1
  "white",                                            # for zero
  colorRampPalette(c("white", "red"))(num_colors)   # for values from 0 to 1
)


# Generate the heatmap
pheatmap(immune_markers_complete, 
         color = my_palette, 
         cluster_rows = FALSE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)

# Define color palette (Blue-White-Red)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 50)

# Plot with custom breaks to center at 0
pheatmap(immune_markers_complete, 
         color = my_palette, 
         cluster_rows = FALSE, 
         breaks = seq(-max(abs(immune_markers_complete)), max(abs(immune_markers_complete)), length.out = 51))



##### Complete dataset Analysis ##########
# Read olink data of only complete subject dataset
olink_data_complete <- read.csv("km_olink_deltas_complete.csv")

# Separate metadata and immune markers
metadata_complete <- olink_data_complete[, 1:6]  # Adjust the column indices based on meta data columns
immune_markers_complete <- olink_data_complete[, 7:84]  # Adjust to include all immune marker columns

# Scale the immune markers
scaled_markers_complete <- scale(immune_markers_complete)
# Set PPID as row names
rownames(scaled_markers_complete) <- olink_data_complete$sample

# Combine metadata with scaled immune markers
scaled_data_complete <- cbind(metadata_complete, scaled_markers_complete)

#set colors
# Define the number of colors for negative and positive values
num_neg_colors <- 50  # Number of colors for negative values
num_pos_colors <- 50  # Number of colors for positive values
# Create the custom color palette
custom_colors <- c(colorRampPalette(c("red", "white"))(num_neg_colors), 
                   "white", 
                   colorRampPalette(c("white", "blue"))(num_pos_colors))


# Generate the heatmap
pheatmap(scaled_markers_complete, 
         main = "olink changes complete",
         #color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)


#subset data to ket_del and med_del
subset_data_complete <- scaled_data_complete[scaled_data_complete$diet.phase %in% c("med_del", "keto_del"), ]
immune_data_for_heatmap_complete <- subset_data_complete[, -c(1:6)]

# Generate the heatmap
pheatmap(immune_data_for_heatmap_complete, 
         main = "olink changes complete",
         color = custom_colors, 
         cluster_rows = FALSE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)

# order by participant ID
subset_data_complete_ordered <- subset_data_complete[order(subset_data_complete$part_id), ]
immune_data_for_heatmap_complete_ordered <- subset_data_complete_ordered[, -c(1:6)]

# Generate the heatmap
pheatmap(immune_data_for_heatmap_complete_ordered, 
         main = "olink changes complete",
         color = custom_colors, 
         cluster_rows = FALSE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)

### Testing whether changes in immune markers are more similar within participant or within diet  #####

# Calculate Euclidean distance
# Select only the immune markers for distance calculation
distance_data <- subset_data_complete_ordered %>%
  select(part_id, diet.phase, 7:84) 

# Calculate the distance matrix
euclidean_distance <- dist(distance_data[, -c(1, 2)], method = "euclidean")  # Exclude participant_id and diet_phase

# Convert the distance object to a matrix
distance_matrix <- as.matrix(euclidean_distance)

# Perform hierarchical clustering
hc <- hclust(euclidean_distance)

# Plot the dendrogram
plot(hc, main = "Hierarchical Clustering of Immune Markers (Euclidean Distance)")

# Calculate within-participant distance
within_participant_distance <- sapply(unique(distance_data$part_id), function(participant) {
  # Get the indices of the current participant
  participant_indices <- which(distance_data$part_id == participant)
  # Extract the distances for the current participant
  participant_distances <- distance_matrix[participant_indices, participant_indices]
  # Calculate the mean distance for the current participant
  mean(participant_distances[upper.tri(participant_distances)])  # Use upper.tri to avoid double counting
})

print(within_participant_distance)

# Convert the distance matrix to a data frame
distance_df <- as.data.frame(distance_matrix)
# Optionally, you can add row names as a new column if you want to keep track of participant IDs
distance_df <- cbind(part_ID = rownames(distance_df), distance_df)
# Export the data frame to a CSV file
write.csv(distance_df, file = "distance_matrix.csv", row.names = FALSE) 

# Calculate between-diet distance
between_diet_distance <- sapply(unique(distance_data$diet.phase), function(diet) {
  # Get the indices of the current diet
  diet_indices <- which(distance_data$diet.phase == diet)
  # Extract the distances for the current diet
  diet_distances <- distance_matrix[diet_indices, diet_indices]
  # Calculate the mean distance for the current diet
  mean(diet_distances[upper.tri(diet_distances)])  # Use upper.tri to avoid double counting
})

# Print the results
print(between_diet_distance)

## Plot distances as box plot ##
olink_distance <- read.csv("olink_distances.csv")
data_long <- pivot_longer(olink_distance, cols = everything(), names_to = "Variable", values_to = "Value")

# Visualize olink distances in box plots
ggplot(data_long, aes(x = Variable, y = Value, fill = Variable)) +
  geom_boxplot(width = 0.6, alpha = 0.3) +  # Adjust position
  geom_jitter(width = 0.1, size = 1, color = "black", alpha = 0.3) +
  labs(title = "Box Plot of Distances", x = "Category", y = "Ave Bray-curtis distance") +
  theme_minimal()+
  theme(axis.text.x = element_blank())


### Generate PCoA plot of olink data. ####
# Perform PCoA
pcoa_result <- cmdscale(distance_matrix, k = 2, eig = TRUE)
# Create a data frame for plotting
pcoa_df <- data.frame(Sample = rownames(pcoa_result$points), 
                      PCoA1 = pcoa_result$points[, 1], 
                      PCoA2 = pcoa_result$points[, 2])
# Add treatment groups 
pcoa_df$diet.phase <- subset_data_complete_ordered$diet.phase[match(rownames(pcoa_df), rownames(subset_data_complete_ordered))]  

# Calculate centroids for each treatment group
Centroids_olink <- pcoa_df %>%
  group_by(diet.phase) %>%
  summarise(
    centroid_PCoA1 = mean(PCoA1),
    centroid_PCoA2 = mean(PCoA2)
  )
# View the centroids to ensure they're calculated
print(Centroids_olink)

# Generate the PCoA plot colored by Treatment
ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = diet.phase)) +
  geom_point() +
  labs(x = "PCoA 1", y = "PCoA 2") +
  theme_classic() + 
  stat_ellipse() +
  #geom_text(aes(label = Sample), vjust = -1, hjust = 0.5, size = 3) +  # Adjust vjust and hjust for label positioning
  # Add centroids with matching colors
  geom_point(data = Centroids_olink, 
             aes(x = centroid_PCoA1, y = centroid_PCoA2, color = diet.phase),  # Use color for centroids
             size = 5, 
             shape = 21,  # Use shape 21 for filled circles with outlines
             fill = "black",  # Fill color for the centroid circles
             stroke = 1.5)  # Adjust stroke width for the outline

# Testing whether olink data is different between diet groups
adonis_result <- adonis2(distance_matrix ~ diet.phase, data = distance_data)
print(adonis_result_BL)



#### Baseline analysis by diagnosis  ####
# Read olink baseline data
olink_BLdata <- read.csv("KM_olink_BL.csv")

# Separate metadata and immune markers
metadata_BL <- olink_BLdata[, 1:5]  # Adjust the column indices based on meta data columns
immune_markers_BL <- olink_BLdata[, 6:78]  # Adjust to include all immune marker columns
# Set PPID as row names
rownames(immune_markers_BL) <- olink_BLdata$part_id
rownames(metadata_BL) <- olink_BLdata$part_id

# Scaling the numeric data
numeric_data_scaled_BL <- scale(immune_markers_BL)

# Re-calculate the distance matrix on the scaled data
distance_matrix_BL <- dist(numeric_data_scaled_BL, method = "euclidean")
# Then perform PCoA 
pcoa_results_BL <- cmdscale(distance_matrix_BL, eig = TRUE, k = 2)

# Create a data frame for plotting
pcoa_df_BL <- data.frame(part_id = rownames(pcoa_results_BL$points), 
                      PCoA1 = pcoa_results_BL$points[, 1], 
                      PCoA2 = pcoa_results_BL$points[, 2])
# Add treatment groups 
pcoa_df_BL$diagnosis <- metadata_BL$diagnosis[match(rownames(pcoa_df_BL), rownames(metadata_BL))]  

# Calculate centroids for each diagnosis group
Centroids_olinkBL <- pcoa_df_BL %>%
  group_by(diagnosis) %>%
  summarise(
    centroid_PCoA1 = mean(PCoA1),
    centroid_PCoA2 = mean(PCoA2)
  )
# View the centroids to ensure they're calculated
print(Centroids_olinkBL)

# Generate the PCoA plot colored by Treatment
ggplot(pcoa_df_BL, aes(x = PCoA1, y = PCoA2, color = diagnosis)) +
  geom_point() +
  labs(x = "PCoA 1", y = "PCoA 2") +
  theme_classic() + 
  stat_ellipse() +
  #geom_text(aes(label = Sample), vjust = -1, hjust = 0.5, size = 3) +  # Adjust vjust and hjust for label positioning
  # Add centroids with matching colors
  geom_point(data = Centroids_olinkBL, 
             aes(x = centroid_PCoA1, y = centroid_PCoA2, color = diagnosis),  # Use color for centroids
             size = 5, 
             shape = 21,  # Use shape 21 for filled circles with outlines
             fill = "black",  # Fill color for the centroid circles
             stroke = 1.5)  # Adjust stroke width for the outline

# Testing whether olink data is different between diagnosis groups
adonis_result_BL <- adonis2(distance_matrix_BL ~ diagnosis, data = metadata_BL)
print(adonis_result_BL)

#### BL only looking at pro-inflam markers ####
# Read olink baseline data
olink_BLdata_pro <- read.csv("KM_olink_BLpro.csv")

# Separate metadata and immune markers
metadata_BLpro <- olink_BLdata_pro[, 1:5]  # Adjust the column indices based on meta data columns
immune_markers_BLpro <- olink_BLdata_pro[, 6:28]  # Adjust to include all immune marker columns
# Set PPID as row names
rownames(immune_markers_BLpro) <- olink_BLdata_pro$part_id
rownames(metadata_BLpro) <- olink_BLdata_pro$part_id

# Scaling the numeric data
numeric_data_scaled_BLpro <- scale(immune_markers_BLpro)

# add diagnosis information
diagnosis_info <- metadata_BLpro$diagnosis  
# Combine part_id and diagnosis for row names
new_row_names <- paste(rownames(numeric_data_scaled_BLpro), diagnosis_info, sep = " - ")
# Update row names of the scaled data
rownames(numeric_data_scaled_BLpro) <- new_row_names

#generate heat map of proinflam markers at BL 
pheatmap(numeric_data_scaled_BLpro, 
         main = "Proinflam at BL",
         color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)

### Looking at signficant markers only from WFKD delta diet data
# Read olink baseline data
olink_sig <- read.csv("olink_K_sig_avedelta_fixed.csv")

# Separate metadata and immune markers
metadata <- olink_sig[, 1:1]  # Adjust the column indices based on meta data columns
immune_markers <- olink_sig[, 2:13]  # Adjust to include all immune marker columns

# Set Sample name as row names
rownames(immune_markers) <- olink_sig$diet.phase
#rownames(metadata) <- olink_sig$diet.phase

# Scaling the numeric data
data_scaled <- scale(immune_markers)

#generate heat map of proinflam marker deltas by diet 
pheatmap(immune_markers, 
         main = "Marker deltas by WFKD",
         color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)



### Looking at proinflamm markers only from delta diet data
# Read olink baseline data
olink_delta_pro <- read.csv("KM_olink_deltas_pro.csv")

# Separate metadata and immune markers
metadata_delta_pro <- olink_delta_pro[, 1:5]  # Adjust the column indices based on meta data columns
immune_markers_deltapro <- olink_delta_pro[, 6:28]  # Adjust to include all immune marker columns

# Set Sample name as row names
rownames(immune_markers_deltapro) <- olink_delta_pro$sample
rownames(metadata_delta_pro) <- olink_delta_pro$sample

# Scaling the numeric data
numeric_data_scaled_deltapro <- scale(immune_markers_deltapro)

#generate heat map of proinflam marker deltas by diet 
pheatmap(numeric_data_scaled_deltapro, 
         main = "Proinflam deltas by all diet phases",
         color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)


# Combine metadata with scaled immune markers
scaled_data_pro <- cbind(metadata_delta_pro, numeric_data_scaled_deltapro)

#subset data to ket_del and med_del
subset_data_pro <- scaled_data_pro[scaled_data_pro$diet.phase %in% c("med_del", "keto_del"), ]
immune_data_for_heatmap_pro <- subset_data_pro[, -c(1:5)]

# Generate the heatmap
pheatmap(immune_data_for_heatmap_pro, 
         main = "olink changes pro only delivery phase",
         color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)

#subset data to keto_del or Med_del
subset_data_pro_M <- scaled_data_pro[scaled_data_pro$diet.phase %in% c("med_del"), ]
immune_data_for_heatmap_pro_M <- subset_data_pro_M[, -c(1:5)]

# Generate the heatmap
pheatmap(immune_data_for_heatmap_pro_M, 
         main = "olink changes pro only Med delivery phase",
         color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)

# box plots of signif pro-inflam markers #
marker <- read.csv("cst5.csv")

# Visualize olink distances in box plots
ggplot(marker, aes(x = diet.phase, y = il18_inflam, fill = diet.phase)) +
  geom_boxplot(width = 0.6, alpha = 0.3) +  # Adjust position
  geom_jitter(width = 0.1, size = 1, color = "black", alpha = 0.3) +
  labs(title = "il18", x = "diet", y = "delta from BL") +
  theme_minimal()



### Looking at antiinflamm markers only from delta diet data
# Read olink baseline data
olink_delta_anti <- read.csv("KM_olink_deltas_anti.csv")

# Separate metadata and immune markers
metadata_delta_anti <- olink_delta_anti[, 1:5]  # Adjust the column indices based on meta data columns
immune_markers_delta_anti <- olink_delta_anti[, 12:13]  # Adjust to include all immune marker columns

# Set Sample name as row names
rownames(immune_markers_delta_anti) <- olink_delta_anti$sample
rownames(metadata_delta_anti) <- olink_delta_anti$sample

# Scaling the numeric data
numeric_data_scaled_delta_anti <- scale(immune_markers_delta_anti)

#generate heat map of proinflam marker deltas by diet 
pheatmap(numeric_data_scaled_delta_anti, 
         main = "Antiinflam deltas by diet phases",
         color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)


### Looking at selected markers only from delta diet data
# Read olink baseline data
olink_delta_select <- read.csv("KM_olink_delta_selected2.csv")

# Separate metadata and immune markers
metadata_delta_select <- olink_delta_select[, 1:5]  # Adjust the column indices based on meta data columns
immune_markers_delta_select <- olink_delta_select[, 6:9]  # Adjust to include all immune marker columns

# Set Sample name as row names
rownames(immune_markers_delta_select) <- olink_delta_select$sample
rownames(metadata_delta_select) <- olink_delta_select$sample

# Scaling the numeric data
numeric_data_scaled_delta_select <- scale(immune_markers_delta_select)

#generate heat map of proinflam marker deltas by diet 
pheatmap(numeric_data_scaled_delta_select, 
         main = "Select deltas by diet phases",
         color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_row = 6)


# Read combined data with SCFA,Olink,Shannon
combined_data <- read.csv("combined_deltas.csv")

# Separate metadata and immune markers
metadata <- combined_data[, 1:2]  # Adjust the column indices based on meta data columns
markers <- combined_data[, 3:34]  # Adjust to include all marker columns

# Scale the markers
scaled_markers <- scale(markers)
# Set PPID as row names
rownames(scaled_markers) <- combined_data$PPID

#set colors
# Define the number of colors for negative and positive values
num_neg_colors <- 50  # Number of colors for negative values
num_pos_colors <- 50  # Number of colors for positive values
# Create the custom color palette
custom_colors <- c(colorRampPalette(c("blue", "white"))(num_neg_colors), 
                   "white", 
                   colorRampPalette(c("white", "red"))(num_pos_colors))


# Generate the heatmap
pheatmap(scaled_markers, 
         main = "marker changes",
         color = custom_colors, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE)

## Generate box plots of sig markers ##
markers <- read.csv("cst5.csv")
data_long <- pivot_longer(markers, cols = everything(), names_to = "Variable", values_to = "Value")

# Visualize olink distances in box plots
ggplot(data_long, aes(x = Variable, y = Value)) +
  geom_boxplot(width = 0.6, alpha = 0.3) +  # Adjust position
  geom_jitter(width = 0.1, size = 1, color = "black", alpha = 0.3) +
  labs(title = "CST5", x = "Category", y = "Ave Bray-curtis distance") +
  theme_minimal()

# Read the CSV file
data <- read.csv("KD_sig.csv")  

# Reshape the data to long format
data_long <- pivot_longer(data, 
                          cols = everything(), 
                          names_to = c("Marker", "Timepoint"), 
                          names_sep = "_")  

library(gridExtra)

# Step 3: Create individual box plots for each marker

plot_list <- list()


# Create a box plot for each marker
markers <- unique(data_long$Marker)

for (marker in markers) {
  p <- ggplot(data_long[data_long$Marker == marker, ], aes(x = Timepoint, y = value)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_jitter(width = 0.2, size = 1, color = "black", alpha = 0.3) +  # Removed position argument
    labs(title = paste(marker), x = NULL, y = NULL) +
    theme_minimal()+
    theme(axis.title.x = element_blank(),  # Remove x-axis title
          
          axis.title.y = element_blank(),  # Remove y-axis title
          
          panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add box around each plot
          
          plot.title = element_text(hjust = 0.5, size = 14),  # Center the title
          
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10))  # Add margin around the plot
  
  
  plot_list[[marker]] <- p  # Store each plot in a list
}

# Combine the plots into one figure
grid.arrange(grobs = plot_list, ncol = 3)  # Arrange plots in a single row



###heatmap###
  
library(tidyverse)
library(RColorBrewer)

md   <- read_csv("keto_med_metadata.csv")
data <- read_csv("KetoMed_Olink_RAWData.csv")
S2   <- read.csv("TableS2.csv", header = FALSE)

#Clean + label
data2 <- data %>%
  mutate(redcap_event_name = str_replace(redcap_event_name, "b_arm", "_arm"),
         visit_number = as.numeric(str_extract(redcap_event_name, "(?<=visit_)\\d+"))) %>%
  rename(participant_id = record_id)

protein_markers <- names(data2)[3:94]   # il8 ... csf1

md <- md %>%
  filter(stool_order >= 2) %>%
  distinct(participant_id, stool_order, .keep_all = TRUE) %>%
  select(participant_id, stool_order, FirstDiet, Diet, Group)

data_merged <- data2 %>%
  left_join(md, by = c("participant_id", "visit_number" = "stool_order")) %>%
  mutate(DietPhase = case_when(
    visit_number == 2 ~ "Baseline",
    FirstDiet == "Med first"  & visit_number == 3 ~ "Med_Delivery",
    FirstDiet == "Med first"  & visit_number == 4 ~ "Med_Own",
    FirstDiet == "Med first"  & visit_number == 5 ~ "Keto_Delivery",
    FirstDiet == "Med first"  & visit_number == 6 ~ "Keto_Own",
    FirstDiet == "Keto first" & visit_number == 3 ~ "Keto_Delivery",
    FirstDiet == "Keto first" & visit_number == 4 ~ "Keto_Own",
    FirstDiet == "Keto first" & visit_number == 5 ~ "Med_Delivery",
    FirstDiet == "Keto first" & visit_number == 6 ~ "Med_Own",
    TRUE ~ NA_character_
  ))

#Change from baseline 
baseline_values <- data_merged %>%
  filter(DietPhase == "Baseline") %>%
  select(participant_id, all_of(protein_markers)) %>%
  rename_with(~ paste0(., "_baseline"), all_of(protein_markers))

data_deltas <- data_merged %>%
  filter(DietPhase != "Baseline", !is.na(DietPhase)) %>%
  left_join(baseline_values, by = "participant_id")

for (marker in protein_markers) {
  data_deltas[[paste0(marker, "_delta")]] <-
    data_deltas[[marker]] - data_deltas[[paste0(marker, "_baseline")]]
}

# Keep only Table S2 markers
markers_keep <- S2 %>%
  mutate(marker = str_to_lower(str_trim(V1))) %>%
  filter(marker != "") %>%
  distinct(marker) %>%
  pull(marker)

delta_cols <- paste0(intersect(protein_markers, markers_keep), "_delta")

# Per-participant mean delta (average Delivery + Own within each diet)
participant_means_long <- data_deltas %>%
  mutate(Diet = if_else(str_detect(DietPhase, "Keto"), "Keto", "Med")) %>%
  select(participant_id, Diet, all_of(delta_cols)) %>%
  pivot_longer(all_of(delta_cols), names_to = "marker", values_to = "delta") %>%
  mutate(marker = str_remove(marker, "_delta$")) %>%
  group_by(participant_id, Diet, marker) %>%
  summarise(mean_delta = mean(delta, na.rm = TRUE), .groups = "drop")

# Within-diet one-sample t-test vs 0
within_diet_tests <- participant_means_long %>%
  filter(!is.na(mean_delta)) %>%
  group_by(Diet, marker) %>%
  summarise(
    n          = n(),
    p_value    = if (n() >= 3) t.test(mean_delta, mu = 0)$p.value else NA_real_,
    mean_delta = mean(mean_delta),
    .groups = "drop"
  ) %>%
  group_by(Diet) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

#Order markers
order_tbl <- within_diet_tests %>%
  select(Diet, marker, mean_delta, p_adj) %>%
  pivot_wider(names_from = Diet, values_from = c(mean_delta, p_adj)) %>%
  mutate(
    sig_keto = !is.na(p_adj_Keto) & p_adj_Keto < 0.05,
    sig_med  = !is.na(p_adj_Med)  & p_adj_Med  < 0.05
  ) %>%
  filter(sig_keto | sig_med) %>%
  mutate(
    sig_group = case_when(sig_keto & sig_med ~ 1,   # both
                          sig_keto           ~ 2,   # Keto only
                          TRUE               ~ 3),  # Med only
    dir_red_first = if_else((mean_delta_Keto + mean_delta_Med) / 2 > 0, 1, 2),
    strength = case_when(
      sig_group == 1 ~ (abs(mean_delta_Keto) + abs(mean_delta_Med)) / 2,
      sig_group == 2 ~ abs(mean_delta_Keto),
      TRUE           ~ abs(mean_delta_Med)
    )
  ) %>%
  arrange(sig_group, dir_red_first, desc(strength))

marker_order <- order_tbl$marker

#heatmap
heatmap_data_sig <- within_diet_tests %>%
  filter(marker %in% marker_order) %>%
  mutate(
    sig_lab = case_when(p_adj < 0.001 ~ "***",
                        p_adj < 0.01  ~ "**",
                        p_adj < 0.05  ~ "*",
                        TRUE          ~ ""),
    Diet   = factor(Diet, levels = c("Keto", "Med"), labels = c("WFKD", "Med")),
    marker = factor(marker, levels = marker_order)
  )


lim <- max(abs(heatmap_data_sig$mean_delta), na.rm = TRUE)


p1 <- ggplot(heatmap_data_sig, aes(x = marker, y = Diet, fill = mean_delta)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_lab), size = 5) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-lim, lim),
                       name = "Mean Change\nfrom baseline") +
  theme_classic() +
  theme(axis.line  = element_line(color = "black", linewidth = 0.8),
        axis.ticks = element_line(color = "black", linewidth = 0.8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.caption = element_text(hjust = 0, size = 10)) +
  labs(title = "Mean Change from Baseline (significant markers only)",
       x = "Protein Marker",
       y = "Diet",
       caption = "* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001 (within-diet vs baseline)")

print(p1)

