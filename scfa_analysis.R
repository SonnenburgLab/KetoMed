library(vegan)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(nlme)
library(pheatmap)
library(reshape2)

# set working directory 
setwd("/Users/ericas1/Desktop/keto_med_diet/SCFA")
# Read SCFA data
SCFA_data <- read.csv("/Users/samgoldman/Documents/keto_med/from_erica/scfa.csv")
# Exclude samples with "baseline" in the "diet" column
SCFA_data_filtered <- subset(SCFA_data, diet != "baseline")
# Set the order of the diet groups
SCFA_data_filtered$diet <- factor(SCFA_data_filtered$diet, levels = c("pre_baseline", "keto_del", "keto_own"))

# Exclude all samples except prebaseline and keto
SCFA_data_keto <- subset(SCFA_data, diet %in% c("pre_baseline", "keto_del", "keto_own"))
# Set the order of the diet groups
SCFA_data_keto$diet <- factor(SCFA_data_keto$diet, levels = c("pre_baseline", "keto_del", "keto_own"))


# Exclude all samples except prebaseline
SCFA_data_bl <- subset(SCFA_data, diet %in% c("pre_baseline"))
# Set the order of the diet groups
SCFA_data_blmed$diet <- factor(SCFA_data_blmed$diet, levels = c("pre_baseline", "med_del"))

# Perform the correlation test between SCFA types
cor_test_result <- cor.test(SCFA_data_filtered$propionic, SCFA_data_filtered$acetate, method = "pearson")
# Print the results of the correlation test
print(cor_test_result)

# Exclude all samples except prebaseline and keto_del
SCFA_data_blketo <- subset(SCFA_data, diet %in% c("pre_baseline", "keto_del"))
# Set the order of the diet groups
SCFA_data_blmed$diet <- factor(SCFA_data_blmed$diet, levels = c("pre_baseline", "med_del"))

custom_colors <- c("pre_baseline" = "#b7b7b7ff","keto_del" = "#c7885fff", "keto_own" = "#f5c09fff", "med_del" = "#669289ff", "med_own" = "#afd3caff")

# Visualize scfa in box plots
ggplot(SCFA_data_blketo, aes(x = diet, y = propionic, fill = diet)) +  # Replace SampleType with actual treatment column
  #geom_violin(trim = TRUE, alpha = 0.4) +  # Fill color based on treatment group
  scale_fill_manual(values = custom_colors) + 
  geom_boxplot(width = 0.4, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, color = "black", alpha = 0.5) +
  labs(title = "Propionic by diet", 
       x = "Diet", 
       y = "propionate") +
  facet_wrap(~ FirstDiet) +  # Create separate panels for each diagnosis
  theme_minimal()


# Visualize scfa in box plots of baseline only by diagnosis
diagnosis_colors <- c("T2_diabetes" = "#741b47ff","Prediabetes" = "#d5a6bdff")

ggplot(SCFA_data_bl, aes(x = diagnosis, y = propionic, fill = diagnosis)) +  # Replace SampleType with actual treatment column
  #geom_violin(trim = TRUE, alpha = 0.4) +  # Fill color based on treatment group
  scale_fill_manual(values = diagnosis_colors) + 
  geom_boxplot(width = 0.4, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, color = "black", alpha = 0.5) +
  labs(title = "Propionic by diet", 
       x = "Diagnosis", 
       y = "propionic") +
  theme_minimal()


# Read SCFA data that contains 16s diversity information
SCFA_div_data <- read.csv("scfa_div.csv")
# Set the order of the diet groups
SCFA_div_data$diet <- factor(SCFA_data_filtered$diet, levels = c("pre_baseline", "keto_del", "keto_own", "med_del", "med_own"))

# Fit a mixed effects model with an interaction term for all data
model_scfa_div <- lmer(butyrate ~ Simpson + (1|PPID) + (1|diet), data = SCFA_div_data)
# Get model summary
summary_scfa_div <- summary(model_scfa_div)
print(summary_scfa_div)
# Run ANOVA to test significance of the interaction
interaction_anova <- anova(model_scfa_div)
print(interaction_anova)

# Create a scatter plot with a regression line
ggplot(SCFA_div_data, aes(x = Shannon, y = butyrate)) +
  geom_point(aes(color = diet), alpha = 0.7) +  # Points colored by arm
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Regression line
  labs(title = "Impact of Observed on Butyrate",
       x = "Observed",
       y = "butyrate") +
  theme_minimal()

## only look at interaction between diversity and scfa for med diet
# Exclude all samples except med_del
SCFA_div_data <- read.csv("scfa_div.csv")
SCFA_div_data_filtered <- subset(SCFA_div_data, diet %in% c("med_del", "med_own"))

# Fit a mixed effects model with an interaction term
model_scfa_div_filtered <- lmer(butyrate ~ Observed + (1|PPID) + (1|diet), data = SCFA_div_data_filtered)
# Get model summary
summary_scfa_div_filtered <- summary(model_scfa_div_filtered)
print(summary_scfa_div_filtered)
# Run ANOVA to test significance of the interaction
interaction_anova_filtered <- anova(model_scfa_div_filtered)
print(interaction_anova_filtered)

# Create a scatter plot with a regression line
ggplot(SCFA_div_data_filtered, aes(x = Simpson, y = butyrate)) +
  geom_point(aes(color = diet), alpha = 0.7) +  # Points colored by arm
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Regression line
  labs(title = "Impact of Observed on Butyrate",
       x = "Observed",
       y = "butyrate") +
  theme_minimal()


## look at interaction between diversity and change in SCFA in med_del diet
SCFA_del_data <- read.csv("scfa_div_delta.csv")

# Calculate correlation between diversity and SCFA
correlation_result <- cor(SCFA_del_data$d_acetate, SCFA_del_data$d_butyrate, method = "pearson")
# Print the correlation result
print(correlation_result)

# Perform the correlation test
cor_test_result <- cor.test(SCFA_del_data$d_succinate, SCFA_del_data$d_butyrate, method = "pearson")
# Print the results of the correlation test
print(cor_test_result)


# Create a scatter plot 
ggplot(SCFA_del_data, aes(x = d_acetate, y = d_butyrate)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a linear regression line
  labs(title = "Relationship between Change in Propionic and Change in Butyrate",
       x = "Change in propionic",
       y = "Change in Butyrate") +
  theme_minimal()

# Fit a linear model
linear_model <- lm(d_succinate ~ Shannon, data = SCFA_del_data)
# Print the summary of the linear model
summary(linear_model)


# Generate heat map of scfa data
# Subset the data to include only specific columns
selected_columns <- SCFA_del_data[, c("d_butyrate", "d_acetate", "d_propionic")]
# Set PPID as row names
rownames(selected_columns) <- SCFA_del_data$PPID
# Scale the data
scaled_data <- scale(selected_columns)

#set colors
# Define the number of colors for negative and positive values
num_neg_colors <- 50  # Number of colors for negative values
num_pos_colors <- 50  # Number of colors for positive values
# Create the custom color palette
custom_colors2 <- c(colorRampPalette(c("red", "white"))(num_neg_colors), 
                    "white", 
                   colorRampPalette(c("white", "blue"))(num_pos_colors))
# Define breaks for the color scale
breaks <- seq(min(scaled_data), max(scaled_data), length.out = length(custom_colors2) + 1)

# Generate the heatmap
pheatmap(scaled_data, 
         main = "SCFA changes from baseline to Med_del",
         color = custom_colors2, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE)

