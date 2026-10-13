library(ape) # Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq
library(ggplot2) # Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.
library(phyloseq) # The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.
library(plotly) # A package to create interactive web graphics of use in 3D plots
library(vegan) # The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.
library(philr) # This package provides functions for the analysis of compositional data 
library(tidyverse) # This package is designed to make it easy to install and load multiple 'tidyverse' packages in a single step
library(adespatial) # Tools for the multiscale spatial analysis of multivariate data
library(devtools) # Make package development easier by providing R functions that simplify and expedite common tasks
library("pander") # provide a minimal and easy tool for rendering R objects into Pandoc's markdown
library(grid) # support data visualization
library(gridExtra)  # support data visualization
library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.
library(png) # Figure download
library("ggdendro") #set of tools for dendrograms and tree plots using 'ggplot2'
library(ggpubr) # publication quality figures, based on ggplot2
library(RColorBrewer) # nice color options
library(ggplotify) # Data visualization
library(breakaway) # alpha diversity estimator
library(magrittr) # read in phyloseq object and filter 
library(microbiome)
library(microbiomeutilities)
library(pairwiseAdonis)
library(ggalluvial)

# read in phyloseq from metagenomic data Sailen sent me
physeq <- readRDS("/Users/samgoldman/Documents/keto_med/from_erica/ps_mpa3_species_ketomed_050125.rds") 
# Filtering out ASVs with less than 10 reads across all phyloseq object
physeq <- prune_taxa(taxa_sums(physeq) > 10, physeq)

# Load sample metadata from CSV metadata file
sample_metadata <- read.csv("/Users/samgoldman/Documents/keto_med/from_erica/keto_med_metadata.csv", header = TRUE, row.names = 1) 
# Transfer metadata to the phyloseq object
# Create a phyloseq sample data object
sample_data_physeq1 <- sample_data(sample_metadata)
# Assign this sample data to existing phyloseq object
sample_data(physeq) <- sample_data_physeq1
# Verify the sample data has been added
print(sample_data(physeq))


# generate distribution plot of filtered reads to determine how to rarefy
# Extract total reads per sample using sample_sums
total_reads <- sample_sums(physeq)

# smallest number of reads per sample is 47050484
# rarefy data to 47050484, Set a seed for reproducibility
set.seed(123)  
rarefied_physeq <- rarefy_even_depth(physeq, sample.size = 47050484, replace = FALSE)

# Transfer metadata to the rarefied phyloseq object
sample_data(rarefied_physeq) <- sample_data_physeq1
# Verify the sample data has been added
print(sample_data(rarefied_physeq))

#### ALPHA DIVERSITY ANALYSIS ######
# Calculate alpha diversity metrics using unpruned data
alpha_diversity <- estimate_richness(rarefied_physeq)  
# Convert alpha diversity metrics to data frame
alpha_diversity_df <- data.frame(alpha_diversity)

# Remove the "X" prefix from the row names, not sure why that happened
rownames(alpha_diversity_df) <- gsub("^X", "", rownames(alpha_diversity_df))

# Combine alpha diversity data with sample metadata
alpha_diversity_combined <- cbind(data.frame(sample_data(rarefied_physeq)), alpha_diversity_df)

#create baseline physeq of unrarefied data
physeq_baseline_unrare <- subset_samples(physeq, Diet == "Baseline")

#Filter the phyloseq object for baseline only
physeq_baseline <- subset_samples(rarefied_physeq, Diet == "Baseline")
# Calculate alpha diversity metrics using unpruned data on baseline only
alpha_diversity_BL <- estimate_richness(physeq_baseline)  
# Convert alpha diversity metrics to data frame
alpha_diversity_df_BL <- data.frame(alpha_diversity_BL)
# Combine alpha diversity data with sample metadata
alpha_diversity_combined_BL <- cbind(data.frame(sample_data(physeq_baseline)), alpha_diversity_df_BL)

#####    PLots used for paper figures    ####
diagnosis_colors <- c("T2_diabetes" = "#741b47ff","Prediabetes" = "#d5a6bdff")

# Visualize alpha diversity at baseline by diagnosis in violin/box plots
ggplot(alpha_diversity_combined_BL, aes(x = Eligibility_diagnosis, y = Shannon, fill = Eligibility_diagnosis)) +  # Replace SampleType with actual treatment column
  scale_fill_manual(values = diagnosis_colors) +
  geom_violin(trim = TRUE, alpha = 0.5) +  # Fill color based on treatment group
  geom_jitter(width = 0.1, size = 2, color = "black", alpha = 0.5) +
  geom_boxplot(width = 0.4, alpha = 0.5) +
  labs(title = "Alpha Diversity by Treatment", 
       x = "Diet", 
       y = "Shannon") +
  theme_minimal()

diet_colors <- c("Baseline" = "#b7b7b7ff","Keto" = "#D59266", "Med" = "#6D9B91")

alpha_diversity_combined$Diet <- factor(alpha_diversity_combined$Diet, 
                         levels = c("Baseline", "Keto", "Med"))


# Visualize alpha diversity by treatment in violin/box plots
ggplot(alpha_diversity_combined, aes(x = Diet, y = Shannon, fill = Diet)) +  # Replace SampleType with actual treatment column
  scale_fill_manual(values = diet_colors) +
  geom_violin(trim = TRUE, alpha = 0.7) +  # Fill color based on treatment group
  geom_jitter(width = 0.1, size = 2, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.4, alpha = 0.5) +
  labs(title = "Alpha Diversity by Treatment", 
       x = "Diet", 
       y = "Shannon") +
  facet_wrap(~ FirstDiet) + 
  theme_minimal()

# Visualize alpha diversity by treatment in violin plots
ggplot(alpha_diversity_combined, aes(x = Diet, y = Simpson, fill = Diet)) +  # Replace SampleType with actual treatment column
  geom_violin(trim = TRUE, alpha = 0.6) +  # Fill color based on treatment group
  geom_jitter(width = 0.1, size = 1, color = "black", alpha = 0.5) +
  labs(title = "Alpha Diversity by Treatment", 
       x = "Diet", 
       y = "Simpson") +
  theme_minimal()

# Visualize alpha diversity by diagnosis in box plots
ggplot(alpha_diversity_combined, aes(x = Diet, y = Simpson, fill = Diet)) +  # Replace SampleType with actual treatment column
  geom_boxplot() +
  labs(title = "Alpha Diversity by Treatment", 
       x = "Diet", 
       y = "Simpson") +
  facet_wrap(~ Eligibility_diagnosis) +  # Create separate panels for each diagnosis
  theme_minimal()

# Visualize alpha diversity by diagnosis in violin plots
ggplot(alpha_diversity_combined, aes(x = Diet, y = Simpson, fill = Eligibility_diagnosis)) +  # Replace SampleType with actual treatment column
  geom_violin(trim = TRUE, alpha = 0.6) +  # Fill color based on treatment group
  geom_boxplot(width = 0.4) +
  geom_jitter(width = 0.1, size = 1, color = "black", alpha = 0.5) +
  labs(title = "Alpha Diversity by Treatment", 
       x = "Diet", 
       y = "Simpson") +
  facet_wrap(~ Eligibility_diagnosis) +  # Create separate panels for each diagnosis
  theme_minimal()

# Visualize alpha diversity by diet sequence in box plots
ggplot(alpha_diversity_combined, aes(x = Diet, y = Simpson, fill = Diet)) +  # Replace SampleType with actual treatment column
  geom_boxplot() +
  labs(title = "Alpha Diversity by Treatment", 
       x = "Diet", 
       y = "Simpson") +
  facet_wrap(~ FirstDiet) +  # Create separate panels for each diagnosis
  theme_minimal()

# Visualize alpha diversity by diet sequence in violin plots
ggplot(alpha_diversity_combined, aes(x = Diet, y = Simpson, fill = Diet)) +  # Replace SampleType with actual treatment column
  geom_violin(trim = TRUE, alpha = 0.6) +  # Fill color based on treatment group
  geom_jitter(width = 0.1, size = 1, color = "black", alpha = 0.5) +
  labs(title = "Alpha Diversity by Treatment", 
       x = "Diet", 
       y = "Simpson") +
  facet_wrap(~ FirstDiet) +  # Create separate panels for each diagnosis
  theme_minimal()

# download alpha_diversity_combined in csv format
file_path <- "/Users/samgoldman/Documents/keto_med/alpha_diversity_combined.csv"
write.csv(alpha_diversity_combined, file = file_path)



####.  TAXONOMIC BOX and BAR PLOT ANALYSIS BY DIFFERENT LEVELS.   ######
# Generate composition bar plots
# plot taxa box plots, phylum
physeq_cla <- aggregate_top_taxa2(physeq, level = "Class", top = 10)
physeq_df_c <- microbiomeutilities::phy_to_ldf(physeq_cla, 
                                               transform.counts = "compositional")

# Create stacked bar  plots
ggplot(physeq_df_c, aes(x = Diet, y = Abundance, fill = Diet)) +
  geom_col(aes(fill = Class), position = "fill") +
  labs(title = "Stacked Bar Plot of Class by Diet",
       x = "Diet",
       y = "Relative Abundance") +
  theme_classic() +
  facet_wrap(~ Eligibility_diagnosis) 

write_csv(physeq_df_p, "KM_phy_abundance.csv")

# plot taxa box plots
physeq_fam <- aggregate_top_taxa2(physeq, level = "Family", top = 6)
physeq_df <- microbiomeutilities::phy_to_ldf(physeq_fam, 
                                             transform.counts = "compositional")
write.csv(relative_abundance_df, "relative_abundance_by_family.csv")

ggplot(physeq_df, aes(x = Family, y = Abundance, fill = Diet)) +
  geom_boxplot() +
  labs(title = "Box Plot of Families by Diet",
       x = "Family",
       y = "Relative Abundance") +
  theme_minimal() +
  facet_wrap(~ Eligibility_diagnosis) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility
  
  # plot taxa box plots, genus
  physeq_gen <- aggregate_top_taxa2(physeq, level = "Genus", top = 50)
  
  physeq_df_g <- microbiomeutilities::phy_to_ldf(physeq_gen, 
                                                 transform.counts = "compositional")
  
  filtered_data <- physeq_df_g %>%
    filter(Diet == "Baseline")  # Adjust the column name if necessary
  
  ggplot(filtered_data, aes(x = Genus, y = Abundance, fill = Eligibility_diagnosis)) +
    geom_boxplot() +
    labs(title = "Box Plot of Species by Diet",
         x = "Genus",
         y = "Relative Abundance") +
    theme_minimal() +
    #facet_wrap(~ Eligibility_diagnosis) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility
  
  write_csv(physeq_df_g, "KM_gen_abundance.csv")

# plot taxa box plots, species
  physeq_spe <- aggregate_top_taxa2(physeq, level = "Species", top = 30)
  
  physeq_df_s <- microbiomeutilities::phy_to_ldf(physeq_spe, 
                                               transform.counts = "compositional")
  
  filtered_data_s <- physeq_df_s %>%
    filter(Diet == "Baseline")  # Adjust the column name if necessary
  
  ggplot(filtered_data_s, aes(x = Species, y = Abundance, fill = Eligibility_diagnosis)) +
    geom_boxplot() +
    labs(title = "Box Plot of Species by Diet",
         x = "Species",
         y = "Relative Abundance") +
    theme_minimal() +
    #facet_wrap(~ Eligibility_diagnosis) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility

  write_csv(physeq_df_s, "KM_spe_abundance.csv")
  
### New taxonomic analysis at baseline with no filtering for top taxa.  ####
#Filter the phyloseq object for baseline only
physeq_baseline <- subset_samples(physeq, Diet == "Baseline")
#Convert phyloseq object into long-format data frame with relative abundances of microbial taxa
physeq_df_bl <- microbiomeutilities::phy_to_ldf(physeq_baseline, 
                                                 transform.counts = "compositional")
# Perform Mann-Whitney U test for each species
results_BL_spe <- physeq_df_bl %>%
  group_by(Genus) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Eligibility_diagnosis)$p.value,  # Mann-Whitney U test
    .groups = 'drop'
  )

#Adjust for multiple comparisons using Benjamini-Hochberg method
  results <- results_BL_spe %>%
    mutate(
      p_adjusted = p.adjust(p_value, method = "BH"),
      significant = p_adjusted < 0.05  # Determine significance
    )
  
#Identify significant species
  significant_species <- results %>%
    filter(significant)
  print(significant_species)
  
# NO significant species identified using above analysis between diagnosis at bl
# Only signficant genus was: Genus         p_value p_adjusted significant
# <chr>           <dbl>      <dbl> <lgl>      
# 1 Clostridium 0.0000314     0.0247 TRUE     
  

  
  
 
#####  TAXONOMIC ANALYSIS BETWEEN DIAGNOSIS AT BASELINE.  #####   
  # Check for signficant phyla between diagnosis arms at baseline
  physeq_phy <- aggregate_top_taxa2(physeq_baseline, level = "Phylum", top = 6)
  physeq_df_p <- microbiomeutilities::phy_to_ldf(physeq_phy, 
                                                 transform.counts = "compositional")                                               
# Filter the data for the specific diet and diagnosis groups
filtered_data <- physeq_df_p %>%
  filter(Diet == "Baseline", Eligibility_diagnosis %in% c("Prediabetes", "T2_diabetes"))
# Perform Mann-Whitney U test for each genus
results <- filtered_data %>%
  group_by(Phylum) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Eligibility_diagnosis)$p.value,  # Mann-Whitney U test
    .groups = 'drop'
  )
#Adjust for multiple comparisons using Benjamini-Hochberg method
results <- results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05  # Determine significance
  )
#Identify significant phyla
significant_phyla <- results %>%
  filter(significant)
print(significant_phyla)

#Create the stacked bar plot
ggplot(physeq_df_p, aes(x = OTUID, y = Abundance)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "Top 10 Most Abundant Phyla at Baseline in Pre",
       x = "Phyla",
       y = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Check for signficant families between diagnosis arms at baseline
  physeq_fam <- aggregate_top_taxa2(physeq_baseline_unrare, level = "Family", top = 6)
  physeq_df_f <- microbiomeutilities::phy_to_ldf(physeq_fam, 
                                               transform.counts = "compositional")                                            
# Filter the data for the specific diet and diagnosis groups
  filtered_data <- physeq_df_f %>%
    filter(Diet == "Baseline", Eligibility_diagnosis %in% c("Prediabetes", "T2_diabetes"))
# Perform Mann-Whitney U test for each genus
  results <- filtered_data %>%
    group_by(Family) %>%
    summarise(
      p_value = wilcox.test(Abundance ~ Eligibility_diagnosis)$p.value,  # Mann-Whitney U test
      .groups = 'drop'
    )
#Adjust for multiple comparisons using Benjamini-Hochberg method
  results <- results %>%
    mutate(
      p_adjusted = p.adjust(p_value, method = "BH"),
      significant = p_adjusted < 0.05  # Determine significance
    )
  #Identify significant family
  significant_families <- results %>%
    filter(significant)
  print(significant_families)


# Check for signficant genera between diagnosis arms at baseline
#Filter the data for the specific diet and diagnosis groups
  physeq_gen <- aggregate_top_taxa2(physeq_baseline_unrare, level = "Genus", top = 20)
  physeq_df_g <- microbiomeutilities::phy_to_ldf(physeq_gen, 
                                                 transform.counts = "compositional")   
  filtered_data <- physeq_df_g %>%
    filter(Diet == "Baseline", Eligibility_diagnosis %in% c("Prediabetes", "T2_diabetes"))
# Perform Mann-Whitney U test for each genus
results <- filtered_data %>%
  group_by(Genus) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Eligibility_diagnosis)$p.value,  # Mann-Whitney U test
    .groups = 'drop'
    )
#Adjust for multiple comparisons using Benjamini-Hochberg method
results <- results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05  # Determine significance
    )
#Identify significant genera
significant_genera <- results %>%
  filter(significant)
print(significant_genera)

# Generate box plot of signficant genus
#Filter to the genus of interest
df_med <- physeq_df_g %>% filter(Genus == "Mediterraneibacter")

ggplot(df_med, aes(x = Eligibility_diagnosis, y = Abundance, fill = Eligibility_diagnosis)) + 
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.7) + 
  geom_violin(trim = TRUE, alpha = 0.5) +  # Fill color based on treatment group
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8, color = "black") + 
  scale_fill_manual(values = diagnosis_colors) + 
  #scale_y_continuous(labels = percent_format(accuracy = 0.1)) + 
  labs( x = NULL, y = "Mediterraneibacter relative abundance", title = "Mediterraneibacter by diagnosis" ) + 
  theme_classic(base_size = 12)



#### PAPER FIGURES. #####
# Check for signficant in species between diagnosis arms at baseline
#Filter the data for the specific diet and diagnosis groups
#physeq_spe <- microbiome::aggregate_rare(physeq_baseline, level = "Species", detection = 75/100, prevalence = 75/100)
physeq_spe <- aggregate_top_taxa2(physeq_baseline_unrare, level = "Species", top = 100)
physeq_df_s <- microbiomeutilities::phy_to_ldf(physeq_spe, 
                                               transform.counts = "compositional")
filtered_data_s <- physeq_df_s %>%
  filter(Diet == "Baseline", Eligibility_diagnosis %in% c("Prediabetes", "T2_diabetes"))
# Perform Mann-Whitney U test for each species
results <- filtered_data_s %>%
  group_by(Species) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Eligibility_diagnosis)$p.value,  # Mann-Whitney U test
    .groups = 'drop'
  )
#Adjust for multiple comparisons using Benjamini-Hochberg method
results <- results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05  # Determine significance
  )
#Identify significant species
significant_species <- results %>%
  filter(significant)
print(significant_species)

write.csv(filtered_data_s, "filtered_data_s_BL.csv", row.names = FALSE)

#plot signficant species
df_CCQ <- physeq_df_s %>% filter(Species == "Candidatus_Cibionibacter_quicibialis")
df_RG <- physeq_df_s %>% filter(Species == "Ruminococcus_gnavus")

ggplot(df_CCQ, aes(x = Eligibility_diagnosis, y = Abundance, fill = Eligibility_diagnosis)) + 
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.7) + 
  geom_violin(trim = TRUE, alpha = 0.5) +  # Fill color based on treatment group
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8, color = "black") + 
  scale_fill_manual(values = diagnosis_colors) + 
  labs( x = NULL, y = "Candidatus_Cibionibacter_quicibialis rel abun") + 
  theme_classic(base_size = 12)


# Check for signficant in species between baseline and each diet
#Filter the data for baseline and keto diet using unrarefied physeq
physeq_spe <- aggregate_top_taxa2(physeq, level = "Species", top = 100)
physeq_df_s <- microbiomeutilities::phy_to_ldf(physeq_spe, 
                                               transform.counts = "compositional")
filtered_data_s <- physeq_df_s %>%
  filter(Diet %in% c("Baseline","Keto"))
# Perform Mann-Whitney U test for each species
results <- filtered_data_s %>%
  group_by(Species) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Diet)$p.value,  # Mann-Whitney U test
    .groups = 'drop'
  )
#Adjust for multiple comparisons using Benjamini-Hochberg method
results <- results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05  # Determine significance
  )
#Identify significant species
significant_species <- results %>%
  filter(significant)
print(significant_species)

#   Species              p_value p_adjusted significant
#<chr>                  <dbl>      <dbl> <lgl>      
#1 Bifidobacterium_longum   0.00142      0.0360  TRUE    
#2 Blautia_wexlerae    0.00110      0.0171 TRUE       
#3 Eubacterium_rectale 0.000379     0.0117 TRUE       
#4 Other               0.00216      0.0223 TRUE 
#5 Streptococcus_salivarius 0.0000532    0.00538 TRUE

#plot signficant species
df_BW <- filtered_data_s %>% filter(Species == "Blautia_wexlerae")
df_ER <- filtered_data_s %>% filter(Species == "Eubacterium_rectale")
df_BL <- filtered_data_s %>% filter(Species == "Bifidobacterium_longum")
df_SS <- filtered_data_s %>% filter(Species == "Streptococcus_salivarius")

ggplot(df_SS, aes(x = Diet, y = Abundance, fill = Diet)) + 
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.7) + 
  geom_violin(trim = TRUE, alpha = 0.5) +  # Fill color based on treatment group
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8, color = "black") + 
  scale_fill_manual(values = diet_colors) + 
  labs( x = NULL, y = "Streptococcus_salivarius rel abun") + 
  theme_classic(base_size = 12)

# Check for signficant in species between baseline and each diet
#Filter the data for baseline and keto diet using unrarefied physeq
physeq_spe <- aggregate_top_taxa2(physeq, level = "Species", top = 100)
physeq_df_s <- microbiomeutilities::phy_to_ldf(physeq_spe, 
                                               transform.counts = "compositional")
filtered_data_s <- physeq_df_s %>%
  filter(Diet %in% c("Baseline","Med"))
# Perform Mann-Whitney U test for each species
results <- filtered_data_s %>%
  group_by(Species) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Diet)$p.value,  # Mann-Whitney U test
    .groups = 'drop'
  )
#Adjust for multiple comparisons using Benjamini-Hochberg method
results <- results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05  # Determine significance
  )
#Identify significant species
significant_species <- results %>%
  filter(significant)
print(significant_species)

#   none significant


# Check for signficant in genus between baseline and each diet
#Filter the data for baseline and keto diet using unrarefied physeq
physeq_gen <- aggregate_top_taxa2(physeq, level = "Genus", top = 50)
physeq_df_g <- microbiomeutilities::phy_to_ldf(physeq_gen, 
                                               transform.counts = "compositional")
filtered_data_g <- physeq_df_g %>%
  filter(Diet %in% c("Baseline","Keto"))
# Perform Mann-Whitney U test for each species
results <- filtered_data_g %>%
  group_by(Genus) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Diet)$p.value,  # Mann-Whitney U test
    .groups = 'drop'
  )
#Adjust for multiple comparisons using Benjamini-Hochberg method
results <- results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05  # Determine significance
  )
#Identify significant genus
significant_genus <- results %>%
  filter(significant)
print(significant_genus)

#   Genus                         p_value p_adjusted significant
# <chr>                           <dbl>      <dbl> <lgl>      
# 1 Bifidobacterium              0.000205    0.00881 TRUE       
# 2 Lachnospiraceae_unclassified 0.000346    0.00881 TRUE       
# 3 Other                        0.00281     0.0477  TRUE   


# Check for signficant in genus between baseline and each diet
#Filter the data for baseline and keto diet using unrarefied physeq
physeq_gen <- aggregate_top_taxa2(physeq, level = "Genus", top = 20)
physeq_df_g <- microbiomeutilities::phy_to_ldf(physeq_gen, 
                                               transform.counts = "compositional")
filtered_data_g <- physeq_df_g %>%
  filter(Diet %in% c("Baseline","Med"))
# Perform Mann-Whitney U test for each species
results <- filtered_data_g %>%
  group_by(Genus) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Diet)$p.value,  # Mann-Whitney U test
    .groups = 'drop'
  )
#Adjust for multiple comparisons using Benjamini-Hochberg method
results <- results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05  # Determine significance
  )
#Identify significant species
significant_genus <- results %>%
  filter(significant)
print(significant_genus)

#   Differences between baseline and keto
#Genus                         p_value p_adjusted significant
#<chr>                           <dbl>      <dbl> <lgl>      
#1 Lachnospiraceae_unclassified 0.000346    0.00726 TRUE       
#2 Other                        0.00450     0.0473  TRUE       
#3 Streptococcus                0.00682     0.0478  TRUE   

# Differences between baseline and med
#   Genus              p_value p_adjusted significant
#<chr>                <dbl>      <dbl> <lgl>      
#1 Bacteroides        0.00569     0.0310 TRUE       
#2 Blautia            0.00302     0.0310 TRUE       
#3 Mediterraneibacter 0.00590     0.0310 TRUE       
#4 Streptococcus      0.00590     0.0310 TRUE 

#plot signficant species
df_Bac <- filtered_data_g %>% filter(Genus == "Bacteroides")
df_Bla <- filtered_data_g %>% filter(Genus == "Blautia")
df_Med <- filtered_data_g %>% filter(Genus == "Mediterraneibacter")
df_Ste <- filtered_data_g %>% filter(Genus == "Streptococcus")

ggplot(df_Ste, aes(x = Diet, y = Abundance, fill = Diet)) + 
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.7) + 
  geom_violin(trim = TRUE, alpha = 0.5) +  # Fill color based on treatment group
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8, color = "black") + 
  scale_fill_manual(values = diet_colors) + 
  labs( x = NULL, y = "Streptococcus rel abun") + 
  theme_classic(base_size = 12)


# Check for signficant in species between diagnosis arms during Keto diet
#Filter the data for the specific diet and diagnosis groups
physeq_spe <- aggregate_top_taxa2(physeq, level = "Species", top = 100)
physeq_df_s <- microbiomeutilities::phy_to_ldf(physeq_spe, 
                                               transform.counts = "compositional")
filtered_data_s <- physeq_df_s %>%
  filter(Diet == "Keto", Eligibility_diagnosis %in% c("Prediabetes", "T2_diabetes"))
# Perform Mann-Whitney U test for each species
results <- filtered_data_s %>%
  group_by(Species) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Eligibility_diagnosis)$p.value,  # Mann-Whitney U test
    .groups = 'drop'
  )
#Adjust for multiple comparisons using Benjamini-Hochberg method
results <- results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05  # Determine significance
  )
#Identify significant species
significant_species <- results %>%
  filter(significant)
print(significant_species)



### Top species for each diagnosis   ###.
physeq_t2d <- subset_samples(physeq_baseline, Eligibility_diagnosis == "T2_diabetes")
physeq_spe <- aggregate_top_taxa2(physeq_t2d, level = "Species", top = 11)
physeq_df_s <- microbiomeutilities::phy_to_ldf(physeq_spe, 
                                               transform.counts = "compositional")
write_csv(physeq_df_s, "t2d_top10spe_abundance.csv")

Create bar plot
ggplot(physeq_df_s, aes(x = OTUID, y = Abundance)) +
    geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "Top 10 Most Abundant Species at Baseline in T2D",
            x = "Species",
       y = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Same now for Pre_diabetes
physeq_pre <- subset_samples(physeq_baseline, Eligibility_diagnosis == "Prediabetes")
physeq_spe <- aggregate_top_taxa2(physeq_pre, level = "Species", top = 11)
physeq_df_s <- microbiomeutilities::phy_to_ldf(physeq_spe, 
                                               transform.counts = "compositional")
write_csv(physeq_df_s, "pre_top10spe_abundance.csv")

#Create the stacked bar plot
ggplot(physeq_df_s, aes(x = OTUID, y = Abundance)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "Top 10 Most Abundant Species at Baseline in Pre",
       x = "Species",
       y = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#####  TAXONOMIC ANALYSIS BETWEEN Diagnosis  ##### 
#Filter the data for diagnosis groups at baseline
physeq_spe <- aggregate_top_taxa2(physeq, level = "Species", top = 30)
physeq_df_s <- microbiomeutilities::phy_to_ldf(physeq_spe, 
                                               transform.counts = "compositional")
# Filter for only baseline
physeq_df_filtered <- physeq_df_s %>%
  filter(Diet %in% c("Baseline"))

# Perform Mann-Whitney U test for each species
results <- physeq_df_filtered %>%
  group_by(Species) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Eligibility_diagnosis)$p.value) %>%
  ungroup()
#Adjust for multiple comparisons using Benjamini-Hochberg method
results <- results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05  # Determine significance
  )
#Identify significant species
significant_species <- results %>%
  filter(significant)
print(significant_species)

# find relative abundance values for signficant species
significant_species_names <- significant_species$Phylum
# Filter the original dataset for significant species
abundance_data <- physeq_df_s %>%
  filter(Phylum %in% significant_species_names) %>%
  group_by(Phylum, Eligibility_diagnosis) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = 'drop')  # Calculate mean abundance
# Print the abundance data
print(abundance_data)


### Top species for each diet   ###.
# Top species in keto diet
physeq_keto <- subset_samples(physeq, Diet == "Keto")
physeq_k_spe <- aggregate_top_taxa2(physeq_keto, level = "Species", top = 11)
physeq_df_ks <- microbiomeutilities::phy_to_ldf(physeq_k_spe, 
                                               transform.counts = "compositional")
write_csv(physeq_df_s, "t2d_top10spe_abundance.csv")

#Create bar plot
ggplot(physeq_df_ks, aes(x = OTUID, y = Abundance)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "Top 10 Most Abundant Species during keto",
       x = "Species",
       y = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Top species for keto diet by diagnosis   ###.
physeq_ketoT <- subset_samples(physeq, Diet == "Keto" & Eligibility_diagnosis == "T2_diabetes")
physeq_k_speT <- aggregate_top_taxa2(physeq_ketoT, level = "Species", top = 11)
physeq_df_ksT <- microbiomeutilities::phy_to_ldf(physeq_k_speT, 
                                                transform.counts = "compositional")
write_csv(physeq_df_ksT, "T2D_top10spe_abundance_keto_only.csv")

#Create bar plot
ggplot(physeq_df_ksT, aes(x = OTUID, y = Abundance)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "Top 10 Most Abundant Species during keto in T2D",
       x = "Species",
       y = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Top species for each diet   ###.
physeq_med <- subset_samples(physeq, Diet == "Med")
physeq_m_spe <- aggregate_top_taxa2(physeq_med, level = "Species", top = 11)
physeq_df_ms <- microbiomeutilities::phy_to_ldf(physeq_m_spe, 
                                                transform.counts = "compositional")
write_csv(physeq_df_ms, "t2d_top10spe_abundance_medonly.csv")

#Create bar plot
ggplot(physeq_df_ms, aes(x = OTUID, y = Abundance)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "Top 10 Most Abundant Species during med",
       x = "Species",
       y = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



####.  PCOA plots baseline data only.  ###
#Perform ordination using Bray-Curtis dissimilarity
physeq_baseline_unrare.ord <- ordinate(physeq_baseline_unrare, "PCoA", "bray")

# Extract sample data from the phyloseq object
sample_data_df <- as.data.frame(sample_data(physeq_baseline_unrare))
# Extract the principal coordinates (PC1 and PC2)
ordination_data <- as.data.frame(physeq_baseline_unrare.ord$vectors)
# Combine sample data with ordination data
combined_data <- cbind(sample_data_df, ordination_data)
# Calculate centroids for each group
centroids <- combined_data %>%
  group_by(Eligibility_diagnosis) %>%
  summarise(PC1 = mean(Axis.1),  # Replace with the actual name of the PC1 column
            PC2 = mean(Axis.2)   # Replace with the actual name of the PC2 column
            )
# View the centroids
print(centroids)

# generate plot
b.div.bray <- plot_ordination(physeq_baseline_unrare, physeq_baseline_unrare.ord, type = "samples", color = "Eligibility_diagnosis") +
  geom_point(size = 3) +  
  scale_color_manual(values = diagnosis_colors) +
  stat_ellipse() +  # Add ellipses
  ggtitle("Bray Curtis - Baseline") +
  #geom_text(aes(label = participant.ID), vjust = -0.5, size = 3) +
  theme_classic() +
  #scale_color_brewer("Diagnosis", palette = "Set2") +  # Color scale for points
# Add centroids with matching colors
  geom_point(data = centroids, 
             aes(x = PC1, y = PC2, color = Eligibility_diagnosis),  # Use color for centroids
             size = 5, 
             shape = 21,  # Use shape 21 for filled circles with outlines
             fill = "#c9daf8ff",  # Fill color for the centroid circles
             stroke = 1.5)  # Adjust stroke width for the outline
# Print the updated plot
print(b.div.bray)

# ck stats on baseline data
otu <- abundances(physeq_baseline_unrare)
meta <- meta(physeq_baseline_unrare)
# adonis2
permanova <- adonis2(t(otu) ~ Eligibility_diagnosis, data = meta, permutations=99, method = "bray")
permanova
capture.output(permanova, file = "permanova.txt")

#Check multivariate dispersion (differences in spread)
bc <- vegan::vegdist(t(otu), method = "bray") 
bd <- vegan::betadisper(bc, meta$Eligibility_diagnosis) 
anova(bd) 
# p > 0.05 suggests dispersion is not the main driver vegan::permutest(bd, permutations = 9999) P was >0.05

###

#####.   PCOA by diet.   #####
#Perform ordination using Bray-Curtis dissimilarity
physeq.ord <- ordinate(physeq, "PCoA", "bray")

# Extract sample data from the phyloseq object
sample_data_df <- as.data.frame(sample_data(physeq))
# Extract the principal coordinates (PC1 and PC2)
ordination_data <- as.data.frame(physeq.ord$vectors)
# Combine sample data with ordination data
combined_data <- cbind(sample_data_df, ordination_data)
# Calculate centroids for each group
centroids <- combined_data %>%
  group_by(Diet) %>%
  summarise(PC1 = mean(Axis.1),  # Replace with the actual name of the PC1 column
            PC2 = mean(Axis.2)   # Replace with the actual name of the PC2 column
  )
# View the centroids
print(centroids)

#Plot the ordination results, comparing diet conditions (bray curtis)
b.div.bray <- plot_ordination(physeq, physeq.ord, type = "samples", color = "Diet") +
  geom_point(size = 3) + stat_ellipse() +
  ggtitle("Bray Curtis - Diet") +
  theme_classic() +
  facet_wrap(~ meta$Eligibility_diagnosis) +
  scale_color_brewer("Diet", palette = "Set2") +
# Add centroids with matching colors
#geom_point(data = centroids, 
#           aes(x = PC1, y = PC2, color = Diet),  # Use color for centroids
#           size = 5, 
#           shape = 21,  # Use shape 21 for filled circles with outlines
#           fill = "black",  # Fill color for the centroid circles
#           stroke = 1.5)  # Adjust stroke width for the outline
print(b.div.bray)

# ck stats on diet compare
otu <- abundances(physeq)
meta <- meta(physeq)
# adonis2
permanova <- adonis2(t(otu) ~ Diet, data = meta, permutations=99, method = "bray")
permanova
capture.output(permanova, file = "permanova_diet.txt")

# Perform pairwise comparisons
distance_matrix_diet <- vegdist(t(otu), method = "bray")
pairwise_results_diet <- pairwise.adonis(distance_matrix_diet, meta$Diet)

# Adjust p-values using Bonferroni correction
p_values <- pairwise_results_diet$p.value
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
results_df <- data.frame(pairwise_results_diet, adjusted_p.value = adjusted_p_values)
print(results_df)

##### ck stats on diet compare by diagnosis ######

#Filter the phyloseq object for prediabetes only
physeq_pre <- subset_samples(physeq, Eligibility_diagnosis == "Prediabetes")
#Perform ordination using Bray-Curtis dissimilarity
physeq_pre.ord <- ordinate(physeq_pre, "PCoA", "bray")

# Extract sample data from the phyloseq object
sample_data_df <- as.data.frame(sample_data(physeq_pre))
# Extract the principal coordinates (PC1 and PC2)
ordination_data <- as.data.frame(physeq_pre.ord$vectors)
# Combine sample data with ordination data
combined_data <- cbind(sample_data_df, ordination_data)
# Calculate centroids for each group
centroids <- combined_data %>%
  group_by(Diet) %>%
  summarise(PC1 = mean(Axis.1),  # Replace with the actual name of the PC1 column
            PC2 = mean(Axis.2)   # Replace with the actual name of the PC2 column
  )
# View the centroids
print(centroids)

# generate plot
b.div.bray <- plot_ordination(physeq_pre, physeq_pre.ord, type = "samples", color = "Diet") +
  geom_point(size = 3) +  
  stat_ellipse() +  # Add ellipses
  ggtitle("Bray Curtis - Prediabetes") +
  #geom_text(aes(label = participant.ID), vjust = -0.5, size = 3) +
  theme_classic() +
  scale_color_brewer("Diet", palette = "Set2") +  # Color scale for points
  # Add centroids with matching colors
  geom_point(data = centroids, 
             aes(x = PC1, y = PC2, color = Diet),  # Use color for centroids
             size = 5, 
             shape = 21,  # Use shape 21 for filled circles with outlines
             fill = "black",  # Fill color for the centroid circles
             stroke = 1.5)  # Adjust stroke width for the outline
# Print the updated plot
print(b.div.bray)

# ck stats on diet compare
otu_pre <- abundances(physeq_pre)
meta_pre <- meta(physeq_pre)
# adonis2
permanova <- adonis2(t(otu_pre) ~ Diet, data = meta_pre, permutations=99, method = "bray")
permanova
capture.output(permanova, file = "permanova_diet_pre.txt")
# Perform pairwise comparisons 
distance_matrix_pre <- vegdist(t(otu_pre), method = "bray")
pairwise_results_pre <- pairwise.adonis(distance_matrix_pre, meta_pre$Diet)
# Adjust p-values using Bonferroni correction
p_values_pre <- pairwise_results_pre$p.value
adjusted_p_values_pre <- p.adjust(p_values_pre, method = "bonferroni")
results_df_pre <- data.frame(pairwise_results_pre, adjusted_p.value = adjusted_p_values_pre)
print(results_df_pre)


#Filter the phyloseq object for T2Ds only
physeq_t2d <- subset_samples(physeq, Eligibility_diagnosis == "T2_diabetes")
#Perform ordination using Bray-Curtis dissimilarity
physeq_t2d.ord <- ordinate(physeq_t2d, "PCoA", "bray")

# Extract sample data from the phyloseq object
sample_data_df <- as.data.frame(sample_data(physeq_t2d))
# Extract the principal coordinates (PC1 and PC2)
ordination_data <- as.data.frame(physeq_t2d.ord$vectors)
# Combine sample data with ordination data
combined_data <- cbind(sample_data_df, ordination_data)
# Calculate centroids for each group
centroids <- combined_data %>%
  group_by(Diet) %>%
  summarise(PC1 = mean(Axis.1),  # Replace with the actual name of the PC1 column
            PC2 = mean(Axis.2)   # Replace with the actual name of the PC2 column
  )
# View the centroids
print(centroids)

# generate plot
b.div.bray <- plot_ordination(physeq_t2d, physeq_t2d.ord, type = "samples", color = "Diet") +
  geom_point(size = 3) +  
  stat_ellipse() +  # Add ellipses
  ggtitle("Bray Curtis - T2D") +
  #geom_text(aes(label = participant.ID), vjust = -0.5, size = 3) +
  theme_classic() +
  scale_color_brewer("Diet", palette = "Set2") +  # Color scale for points
  # Add centroids with matching colors
  geom_point(data = centroids, 
             aes(x = PC1, y = PC2, color = Diet),  # Use color for centroids
             size = 5, 
             shape = 21,  # Use shape 21 for filled circles with outlines
             fill = "black",  # Fill color for the centroid circles
             stroke = 1.5)  # Adjust stroke width for the outline
# Print the updated plot
print(b.div.bray)

# ck stats on diet compare
otu_t2d <- abundances(physeq_t2d)
meta_t2d <- meta(physeq_t2d)
# adonis2
permanova <- adonis2(t(otu_t2d) ~ Diet, data = meta_t2d, permutations=99, method = "bray")
permanova
capture.output(permanova, file = "permanova_diet_t2d.txt")
# Perform pairwise comparisons 
distance_matrix_t2d <- vegdist(t(otu_t2d), method = "bray")
pairwise_results_t2d <- pairwise.adonis(distance_matrix_t2d, meta_t2d$Diet)
# Adjust p-values using Bonferroni correction
p_values_t2d <- pairwise_results_t2d$p.value
adjusted_p_values_t2d <- p.adjust(p_values_t2d, method = "bonferroni")
results_df_t2d <- data.frame(pairwise_results_t2d, adjusted_p.value = adjusted_p_values_t2d)
print(results_df_t2d)

### Compare diagnosis on med diet only
#Perform ordination using Bray-Curtis dissimilarity
physeq_med.ord <- ordinate(physeq_med, "PCoA", "bray")
#Plot the ordination results, comparing diet conditions (bray curtis)
b.div.bray <- plot_ordination(physeq_med, physeq_med.ord, type = "samples", color = "Eligibility_diagnosis") +
  geom_point(size = 3)
#Add ellipses and customize the plot
b.div.bray <- b.div.bray +
  stat_ellipse() +
  ggtitle("Bray Curtis - Med") +
  geom_text(aes(label = participant.ID), vjust = -0.5, size = 3) +
  theme_classic() +
  scale_color_brewer("Eligibility_diagnosis", palette = "Set2")
print(b.div.bray)

# ck stats on diagnosis compare
otu_med <- abundances(physeq_med)
meta_med <- meta(physeq_med)
# adonis2
permanova <- adonis2(t(otu_med) ~ Eligibility_diagnosis, data = meta_med, permutations=99, method = "bray")
permanova
capture.output(permanova, file = "permanova_diet_med.txt")
# Perform pairwise comparisons 
distance_matrix_med <- vegdist(t(otu_med), method = "bray")
pairwise_results_med <- pairwise.adonis(distance_matrix_med, meta_med$Eligibility_diagnosis)
# Adjust p-values using Bonferroni correction
p_values_med <- pairwise_results_med$p.value
adjusted_p_values_med <- p.adjust(p_values_med, method = "bonferroni")
results_df_med <- data.frame(pairwise_results_med, adjusted_p.value = adjusted_p_values_med)
print(results_df_med)

### Compare diagnosis on keto diet only
#Perform ordination using Bray-Curtis dissimilarity
physeq_keto.ord <- ordinate(physeq_keto, "PCoA", "bray")
#Plot the ordination results, comparing diet conditions (bray curtis)
b.div.bray <- plot_ordination(physeq_keto, physeq_keto.ord, type = "samples", color = "Eligibility_diagnosis") +
  geom_point(size = 3)
#Add ellipses and customize the plot
b.div.bray <- b.div.bray +
  stat_ellipse() +
  ggtitle("Bray Curtis - Keto") +
  geom_text(aes(label = participant.ID), vjust = -0.5, size = 3) +
  theme_classic() +
  scale_color_brewer("Eligibility_diagnosis", palette = "Set2")
print(b.div.bray)


###linear mixed effect models
library(lmerTest)   
library(emmeans)
library(dplyr)

mg <- alpha_diversity_combined %>%
  mutate(participant = factor(participant.ID),
         phase       = factor(Diet, levels = c("Baseline", "Keto", "Med")),
         diagnosis   = factor(Eligibility_diagnosis,
                              levels = c("Prediabetes", "T2_diabetes")),
         sequence    = factor(FirstDiet))



metrics <- c("Shannon", "Simpson", "Observed")

## Overall: phase + diagnosis, participant random intercept

for (m in metrics) {
  fit <- lmer(as.formula(paste(m, "~ phase + diagnosis + (1 | participant)")),
              data = mg)

  print(contrast(emmeans(fit, ~ phase), method = "trt.vs.ctrl",
                 ref = "Baseline", adjust = "dunnett"))
}

##  Does the diet effect differ by diagnosis

for (m in metrics) {
  fit <- lmer(as.formula(paste(m, "~ phase * diagnosis + (1 | participant)")),
              data = mg)

  print(anova(fit))
  print(contrast(emmeans(fit, ~ phase | diagnosis), method = "trt.vs.ctrl",
                 ref = "Baseline", adjust = "dunnett"))
}

##Diet sequence / carryover sensitivity 

for (m in metrics) {
  fit <- lmer(as.formula(paste(m, "~ phase * sequence + diagnosis + (1 | participant)")),
              data = mg)
  cat("\n=====", m, "| phase x diet sequence =====\n")
  print(anova(fit))
  print(contrast(emmeans(fit, ~ phase | sequence), method = "trt.vs.ctrl",
                 ref = "Baseline", adjust = "dunnett"))
}



###differential abundance


TOP_SWEEP   <- c(20, 30, 50, 100)
DX_GROUPS   <- c("Prediabetes", "T2_diabetes")
DIET_PHASES <- c("Keto", "Med")
PID_COL     <- "participant.ID"


DROP_UNNAMED_BINS <- FALSE #(species that are named, just not with latin name)

#helper functions

is_unresolved <- function(lab, drop_bins = DROP_UNNAMED_BINS) {
  bad <- is.na(lab) | lab == "" |
    grepl("_unclassified$", lab) |
    grepl("^(unassigned|other)$|incertae.sedis", lab, ignore.case = TRUE)
  if (drop_bins) bad <- bad | grepl("^(GGB|SGB)[0-9]", lab)
  bad
}

## Collapse to a rank and label taxa
ra_at_level <- function(ps, level) {
  ps_lvl <- tax_glom(ps, taxrank = level, NArm = FALSE)
  lab <- as.character(tax_table(ps_lvl)[, level]); lab[is.na(lab)] <- ""
  taxa_names(ps_lvl) <- make.unique(ifelse(lab == "", "Unassigned", lab))
  ps_lvl
}

## Long-format relative abundance at a rank, for plotting.
ra_long <- function(ps, level) {
  psmelt(transform_sample_counts(ra_at_level(ps, level),
                                 function(x) x / sum(x)))
}



##baseline differential abundance by diagnosis



da_baseline <- function(ps, level, top,
                        group       = "Eligibility_diagnosis",
                        groups_keep = DX_GROUPS) {
  
  ps_lvl <- ra_at_level(ps, level)
  lab    <- taxa_names(ps_lvl)
  
  ## relative abundance on the full community, before taxon selection
  ps_ra <- transform_sample_counts(ps_lvl, function(x) x / sum(x))
  
  ## the two diagnosis groups only
  grp_vec <- as.character(sample_data(ps_ra)[[group]])
  ps_ra   <- prune_samples(sample_names(ps_ra)[grp_vec %in% groups_keep], ps_ra)
  
  ## candidate set: resolved taxa, ranked by total abundance, top N
  resolved <- lab[!is_unresolved(lab)]
  ranked   <- sort(taxa_sums(ps_lvl)[resolved], decreasing = TRUE)
  tested   <- names(ranked)[seq_len(min(top, length(ranked)))]
  
  df <- psmelt(prune_taxa(tested, ps_ra))
  df$grp <- factor(df[[group]], levels = groups_keep)
  
  df %>%
    group_by(taxon = OTU) %>%
    summarise(
      n_pre        = sum(grp == groups_keep[1]),
      n_t2d        = sum(grp == groups_keep[2]),
      prev_pre_pct = 100 * mean(Abundance[grp == groups_keep[1]] > 0),
      prev_t2d_pct = 100 * mean(Abundance[grp == groups_keep[2]] > 0),
      mean_pre     = 100 * mean(Abundance[grp == groups_keep[1]]),
      mean_t2d     = 100 * mean(Abundance[grp == groups_keep[2]]),
      med_pre      = 100 * median(Abundance[grp == groups_keep[1]]),
      med_t2d      = 100 * median(Abundance[grp == groups_keep[2]]),
      ## Mann-Whitney U (= Wilcoxon rank-sum), unpaired
      p_value      = tryCatch(suppressWarnings(
        wilcox.test(Abundance ~ grp)$p.value),
        error = function(e) NA_real_),
      .groups = "drop") %>%
    mutate(level    = level,
           top      = top,
           n_tested = n(),
           ## Benjamini-Hochberg across the taxa tested at this rank
           p_BH_within_level = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}

##run both ranks at all four cutoffs 
da_sweep <- bind_rows(lapply(TOP_SWEEP, function(n) {
  bind_rows(da_baseline(physeq_baseline_unrare, "Genus",   top = n),
            da_baseline(physeq_baseline_unrare, "Species", top = n)) %>%
    ## BH pooled across both ranks, within this cutoff
    mutate(p_BH_pooled = p.adjust(p_value, method = "BH"))
}))



padj_matrix <- da_sweep %>%
  select(taxon, level, top, p_BH_within_level) %>%
  pivot_wider(names_from = top, values_from = p_BH_within_level,
              names_prefix = "padj_N") %>%
  left_join(da_sweep %>%
              group_by(taxon, level) %>%
              summarise(across(c(p_value, med_pre, med_t2d,
                                 mean_pre, mean_t2d,
                                 prev_pre_pct, prev_t2d_pct), first),
                        .groups = "drop"),
            by = c("taxon", "level")) %>%
  select(taxon, level, p_value, starts_with("padj_N"),
         med_pre, med_t2d, mean_pre, mean_t2d,
         prev_pre_pct, prev_t2d_pct) %>%
  arrange(level, p_value)

#base R to avoid dplyr error
padj_cols    <- grep("^padj_N", names(padj_matrix), value = TRUE)
min_padj_f2  <- do.call(pmin, c(padj_matrix[padj_cols], na.rm = TRUE))





#######baseline vs each diet delivery phase

da_diet <- function(ps, level, phase,
                    top     = 100,
                    rank_on = c("baseline", "all"),
                    paired  = TRUE,
                    pid     = PID_COL) {
  
  rank_on <- match.arg(rank_on)
  
  ps_lvl <- ra_at_level(ps, level)
  lab    <- taxa_names(ps_lvl)

  ps_ra <- transform_sample_counts(ps_lvl, function(x) x / sum(x))

  diet_vec <- as.character(sample_data(ps_ra)$Diet)
  ps_ra    <- prune_samples(
    sample_names(ps_ra)[diet_vec %in% c("Baseline", phase)], ps_ra)

  
  resolved <- lab[!is_unresolved(lab)]
  rank_ps  <- if (rank_on == "baseline") {
    prune_samples(sample_names(ps_lvl)[
      as.character(sample_data(ps_lvl)$Diet) == "Baseline"], ps_lvl)
  } else ps_lvl
  ranked <- sort(taxa_sums(rank_ps)[resolved], decreasing = TRUE)
  tested <- names(ranked)[seq_len(min(top, length(ranked)))]
  
  df <- psmelt(prune_taxa(tested, ps_ra))
  df$phase_f <- factor(df$Diet, levels = c("Baseline", phase))
  df$pid     <- as.character(df[[pid]])
  
  df %>%
    group_by(taxon = OTU) %>%
    group_modify(function(d, key) {
      
      bl <- d[d$phase_f == "Baseline", ]
      on <- d[d$phase_f == phase, ]
      
      
      if (paired) {
        ids <- intersect(names(bl), names(on))
        p <- if (length(ids) >= 3) {
          tryCatch(suppressWarnings(
            wilcox.test(on[ids], bl[ids], paired = TRUE)$p.value),
            error = function(e) NA_real_)
        } else NA_real_
        n_used     <- length(ids)
        med_change <- 100 * median(on[ids] - bl[ids])
      } else {
        p <- tryCatch(suppressWarnings(
          wilcox.test(d$Abundance ~ d$phase_f)$p.value),
          error = function(e) NA_real_)
        n_used     <- NA_integer_
        med_change <- NA_real_
      }
      
      data.frame(
        n_baseline        = length(bl),
        n_phase           = length(on),
        n_pairs           = n_used,
        med_baseline      = 100 * median(bl),
        med_phase         = 100 * median(on),
        mean_baseline     = 100 * mean(bl),
        mean_phase        = 100 * mean(on),
        med_paired_change = med_change,
        p_value           = p)
    }) %>%
    ungroup() %>%
    mutate(level    = level,
           phase    = phase,
           test     = ifelse(paired, "paired_signedrank", "unpaired_MWU"),
           rank_on  = rank_on,
           top      = top,
           n_tested = n(),
           p_BH     = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}

##both phases, both ranks, four cutoffs, both tests
grid4 <- expand.grid(phase  = DIET_PHASES,
                     level  = c("Genus", "Species"),
                     top    = TOP_SWEEP,
                     paired = c(TRUE, FALSE),
                     stringsAsFactors = FALSE)

da_diet_all <- bind_rows(lapply(seq_len(nrow(grid4)), function(i)
  da_diet(physeq,
          level   = grid4$level[i],
          phase   = grid4$phase[i],
          top     = grid4$top[i],
          rank_on = "baseline",
          paired  = grid4$paired[i])))



##what is significant
diet_summary <- da_diet_all %>%
  group_by(phase, level, test, top) %>%
  summarise(n_tested = first(n_tested),
            n_sig    = sum(p_BH < 0.05, na.rm = TRUE),
            sig_taxa = paste(taxon[which(p_BH < 0.05)], collapse = "; "),
            .groups  = "drop") %>%
  arrange(phase, level, test, top)
print(as.data.frame(diet_summary))
write.csv(diet_summary, "fig4_DA_diet_summary.csv", row.names = FALSE)


## padj at every cutoff, one row per taxon x phase x rank x test
diet_padj <- da_diet_all %>%
  select(taxon, level, phase, test, top, p_BH) %>%
  pivot_wider(names_from = top, values_from = p_BH,
              names_prefix = "padj_N") %>%
  left_join(da_diet_all %>%
              group_by(taxon, level, phase, test) %>%
              summarise(across(c(p_value, n_pairs,
                                 med_baseline, med_phase,
                                 mean_baseline, mean_phase,
                                 med_paired_change), first),
                        .groups = "drop"),
            by = c("taxon", "level", "phase", "test")) %>%
  arrange(phase, level, test, p_value)

padj_cols_f4 <- grep("^padj_N", names(diet_padj), value = TRUE)
min_padj_f4  <- do.call(pmin, c(diet_padj[padj_cols_f4], na.rm = TRUE))






