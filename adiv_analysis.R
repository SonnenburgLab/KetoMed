# Load libraries
library(ape) # Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq
library(ggplot2) # Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.
library(phyloseq) # The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.
library(plotly) # A package to create interactive web graphics of use in 3D plots
library(vegan) # The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.
#library(philr) # This package provides functions for the analysis of compositional data 
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
library(microbiome) # Data analysis and visualization
library(ggplotify) # Data visualization
library(magrittr)

# read in phyloseq object and filter 
physeq <- readRDS("/Users/samgoldman/Documents/keto_med/from_erica/dada2_phyloseq.rds")
# Filtering out ASVs with less than 10 reads across all phyloseq object
physeq <- prune_taxa(taxa_sums(physeq) > 10, physeq)
# filtering out ASVs not present in more than 5 samples
physeq <- filter_taxa(physeq, function(x) mean(x) > 5, TRUE)
# Filtering out samples with no reads
physeq <- prune_samples(sample_sums(physeq) > 0, physeq)

# Load sample metadata from CSV file using sailen's metadata file
sample_metadata <- read.csv("/Users/samgoldman/Documents/keto_med/from_erica/ketomed_meta_sailen1.csv", header = TRUE, row.names = 1) 

# Transfer metadata to the phyloseq object
# Create a phyloseq sample data object
sample_data_physeq1 <- sample_data(sample_metadata)
# Assign this sample data to existing phyloseq object
sample_data(physeq) <- sample_data_physeq1
# Verify the sample data has been added
print(sample_data(physeq))

# distribution plot of filtered reads to determine how to rarefy
# Extract total reads per sample using sample_sums
total_reads <- sample_sums(physeq)
# Create a data frame for plotting
total_reads_df <- data.frame(Sample = names(total_reads), TotalReads = total_reads)
# Visualize the distribution of total read counts
read_distribution_plot <- ggplot(total_reads_df, aes(x = TotalReads)) +
  geom_histogram(binwidth = 1000, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Total Reads per Sample", 
       x = "Total Reads", 
       y = "Number of Samples") +
  theme_minimal()
# Print the plot
print(read_distribution_plot)
write.csv(total_reads_df, file='read_dist.csv')

# Rarefy the data to 5,000 reads
set.seed(123)   #essentially does nothing since rngseed sets seed to 1
rarefied_physeq <- rarefy_even_depth(physeq, sample.size = 5000, rngseed = TRUE)

# Transfer metadata to the rarefied phyloseq object
# Create a phyloseq sample data object
sample_data_physeq <- sample_data(sample_metadata)
# Assign this sample data to existing phyloseq object
sample_data(rarefied_physeq) <- sample_data_physeq
# Verify the sample data has been added
print(sample_data(rarefied_physeq))


# Calculate alpha diversity metrics
alpha_diversity <- estimate_richness(rarefied_physeq)  # Estimates several metrics
# Convert alpha diversity metrics to data frame
alpha_diversity_df <- data.frame(alpha_diversity)
# Combine alpha diversity data with sample metadata
alpha_diversity_combined <- cbind(data.frame(sample_data(rarefied_physeq)), alpha_diversity_df)
# Check the combined dataset
head(alpha_diversity_combined)

# read in alpha div csv with averaged baselines
a_div_avebl <- read.csv("alpha_rare5k_blave.csv")

# set custom colors and order of graph
custom_colors <- c("baseline" = "#b7b7b7ff","Keto-halfway" = "#c7885fff", "Keto-complete" = "#f5c09fff", "Med-halfway" = "#669289ff", "Med-complete" = "#afd3caff")

a_div_avebl$Group <- factor(a_div_avebl$Group, 
                            levels = c("baseline", "Keto-halfway", "Keto-complete", "Med-halfway", "Med-complete"))

# Visualize alpha diversity by diagnosis in box plots
ggplot(a_div_avebl, aes(x = Group, y = Simpson, fill = Group)) +  # Replace SampleType with actual treatment column
  scale_fill_manual(values = custom_colors) +
  geom_violin(trim = TRUE, alpha = 0.5) +  # Fill color based on treatment group
  geom_boxplot(width = 0.4, alpha = 0.5) +
  geom_jitter(width = 0.1, size = 1, color = "black", alpha = 0.3) +
  labs(title = "Alpha Diversity by Treatment", 
       x = "Group", 
       y = "Simpson") +
  #facet_wrap(~ Eligibility_diagnosis) +  # Create separate panels for each diagnosis
  theme_minimal()

## bray curtis distance comparisons ###

### Baseline Diagnosis Analysis
#create baseline physeq of unrarefied data
physeq_baseline_unrare <- subset_samples(physeq, Group == "baseline")
#define custom colors
diagnosis_colors <- c("T2_diabetes" = "#741b47ff","Prediabetes" = "#d5a6bdff")

# Calculate Bray-Curtis distance
physeq_bl.ord <- ordinate(physeq_baseline_unrare, "PCoA", "bray")
# Extract sample data from the phyloseq object
sample_data_bl <- as.data.frame(sample_data(physeq_baseline_unrare))
# Extract the principal coordinates (PC1 and PC2)
ordination_data_bl <- as.data.frame(physeq_bl.ord$vectors)
# Combine sample data with ordination data
combined_data_bl <- cbind(sample_data_bl, ordination_data_bl)

#Calculate centroids
centroids_bl <- combined_data_bl %>%
  group_by(Eligibility_diagnosis) %>%
  summarise(PC1 = mean(Axis.1),  # Replace with the actual name of the PC1 column
            PC2 = mean(Axis.2)   # Replace with the actual name of the PC2 column
  )
# View the centroids
print(centroids_bl)

# generate plot
b.div.bray <- plot_ordination(physeq_baseline_unrare, physeq_bl.ord, type = "samples", color = "Eligibility_diagnosis") +
  geom_point(size = 3) +  
  scale_color_manual(values = diagnosis_colors) +
  stat_ellipse() +  # Add ellipses
  ggtitle("Bray Curtis") +
  theme_classic() +
  # Add centroids with matching colors
  geom_point(data = centroids_bl, 
             aes(x = PC1, y = PC2, color = Eligibility_diagnosis),  # Use color for centroids
             size = 5, 
             shape = 21,  # Use shape 21 for filled circles with outlines
             fill = "black",
             #fill = "#c9daf8ff",  # Fill color for the centroid circles
             stroke = 1.5)  # Adjust stroke width for the outline
# Print the updated plot
print(b.div.bray)

######################

#Perform ordination using Bray-Curtis dissimilarity of unrarefied 16s data for all timepoints
physeq.ord <- ordinate(physeq, "PCoA", "bray")

# Extract sample data from the phyloseq object
sample_data_df <- as.data.frame(sample_data(physeq))
# Extract the principal coordinates (PC1 and PC2)
ordination_data <- as.data.frame(physeq.ord$vectors)
# Combine sample data with ordination data
combined_data <- cbind(sample_data_df, ordination_data)

# Calculate centroids for each group
centroids <- combined_data %>%
  group_by(Group2) %>%
  summarise(PC1 = mean(Axis.1),  # Replace with the actual name of the PC1 column
            PC2 = mean(Axis.2)   # Replace with the actual name of the PC2 column
  )
# View the centroids
print(centroids)



custom_colors2 <- c("B" = "#b7b7b7ff","KD" = "#c7885fff", "KO" = "#f5c09fff", "MD" = "#669289ff", "MO" = "#afd3caff")


# generate plot
b.div.bray <- plot_ordination(physeq, physeq.ord, type = "samples", color = "Group2") +
  geom_point(size = 3) +  
  scale_color_manual(values = custom_colors2) +
  stat_ellipse() +  # Add ellipses
  ggtitle("Bray Curtis") +
  #geom_text(aes(label = participant.ID), vjust = -0.5, size = 3) +
  theme_classic() +
  #scale_color_brewer("Diagnosis", palette = "Set2") +  # Color scale for points
  # Add centroids with matching colors
  geom_point(data = centroids, 
             aes(x = PC1, y = PC2, color = Group2),  # Use color for centroids
             size = 5, 
             shape = 21,  # Use shape 21 for filled circles with outlines
             fill = "black",
             #fill = "#c9daf8ff",  # Fill color for the centroid circles
             stroke = 1.5)  # Adjust stroke width for the outline
# Print the updated plot
print(b.div.bray)

# ck stats 
otu <- abundances(physeq)
meta <- meta(physeq)
# adonis2
permanova <- adonis2(t(otu) ~ Group2, data = meta, permutations=99, method = "bray")
permanova
capture.output(permanova, file = "permanova.txt")

# Perform pairwise comparisons
distance_matrix_diet <- vegdist(t(otu), method = "bray")
pairwise_results_diet <- pairwise.adonis(distance_matrix_diet, meta$Group2)

# Adjust p-values using Bonferroni correction
p_values <- pairwise_results_diet$p.value
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
results_df <- data.frame(pairwise_results_diet, adjusted_p.value = adjusted_p_values)
print(results_df)



################## get bray curtis distances ################
# Calculate the Bray-Curtis distance matrix
bray_curtis_dist <- vegdist(t(otu), method = "bray")
# Convert the distance matrix to a data frame for visualization
bray_curtis_df <- as.data.frame(as.matrix(bray_curtis_dist))
# Add metadata to the distance data frame
# Assuming the row names of bray_curtis_df match the SampleID in metadata
bray_curtis_df$sample_id <- rownames(bray_curtis_df)
merged_data <- merge(bray_curtis_df, sample_data_df, by = rownames)


# Convert the distance matrix to long format
bray_curtis_long <- melt(as.matrix(bray_curtis_dist), varnames = c("Sample1", "Sample2"), value.name = "Distance")
write_csv(bray_curtis_long, "bray_dist_16s_KM.csv")

# read in analyzed bray curtis distance csv
bc_dist <- read_csv("bc_dist.csv")

custom_colors3 <- c("within" = "#b7b7b7ff","bet_KD" = "#c7885fff", "bet_KO" = "#f5c09fff", "bet_MD" = "#669289ff", "bet_MO" = "#afd3caff")

# Create the bar plot with error bars
ggplot(bc_dist, aes(x = arm, y = Overall_Average)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
  geom_errorbar(aes(ymin = Overall_Average - SD_Average, 
                    ymax = Overall_Average + SD_Average), 
                width = 0.2) +
  labs(title = "Average Bray-Curtis Distance by Arm at P4",
       x = "Arm",
       y = "Average Bray-Curtis Distance") +
  theme_minimal()



# looking for ASVs lost in the keto-own phase from baseline
get_asvs_in_both <- function(physeq, output_file) {
  # Initialize a data frame to store results
  results <- data.frame()
  
  # Get unique participant IDs
  participant_ids <- unique(sample_data(physeq)$participant_id)
  
  
  # taxonomy work
  # Calculate relative abundance
  physeq_rel_abund <- transform_sample_counts(physeq, function(x) x / sum(x))
  # Extract the relative abundance data
  rel_abund_data <- psmelt(physeq_rel_abund)
  # Create the plot
  ggplot(rel_abund_data, aes(x = Sample, y = Abundance, fill = Species)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Group2, scales = "free_x") +  
    labs(title = "Relative Abundance of Taxa by Group",
         x = "Group2",
         y = "Relative Abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Group by treatment and taxon, then calculate the average relative abundance
  avg_rel_abund <- rel_abund_data %>%
    group_by(Group2, Family) %>%  
    summarise(AvgAbundance = mean(Abundance), .groups = 'drop')
  
  
  physeq_fam <- microbiome::aggregate_rare(physeq, level = "Family", detection = 50/100, prevalence = 70/100)
  
  # Alternative
  physeq_fam <- aggregate_top_taxa2(physeq, level = "Family", top = 25)
  
  physeq.fam.rel <- microbiome::transform(physeq_fam, "compositional")
  
  physeq.fam.rel <- physeq %>%
    aggregate_rare(level = "Family", detection = 50/100, prevalence = 70/100) %>%
    microbiome::transform(transform = "compositional")
  
  plot_composition(physeq.fam.rel,sample.sort = "BodySite", x.label = "SampleID") + theme(legend.position = "bottom") + scale_fill_brewer("Family", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))
  
  # b diversity PCOA plots
  dist = phyloseq::distance(physeq, method="bray")
  ordination = ordinate(physeq, method="PCoA", distance=dist)
  #Generate plot
  plot_ordination(physeq, ordination, color="participant_id") + 
    theme_classic() +
    theme(strip.background = element_blank())
  
  # Only look at baseline by diagnosis
  desired_populations <- c("T2_diabetes", "Prediabetes")  
  desired_time_points <- c("B1", "B2")
  # Filter the phyloseq object
  physeq_filtered_base <- subset_samples(physeq, 
                                         Eligibility_diagnosis %in% desired_populations & 
                                           Group1 %in% desired_time_points)
  # Calculate Bray-Curtis distance
  dist = phyloseq::distance(physeq_filtered_base, method = "bray")
  ordination = ordinate(physeq_filtered_base, method = "PCoA", distance = dist)
  #Generate plot
  plot_ordination(physeq_filtered_base, ordination, color = "participant_id") + 
    theme_classic() +
    theme(strip.background = element_blank())
  
  #####  TAXONOMIC ANALYSIS BETWEEN Diagnosis  ##### 
  
  physeq_fam <- aggregate_top_taxa2(physeq, level = "Genus", top = 15)
  physeq_df_s <- microbiomeutilities::phy_to_ldf(physeq_fam, 
                                                 transform.counts = "compositional")
  # Filter for only baseline
  physeq_df_filtered <- physeq_df_s %>%
    filter(Diet %in% c("Baseline"))
  
  # Perform Mann-Whitney U test for each species
  results <- physeq_df_filtered %>%
    group_by(Genus) %>%
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
  
  
  
  
  
  
  
  
  ### linear mixed models ###
  
  library(lmerTest)     # lmer() with Satterthwaite p-values
  library(emmeans)      # baseline (Dunnett) contrasts
  
  #build modeling frame
  adiv_lmm <- alpha_diversity_combined
  adiv_lmm$phase <- dplyr::case_when(
    adiv_lmm$Group1 %in% c("B1", "B2") ~ "Baseline",
    adiv_lmm$Group1 == "KD" ~ "WFKD_del",
    adiv_lmm$Group1 == "KO" ~ "WFKD_SP",
    adiv_lmm$Group1 == "MD" ~ "Med_del",
    adiv_lmm$Group1 == "MO" ~ "Med_SP"
  )
  adiv_lmm$diagnosis <- factor(adiv_lmm$Eligibility_diagnosis,
                               levels = c("Prediabetes", "T2_diabetes"))
  adiv_lmm$order  <- factor(adiv_lmm$FirstDiet)
  
  ##average two baseline samples
  dat_avg <- adiv_lmm %>%
    dplyr::group_by(participant_id, phase, diagnosis, order) %>%
    dplyr::summarise(dplyr::across(c(Observed, Shannon, Simpson),
                                   ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::mutate(
      phase = factor(phase, levels = c("Baseline", "WFKD_del", "WFKD_SP",
                                       "Med_del", "Med_SP")),   # Baseline = reference
      participant_id = factor(participant_id))
  
  cat("\nRows per phase after collapsing:\n"); print(table(dat_avg$phase))
  
  metrics <- c("Observed", "Shannon", "Simpson")
  
  #overall, each phase vs baseline
  fit_overall <- function(metric, data) {
    f   <- as.formula(paste0(metric, " ~ phase + diagnosis + (1 | participant_id)"))
    mod <- lmer(f, data = data)
    ct  <- contrast(emmeans(mod, ~ phase),
                    method = "trt.vs.ctrl", ref = "Baseline", adjust = "dunnettx")
    data.frame(metric = metric, model = "overall", as.data.frame(ct))
  }
  overall_tab <- do.call(rbind, lapply(metrics, fit_overall, data = dat_avg))
  cat("LMM contrasts vs baseline \n")
  print(overall_tab, digits = 3)
  
  ##by diagnosis, phase vs baseline

  fit_by_dx <- function(metric, data) {
    f   <- as.formula(paste0(metric, " ~ phase * diagnosis + (1 | participant_id)"))
    mod <- lmer(f, data = data)
    ct  <- contrast(emmeans(mod, ~ phase | diagnosis),
                    method = "trt.vs.ctrl", ref = "Baseline", adjust = "dunnettx")
    data.frame(metric = metric, model = "by_diagnosis", as.data.frame(ct))
  }
  bydx_tab <- do.call(rbind, lapply(metrics, fit_by_dx, data = dat_avg))
  cat("\ nby diagnosis LMM contrasts vs baseline =====\n")
  print(bydx_tab, digits = 3)
  
  
  
  ## beta div
  library(vegan); library(permute)
  
  bray_np <- phyloseq::distance(physeq, method = "bray")   # unrarefied, 203 samples
  meta_np <- data.frame(sample_data(physeq))
  meta_np <- meta_np[labels(bray_np), , drop = FALSE]
  meta_np$Group2         <- factor(meta_np$Group2)
  meta_np$participant_id <- factor(meta_np$participant_id)
  
  ## Fig 3C: individuality (composition ~ participant)
  set.seed(123)
  perm_individual <- adonis2(bray_np ~ participant_id, data = meta_np, permutations = 999)
  cat(" composition ~ participant\n"); print(perm_individual)
  
  bm <- as.matrix(bray_np); ids <- meta_np$participant_id
  same <- outer(ids, ids, "=="); upr <- upper.tri(bm)
  cat(sprintf("Mean within-participant Bray:  %.3f (%d pairs)\n",
              mean(bm[upr & same]),  sum(upr & same)))
  cat(sprintf("Mean between-participant Bray: %.3f (%d pairs)\n",
              mean(bm[upr & !same]), sum(upr & !same)))
  
  ## within-participant restricted permutation
  ctrl <- how(blocks = meta_np$participant_id, nperm = 999)
  set.seed(123)
  permanova_phase <- adonis2(bray_np ~ Group2, data = meta_np, permutations = ctrl)
  cat("\n composition ~ phase, within-participant perms =====\n")
  print(permanova_phase)
  
  ## dispersion check (location vs spread)
  set.seed(123)
  bd_test <- permutest(betadisper(bray_np, meta_np$Group2), permutations = 999)
  cat("\n dispersion (betadisper) =====\n"); print(bd_test)
  
  # Pairwise phase comparisons, within-participant perm

  phases <- levels(meta_np$Group2)
  pairs_df <- do.call(rbind, combn(phases, 2, simplify = FALSE, FUN = function(pr) {
    idx  <- which(meta_np$Group2 %in% pr)
    d    <- as.dist(as.matrix(bray_np)[idx, idx])
    md   <- droplevels(meta_np[idx, ])
    set.seed(123)
    a <- adonis2(d ~ Group2, data = md,
                 permutations = how(blocks = factor(md$participant_id), nperm = 999))
    data.frame(pair = paste(pr, collapse = " vs "),
               R2 = a$R2[1], F = a$F[1], p = a$`Pr(>F)`[1])
  }))
  pairs_df$p_bonferroni <- p.adjust(pairs_df$p, method = "bonferroni")
  cat("\npairwise phase comparisons (within participant)")
  print(pairs_df, row.names = FALSE, digits = 3)
  
  
  ## Diet order 
  fit_order <- function(metric, data) {
    mod <- lmerTest::lmer(as.formula(paste0(metric,
                                            " ~ phase * order + diagnosis + (1 | participant_id)")),
                          data = data)
    aov_tab <- as.data.frame(anova(mod))
    if (!"NumDF" %in% names(aov_tab))
      stop("anova() has no df/p-value columns; model is not an lmerTest fit.")
    row <- grep("phase:order|order:phase", rownames(aov_tab))
    a   <- aov_tab[row, ]
    omn <- data.frame(metric = metric, NumDF = a[["NumDF"]], DenDF = a[["DenDF"]],
                      F = a[["F value"]], p.value = a[["Pr(>F)"]])
    within <- contrast(emmeans(mod, ~ phase | order),
                       method = "trt.vs.ctrl", ref = "Baseline", adjust = "dunnettx")
    did <- contrast(emmeans(mod, ~ phase * order),
                    interaction = c(phase = "trt.vs.ctrl", order = "pairwise"))
    list(omnibus = omn,
         within  = data.frame(metric = metric, as.data.frame(within)),
         did     = data.frame(metric = metric, as.data.frame(did)))
  }
  
  order_res <- lapply(metrics, fit_order, data = dat_avg)
  order_omnibus <- do.call(rbind, lapply(order_res, `[[`, "omnibus"))
  order_within  <- do.call(rbind, lapply(order_res, `[[`, "within"))
  order_did     <- do.call(rbind, lapply(order_res, `[[`, "did"))
  
  
  ######sp vs del
  
  
  metrics <- c("Observed", "Shannon", "Simpson")
  

  PHLV  <- levels(dat_avg$phase)
  spdel <- list(
    "Keto: SP - del" = as.integer(PHLV == "WFKD_SP") - as.integer(PHLV == "WFKD_del"),
    "Med: SP - del"  = as.integer(PHLV == "Med_SP")  - as.integer(PHLV == "Med_del"))
  
  ## OVERALL: SP vs del, adjusted for diagnosis
  spdel_fit_overall <- function(metric, data) {
    mod <- lmerTest::lmer(as.formula(paste0(metric, " ~ phase + diagnosis + (1 | participant_id)")),
                          data = data)
    ct  <- contrast(emmeans(mod, ~ phase), method = spdel, adjust = "bonferroni")
    data.frame(metric = metric, model = "overall", as.data.frame(ct))
  }
  spdel_overall_tab <- do.call(rbind, lapply(metrics, spdel_fit_overall, data = dat_avg))
  cat("\n===== SP - del: overall (Bonferroni across the 2 diets) =====\n")
  print(spdel_overall_tab, digits = 3)
  
  ##order
  spdel_diet_only <- droplevels(subset(dat_avg, phase != "Baseline"))
  spdel_carryover_tab <- do.call(rbind, lapply(metrics, function(m) {
    mod <- lmerTest::lmer(as.formula(paste0(m, " ~ phase + diagnosis + order + (1 | participant_id)")),
                          data = spdel_diet_only)
    a <- as.data.frame(anova(mod))["order", ]
    data.frame(metric = m, F = a[["F value"]], p.value = a[["Pr(>F)"]])
  }))
  cat("\n===== order (carryover) test =====\n")
  print(spdel_carryover_tab, digits = 3)
  
  ##  Diet-order interaction
  
  spdel_fit_order <- function(metric, data) {
    mod <- lmerTest::lmer(as.formula(paste0(metric,
                                            " ~ phase * order + diagnosis + (1 | participant_id)")), data = data)
    aov_tab <- as.data.frame(anova(mod))
    a   <- aov_tab[grep("phase:order|order:phase", rownames(aov_tab)), ]
    omn <- data.frame(metric = metric, NumDF = a[["NumDF"]], DenDF = a[["DenDF"]],
                      F = a[["F value"]], p.value = a[["Pr(>F)"]])
    within <- contrast(emmeans(mod, ~ phase | order), method = spdel, adjust = "bonferroni")
    diff   <- as.data.frame(contrast(within, method = "pairwise", by = "contrast", adjust = "none"))
    diff$p.bonferroni <- pmin(1, diff$p.value * nrow(diff))
    list(omnibus = omn,
         within  = data.frame(metric = metric, as.data.frame(within)),
         diff    = data.frame(metric = metric, diff))
  }
  spdel_order_res     <- lapply(metrics, spdel_fit_order, data = dat_avg)
  spdel_order_omnibus <- do.call(rbind, lapply(spdel_order_res, `[[`, "omnibus"))
  spdel_order_within  <- do.call(rbind, lapply(spdel_order_res, `[[`, "within"))
  spdel_order_diff    <- do.call(rbind, lapply(spdel_order_res, `[[`, "diff"))
 
  #paired t test sp vs del 
  spdel_paired <- function(metric, data) {
    wide <- tidyr::pivot_wider(data[, c("participant_id", "phase", metric)],
                               names_from = phase, values_from = dplyr::all_of(metric))
    do.call(rbind, lapply(list(c("Keto", "WFKD_SP", "WFKD_del"),
                               c("Med",  "Med_SP",  "Med_del")), function(p) {
                                 ok <- complete.cases(wide[[p[2]]], wide[[p[3]]])
                                 tt <- t.test(wide[[p[2]]][ok], wide[[p[3]]][ok], paired = TRUE)   # SP - del
                                 data.frame(metric = metric, diet = p[1], contrast = "SP - del",
                                            n_pairs = sum(ok), p_paired_t = tt$p.value)
                               }))
  }
  spdel_paired_tab <- do.call(rbind, lapply(metrics, spdel_paired, data = dat_avg))
  spdel_paired_tab$p_bonferroni <- ave(spdel_paired_tab$p_paired_t, spdel_paired_tab$metric,
                                       FUN = function(p) pmin(1, p * length(p)))



  
  
  ###longitudinal change in shannon split by diet order
  p_to_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                                            ifelse(p < 0.05, "*", "")))
  
  violin_by_order <- function(seq_group) {
    lv <- chrono_order[[seq_group]]
    d  <- dat_avg[dat_avg$order == seq_group, ]
    d$phase_chrono <- factor(as.character(d$phase), levels = lv)
    
    ## significant del -> SP comparisons only
    ann <- data.frame(diet = c("WFKD", "Med"),
                      x1   = c(match("WFKD_del", lv), match("Med_del", lv)),
                      x2   = c(match("WFKD_SP",  lv), match("Med_SP",  lv)))
    ann$p     <- spdel_sh$p.value[match(paste(seq_group, ann$diet),
                                        paste(spdel_sh$order, spdel_sh$contrast))]
    ann$stars <- p_to_stars(ann$p)
    ann       <- ann[ann$stars != "", ]
    span  <- diff(y_range)
    ann$y <- rep(y_range[2] + 0.05 * span, nrow(ann))
    
    ggplot(d, aes(x = phase_chrono, y = Shannon, fill = phase)) +
      scale_fill_manual(values = phase_cols) +
      geom_violin(trim = TRUE, alpha = 0.5) +
      geom_boxplot(width = 0.4, alpha = 0.5, outlier.shape = NA) +
      geom_jitter(width = 0.1, size = 1, colour = "black", alpha = 0.3) +
      geom_segment(data = ann, aes(x = x1, xend = x2, y = y, yend = y),
                   inherit.aes = FALSE) +
      geom_text(data = ann, aes(x = (x1 + x2) / 2, y = y + 0.02 * span, label = stars),
                inherit.aes = FALSE, size = 6) +
      scale_x_discrete(labels = phase_labels) +
      coord_cartesian(ylim = c(y_range[1], y_range[2] + 0.12 * span)) +
      labs(title = seq_group, x = NULL, y = "Shannon") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            axis.title.y = element_text(margin = margin(r = 8)),
            plot.margin  = margin(10, 20, 10, 10))
  }
  
  p_kf <- violin_by_order("Keto first")
  p_mf <- violin_by_order("Med first")
  p_kf + p_mf + plot_layout(widths = c(1, 1))
  