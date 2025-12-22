### Analysis on bacterial species level ###
# ---- a-diversity: all samples combined to one boxplot ----

# load packages
library(tidyverse)
library(ggplot2)
library(dplyr)

# read CSV
df <- read.csv("species_table_clean_final_binary.csv", header = TRUE)

# select columns: only numeric, but without 'total'
num_cols <- sapply(df, is.numeric)
num_cols["total"] <- FALSE
df_num <- df[, num_cols]

# Species richness for all samples combined
richness <- colSums(df_num)

richness_df <- data.frame(
  Sample = names(richness),
  Richness = as.numeric(richness)
)

# define category-columns
richness_df$Category <- ifelse(grepl("PC_sample", richness_df$Sample), "PC", "Sample")

# boxplot: all samples combined
ggplot(richness_df, aes(x = Category, y = Richness, fill = Category)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("Sample" = "#d95f02", "PC" = "#1b9e77")) +
  labs(title = NULL, x = NULL, y = "number of detected bacteria species") +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )

# ---- a-diversity per nematode species ---- 

sample_names <- colnames(df_num)

# Assign species to column names (adapted to naming convention)
sample_to_species <- sapply(sample_names, function(x) {
  if(grepl("Acrobeles_emmatus", x)) return("A. emmatus")
  if(grepl("Acrobeles_mariannae", x)) return("A. mariannae")
  if(grepl("Acrobeles_sp", x)) return("Acrobeles sp.")
  if(grepl("Acrobeloides", x)) return("Acrobeloides")
  if(grepl("Cervidellus", x)) return("Cervidellus")
  if(grepl("Mesorhabditis", x)) return("Mesorhabditis")
  if(grepl("Nothacrobeles", x)) return("Nothacrobeles")
  if(grepl("Panagrolaimus", x)) return("Panagrolaimus")
  if(grepl("PC", x)) return("PC")
  if(grepl("Stegelleta", x)) return("Stegelleta")
  if(grepl("Stegelletina", x)) return("Stegelletina")
  return(NA)
})

species_long <- df_num %>%
  rownames_to_column(var = "Bacteria") %>%
  pivot_longer(-Bacteria, names_to = "Sample", values_to = "Count") %>%
  mutate(Art = sample_to_species[Sample])



# number of different bacterial species per nematode species
unique_species_per_art <- species_long %>%
  filter(!Art %in% c("Nothacrobeles", "Panagrolaimus")) %>%  
  group_by(Art, Bacteria) %>%
  summarise(n_occurrences = sum(Count > 0), .groups = "drop") %>% 
  filter(n_occurrences == 1) %>%          
  group_by(Art) %>%
  summarise(UniqueSpecies = n(), .groups = "drop")


# species richness
richness_different <- species_long %>%
  group_by(Sample, Art) %>%
  summarise(Richness = sum(Count, na.rm = TRUE), .groups = "drop")


richness_different$Art <- factor(richness_different$Art, 
                                 levels = c("PC", "Cervidellus", "Acrobeles sp.", "Mesorhabditis", "Acrobeloides",  
                                             "A. emmatus", "A. mariannae", "Stegelleta", "Stegelletina"))
                                              

unique_species_per_art$Art <- factor(unique_species_per_art$Art,
                                     levels = levels(richness_different$Art))

# Boxplot
ggplot(richness_different, aes(x = Art, y = Richness, fill = Art)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  geom_text(data = unique_species_per_art, 
            aes(x = Art, y = max(richness_different$Richness)*1.5, 
                label = paste0("n=", UniqueSpecies)), 
            inherit.aes = FALSE, size = 5) +
  labs(x = NULL, y = "Number of detected bacterial species") +
  scale_y_log10() +
  scale_fill_manual(values = c(
    "PC" = "grey80",
    "Cervidellus"   = "grey80",
    "Acrobeloides"  = "grey80",
    "Stegelletina"  = "grey80",
    "Acrobeles sp." = "grey80",
    "A. emmatus"    = "grey80",
    "A. mariannae"  = "grey80",
    "Mesorhabditis" = "grey80",
    "Nothacrobeles" = "grey80",
    "Panagrolaimus" = "grey80",
    "Stegelleta"    = "grey80"
  )) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.y = element_text(size = 15),
    
    axis.line.x = element_line(color = "black", size = 0.8),  # X-Achse Linie
    axis.line.y = element_line(color = "black", size = 0.8),
    
    axis.ticks.x = element_line(color = "black", size = 0.8), # X-Achsen-Ticks
    axis.ticks.y = element_line(color = "black", size = 0.8),
    
    axis.ticks.length = unit(0.25, "cm")
  )




  


# ---- a-diversity: statistics (Shapiro, Kruskal-Wallis, pairwise wilcoxon) ----

# test for normal distribution of data (p<0.05 = not normally distributed -> non-parametric tests needed)
shapiro_species <- richness_different %>%
  group_by(Art) %>%
  summarise(
    W = if(n() >= 3) shapiro.test(Richness)$statistic else NA,
    p_value = if(n() >= 3) shapiro.test(Richness)$p.value else NA,
    n = n()
  )

# display the result
shapiro_species

# since there are groups with <3 samples, and also the other groups only have max. 6 samples,
# it is very advisable to use a non-parametric test
# Kruskal-Wallis (non-parametric) is more stable with low number of samples

# number of samples per group
group_counts <- table(richness_different$Art)

# only keep groups with >2 samples
valid_groups <- names(group_counts[group_counts >= 2])

# filter dataset
richness_valid <- richness_different %>%
  filter(Art %in% valid_groups)

# Kruskal-Wallis test
kruskal_species <- kruskal.test(Richness ~ Art, data = richness_valid)

# Display result
kruskal_species


# To see, which groups differ, do a Wilcoxon test
pairwise_results <- pairwise.wilcox.test(richness_valid$Richness,
                                         richness_valid$Art,
                                         p.adjust.method = "BH")  # Benjamini-Hochberg correction reduces the risk of including false positives
pairwise_results



# ---- PERMANOVA , pairwise PERMANOVA & Heatmap ----

library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

comm_mat <- df_num %>%
  rownames_to_column("Bacteria") %>%
  pivot_longer(-Bacteria, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = Sample, values_from = Count, values_fill = 0) %>%
  column_to_rownames("Bacteria") %>%
  t()

group <- factor(sample_to_species[colnames(df_num)])

valid_samples <- names(group)[group %in% names(table(group)[table(group) >= 2])]
group <- droplevels(group[valid_samples])
comm_mat <- comm_mat[valid_samples, ]

dist_mat <- vegdist(comm_mat, method = "jaccard")

adonis2(dist_mat ~ group, permutations = 999)

betadisp_res <- betadisper(dist_mat, group)
anova(betadisp_res)

group_pairs <- combn(levels(group), 2, simplify = FALSE)
pairwise_permanova <- list()

for (pair in group_pairs) {
  subset_ids <- which(group %in% pair)
  dist_sub <- as.dist(as.matrix(dist_mat)[subset_ids, subset_ids])
  group_sub <- droplevels(group[subset_ids])
  adonis_res <- adonis2(dist_sub ~ group_sub, permutations = 999)
  pairwise_permanova[[paste(pair, collapse = "_vs_")]] <- adonis_res$`Pr(>F)`[1]
}

pairwise_permanova_df <- data.frame(
  Comparison = names(pairwise_permanova),
  p_value = unlist(pairwise_permanova)
)

# ----NMDS plot ----

set.seed(42)
nmds <- metaMDS(comm_mat, distance = "jaccard", k = 2, trymax = 100)
nmds$stress

# extract sample scores
scores_df <- as.data.frame(scores(nmds, display = "sites"))
scores_df$Sample <- rownames(scores_df)
scores_df$Nematode <- factor(sample_to_species[scores_df$Sample],
                             levels = c("A. emmatus", "A. mariannae", "Acrobeles sp.",
                                        "Acrobeloides", "Cervidellus", "Mesorhabditis",
                                        "Nothacrobeles", "Panagrolaimus", "PC",
                                        "Stegelleta", "Stegelletina"))


# group sizes for ellipses
group_sizes <- table(scores_df$Nematode)
large_groups <- names(group_sizes[group_sizes >= 4])
small_groups <- names(group_sizes[group_sizes < 4])

# 
desired_levels <- c("PC", "Cervidellus", "Acrobeles sp.", "Mesorhabditis", 
                    "Acrobeloides", "A. emmatus", "A. mariannae", 
                    "Stegelleta", "Stegelletina", "Nothacrobeles")

# 
scores_df$Nematode <- factor(scores_df$Nematode, levels = desired_levels)

# 
nem_colors <- c(
  "PC" = "#a1d99b",
  "Cervidellus"   = "#E69F00",
  "Acrobeloides"  = "#F0E442",
  "Stegelletina"  = "#999999",
  "Acrobeles sp." = "#56B4E9",
  "A. emmatus"    = "#0072B2",
  "A. mariannae"  = "#D55E00",
  "Mesorhabditis" = "#1b9e77",
  "Stegelleta"    = "#CC79A7"
)

library(ggforce)

# NMDS plot
ggplot(scores_df, aes(x = NMDS1, y = NMDS2, color = Nematode, fill = Nematode)) +
  geom_point(size = 5) +

  stat_ellipse(
    data = subset(scores_df, Nematode %in% large_groups),
    aes(group = Nematode),
    type = "t",  
    geom = "polygon",
    alpha = 0.2,
    linewidth = 1
  ) +
  
  geom_mark_ellipse(
    data = subset(scores_df, Nematode %in% small_groups),
    aes(group = Nematode),
    alpha = 0.2,
    expand = 0.05,
    linewidth = 1,
    linetype = "dashed"
  ) +
  
  scale_color_manual(values = nem_colors) +
  scale_fill_manual(values = nem_colors) +
  
  theme_minimal(base_size = 25) +
  theme(
    axis.title = element_text(size = 35),
    axis.text  = element_text(size = 25),
    legend.position = "none",  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black", size = 0.8)
  ) +
  
  labs(x = "NMDS1", y = "NMDS2", title = NULL)




# ---- Cluster dendrogram ---- 


library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggdendro)

comm_mat_original <- df_num %>%
  rownames_to_column("Sample") %>%
  column_to_rownames("Sample") %>%
  t()

# exclude positive control
valid_samples <- setdiff(rownames(comm_mat_original), names(sample_to_species)[sample_to_species == "PC"])
comm_mat_filtered <- comm_mat_original[valid_samples, , drop = FALSE]
group_filtered <- sample_to_species[valid_samples]

# include Cervidellus
if(!"Cervidellus" %in% group_filtered) {
  # Angenommen das Sample heißt "Cervidellus1"
  comm_mat_filtered <- rbind(comm_mat_filtered, comm_mat_original["Cervidellus1", , drop = FALSE])
  group_filtered <- c(as.character(group_filtered), "Cervidellus")
}

group_filtered <- factor(group_filtered)

# aggregate
comm_mat_grouped <- comm_mat_filtered %>%
  as.data.frame() %>%
  mutate(Nematode = group_filtered) %>%
  group_by(Nematode) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
  column_to_rownames("Nematode")

# distance matrix based on aggregated data
dist_mat_grouped <- vegdist(comm_mat_grouped, method = "jaccard")

# hierarchical clustering
hc_grouped <- hclust(dist_mat_grouped, method = "average")

dend_data_grouped <- ggdendro::dendro_data(hc_grouped, type = "rectangle")

y_min <- min(dend_data_grouped$segments$y)

# plot Dendrogramm
ggplot() +
  geom_segment(data = dend_data_grouped$segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               size = 1.2) +
  geom_text(data = dend_data_grouped$labels,
            aes(x = x, y = y_min - 0.9, label = label),
            angle = 90, hjust = 1, size = 5) +
  labs(title = NULL,
       x = NULL,
       y = "Jaccard distance") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 80, l = 10),
    axis.line.y = element_line(color = "black", size = 0.8),
    axis.line.x = element_blank(),
    axis.ticks.y = element_line(color = "black", size = 0.8),   # kleine Striche
    axis.ticks.length = unit(0.25, "cm")                        # Länge der Striche
  ) +
  coord_cartesian(clip = "off", ylim = c(0, 1))  # y-Achse nur von 0 bis 1




# ---- barplot: ratio PC species in samples ----

library(tidyverse)
library(dplyr)

df <- read.csv("species_table_clean_final_binary.csv", header = TRUE)

# only numeric columns
num_cols <- sapply(df, is.numeric)
num_cols["total"] <- FALSE
df_num <- df[, num_cols]

# Only bacteria that are presentin at least one PC
pc_cols <- grep("PC", colnames(df_num), value = TRUE)
pc_rows <- rowSums(df_num[, pc_cols, drop=FALSE]) > 0
df_num <- df_num[pc_rows, ]

# delete columns that are only 0, if there are any
df_num <- df_num[, colSums(df_num) != 0]

sample_names <- colnames(df_num)
sample_to_species <- sapply(sample_names, function(x) {
  if(grepl("Acrobeles_emmatus", x)) return("A. emmatus")
  if(grepl("Acrobeles_mariannae", x)) return("A. mariannae")
  if(grepl("Acrobeles_sp", x)) return("Acrobeles sp.")
  if(grepl("Acrobeloides", x)) return("Acrobeloides")
  if(grepl("Cervidellus", x)) return("Cervidellus")
  if(grepl("Mesorhabditis", x)) return("Mesorhabditis")
  if(grepl("Nothacrobeles", x)) return("Nothacrobeles")
  if(grepl("Panagrolaimus", x)) return("Panagrolaimus")
  if(grepl("PC", x)) return("PC")
  if(grepl("Stegelleta", x)) return("Stegelleta")
  if(grepl("Stegelletina", x)) return("Stegelletina")
  return(NA)
})

# long-format
species_long <- df_num %>%
  rownames_to_column(var="Bacteria") %>%
  pivot_longer(-Bacteria, names_to="Sample", values_to="Count") %>%
  mutate(Taxon = sample_to_species[Sample])

# presence matrix per taxon
taxon_presence <- species_long %>%
  group_by(Taxon, Bacteria) %>%
  summarise(Present = any(Count>0), .groups="drop") %>%
  filter(!is.na(Taxon))

# List of PC-bacteria
pc_bac <- taxon_presence %>% 
  filter(Taxon=="PC" & Present==TRUE) %>% 
  pull(Bacteria)

# bacteria per nematode taxon
taxon_bac <- taxon_presence %>% 
  filter(Taxon != "PC" & Present==TRUE) %>%
  group_by(Taxon) %>%
  summarise(Bacteria = list(Bacteria), .groups="drop")

# number f samples per taxon
samples_per_taxon <- species_long %>%
  filter(Taxon != "PC") %>%
  group_by(Taxon) %>%
  summarise(n_samples = n_distinct(Sample), .groups="drop")

# calculate number of PCnbacteria
pc_total <- length(pc_bac)

pc_overlap <- taxon_bac %>%
  mutate(
    Shared = sapply(Bacteria, function(b) sum(b %in% pc_bac))
  ) %>%
  left_join(samples_per_taxon, by="Taxon") %>%
  mutate(
    Percent_of_PC = (Shared / pc_total) * 100,
    Ratio_per_sample = Percent_of_PC / n_samples
  ) %>%
  select(Taxon, Shared, Percent_of_PC, n_samples, Ratio_per_sample)


pc_overlap$Taxon <- factor(pc_overlap$Taxon, levels = c(
  "PC", "Cervidellus", "Acrobeles sp.", "Mesorhabditis", "Acrobeloides",  
  "A. emmatus", "A. mariannae", "Stegelleta", "Stegelletina"
))

ggplot(pc_overlap, aes(x=Taxon, y=Ratio_per_sample, fill=Taxon)) +
  geom_col(alpha=0.8, color="black") +  # <-- schwarzer Rand
  geom_text(aes(label=round(Ratio_per_sample,2)), vjust=-0.5, size=5) +
  labs(x=NULL, y="PC bacteria per sample [%]") +
  theme_minimal(base_size=15) +
  theme(
    legend.position="none",
    axis.text.x = element_text(angle=45, hjust=1),
    axis.text.y = element_text(size=15),
    axis.title.y = element_text(size=15),
    panel.grid = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.8),
    axis.line.x = element_blank(),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.25, "cm")
  ) +
  scale_fill_manual(values = c(
    "Cervidellus"   = "grey80",
    "Acrobeles sp." = "grey80",
    "Mesorhabditis" = "grey80",
    "Acrobeloides"  = "grey80",
    "A. emmatus"    = "grey80",
    "A. mariannae"  = "grey80",
    "Stegelleta"    = "grey80",
    "Stegelletina"  = "grey80"
  ))
