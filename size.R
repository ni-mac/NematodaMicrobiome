library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)       
library(ggforce)     
library(ggdendro)    
library(tibble) 

df <- read.csv("measurements_table_clean_final_binary.csv", header = TRUE)

# only select numeric columns (without "total", "cell.width" and "cell.length") and "shape"
num_cols <- sapply(df, is.numeric)
num_cols[c("total", "ratio", "cell.width", "cell.length")] <- FALSE
num_cols["shape"] <- TRUE

# delete column "ratio"
df_num <- df[, num_cols]

zero_cols <- sapply(df_num, function(x)
  is.numeric(x) && all(x == 0))
names(df_num)[zero_cols]
df_num <- df_num[, !zero_cols]  

df_num_measure <- subset(df_num, select = -shape)


# Assign samples to nematode taxa
sample_names <- colnames(df_num_measure)[colnames(df_num_measure) != "size.class"]
sample_to_measurement <- sapply(sample_names, function(x) {
  if(grepl("Acrobeles_emmatus", x)) return("A. emmatus")
  if(grepl("Acrobeles_mariannae", x)) return("A. mariannae")
  if(grepl("Acrobeles_sp", x)) return("Acrobeles sp.")
  if(grepl("Acrobeloides", x)) return("Acrobeloides")
  if(grepl("Cervidellus", x)) return("Cervidellus")
  if(grepl("Mesorhabditis", x)) return("Mesorhabditis")
  if(grepl("PC", x)) return("PC")
  if(grepl("Stegelleta", x)) return("Stegelleta")
  if(grepl("Stegelletina", x)) return("Stegelletina")
  return(NA)
})

# bring dataset into long format
measurement_long <- df_num_measure %>%
  pivot_longer(
    cols = all_of(sample_names),
    names_to = "Sample",
    values_to = "presence.absence"
  ) %>%
  mutate(nematode = sample_to_measurement[Sample])


# ---- stacked barplot ----

size_summary <- measurement_long %>%
  group_by(nematode, size.class) %>%
  summarise(count = sum(presence.absence), .groups = "drop")

size_relative <- size_summary %>%
  group_by(nematode) %>%
  mutate(freq = count / sum(count))

size_labels <- c(
  "1" = "1.0 – <1.5",
  "2" = "1.5 – <2.0",
  "3" = "2.0 – <2.5",
  "4" = "2.5 – <3.0",
  "5" = "3.0 – <3.5",
  "6" = "3.5 – <4.0",
  "7" = "4.0 – <5.0",
  "8" = "5.0 – <6.0",
  "9" = "6.0 – <10.0",
  "10" = "10 – 61.5"
)


library(ggplot2)

# desired order
desired_order <- c(
  "PC", "Acrobeloides", "Cervidellus", "Acrobeles sp.",
  "Mesorhabditis", "Stegelleta", "Stegelletina",
  "A. emmatus", "A. mariannae"
)


size_relative$nematode <- factor(size_relative$nematode, levels = desired_order)

# Plot
ggplot(size_relative, aes(x = nematode, y = freq, fill = factor(size.class))) +
  geom_col() +
  scale_fill_discrete(name = "length/width ratio", labels = size_labels) +
  labs(
    x = NULL,
    y = "Proportion of bacterial size classes"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.8),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 25),
    axis.text.y = element_text(size = 22),
    axis.title.y = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 20),
    legend.position = "none"  # Legende komplett entfernen
  )


# ---- statistics for number of shapes ----

library(dplyr)

# Step 1: Count bacteria per sample and size.class
bacteria_counts_by_class <- measurement_long %>%
  group_by(Sample, nematode, size.class) %>%
  summarise(count = sum(presence.absence), .groups = "drop")

# Step 2: Shapiro-Wilk test per nematode and size.class
shapiro_results <- bacteria_counts_by_class %>%
  group_by(size.class, nematode) %>%
  summarise(
    n = n(),
    unique_count_values = length(unique(count)),
    shapiro_p = if(n >= 3 & unique_count_values > 1) {
      shapiro.test(count)$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  )


# number of bacteria per sample and size class
bacteria_counts_by_class <- measurement_long %>%
  group_by(Sample, nematode, size.class) %>%
  summarise(count = sum(presence.absence), .groups = "drop")

# Kruskal Wallis test 
kruskal_results <- bacteria_counts_by_class %>%
  group_by(size.class) %>%
  summarise(
    kruskal_p = kruskal.test(count ~ nematode)$p.value
  )

kruskal_results

pairwise_results <- bacteria_counts_by_class %>%
  group_by(size.class) %>%
  summarise(
    wilcox = list(pairwise.wilcox.test(count, nematode, p.adjust.method = "BH")),
    .groups = "drop"
  )

library(tidyr)


pairwise_long <- pairwise_results %>%
  mutate(p_matrix = lapply(wilcox, function(x) {
    as.data.frame(as.table(x$p.value))
  })) %>%
  select(size.class, p_matrix) %>%
  unnest(p_matrix) %>%
  rename(nematode1 = Var1, nematode2 = Var2, p_value = Freq)

# Add significance column
pairwise_long <- pairwise_long %>%
  mutate(signif = ifelse(!is.na(p_value) & p_value < 0.05, "*", ""))

# View result
pairwise_long



# ---- statistics for composition of shapes -----

bacteria_wide <- bacteria_counts_by_class %>%
  select(Sample, size.class, count) %>%
  pivot_wider(names_from = size.class, values_from = count, values_fill = 0)

sample_metadata <- bacteria_counts_by_class %>%
  select(Sample, nematode) %>%
  distinct()

count_matrix <- as.matrix(bacteria_wide[,-1])  # Remove Sample column

# PERMANOVA (Jaccard) 
permanova_result <- adonis2(count_matrix ~ nematode, data = sample_metadata,
                            method = "jaccard", permutations = 999)
print(permanova_result)

# Betadispersion (homogeneity of dispersion)
dist_matrix <- vegdist(count_matrix, method = "jaccard")
beta_disp <- betadisper(dist_matrix, sample_metadata$nematode)
anova(beta_disp)

# Pairwise PERMANOVA 
groups <- as.character(unique(sample_metadata$nematode))

pairwise_grid <- expand.grid(g1 = groups, g2 = groups) %>%
  filter(g1 != g2)

rownames(count_matrix) <- bacteria_wide$Sample

groups <- as.character(unique(sample_metadata$nematode))

pairwise_grid <- expand.grid(g1 = groups, g2 = groups) %>%
  filter(g1 != g2)

# calculation of p-values
pairwise_permanova <- pairwise_grid %>%
  rowwise() %>%
  mutate(
    adonis_p = {
      sel_samples <- sample_metadata %>% filter(nematode %in% c(g1, g2))
      sel_matrix <- count_matrix[sel_samples$Sample, , drop = FALSE]
      if(length(unique(sel_samples$nematode)) == 2) {
        adonis2(sel_matrix ~ nematode, data = sel_samples,
                method = "jaccard", permutations = 999)$`Pr(>F)`[1]
      } else {
        NA_real_
      }
    }
  )

print(pairwise_permanova)




# data matrix for NMDS
count_matrix_nmds <- as.matrix(bacteria_wide[,-1])  
rownames(count_matrix_nmds) <- bacteria_wide$Sample

# Jaccard distance matrix
dist_matrix_nmds <- vegdist(count_matrix_nmds, method = "jaccard")

# NMDS 
set.seed(123)  
nmds <- metaMDS(dist_matrix_nmds, k = 2, trymax = 100)

nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Sample <- rownames(nmds_scores)

nmds_scores <- nmds_scores %>%
  left_join(sample_metadata, by = "Sample")

group_sizes <- table(nmds_scores$nematode)
large_groups <- names(group_sizes[group_sizes >= 3])
small_groups <- names(group_sizes[group_sizes <= 3])

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

# NMDS plot
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = nematode, fill = nematode)) +
  geom_point(size = 5) +
  stat_ellipse(
    data = subset(nmds_scores, nematode %in% large_groups),
    aes(group = nematode),
    type = "t",
    geom = "polygon",
    alpha = 0.2,
    linewidth = 1
  ) +
  geom_mark_ellipse(
    data = subset(nmds_scores, nematode %in% small_groups),
    aes(group = nematode),
    alpha = 0.2,
    expand = 0.025,
    linewidth = 1,
    linetype = "dashed"
  ) +
  scale_color_manual(values = nem_colors) +
  scale_fill_manual(values = nem_colors) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.line = element_line(color = "black", linewidth = 0.8)
  )

nmds$stress


# ---- dendrogram ----
count_matrix_grouped <- bacteria_counts_by_class %>%
  group_by(nematode, size.class) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = size.class, values_from = count, values_fill = 0) %>%
  column_to_rownames("nematode") %>%
  as.matrix()

dist_group <- vegdist(count_matrix_grouped, method = "jaccard")
hc_group <- hclust(dist_group, method = "average")  # UPGMA

dendro_group <- dendro_data(hc_group, type = "rectangle")

ggplot() +
  geom_segment(
    data = dendro_group$segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    size = 1.5,   # Linien dicker
    color = "black"
  ) +
  geom_text(
    data = dendro_group$labels,
    aes(x = x, y = y - 0.02, label = label),
    angle = 90,
    hjust = 1,
    size = 6,
    color = "black"
  ) +
  labs(x = NULL, y = "Jaccard distance") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 18, color = "black"),
    axis.ticks.y = element_line(color = "black", size = 1.2),  # dickere Ticks
    axis.ticks.length = unit(0.4, "cm"),                        # längere Ticks
    axis.line.y  = element_line(color = "black", size = 1.5),   # dickere Y-Achse
    panel.grid = element_blank(),
    legend.position = "none"
  )


