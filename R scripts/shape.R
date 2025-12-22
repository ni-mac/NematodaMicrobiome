# load .csv file
df <- read.csv("shape_table_clean_final_binary.csv", header = TRUE)

# only select numeric columns (without "total", "cell.width" and "cell.length") and "shape"
num_cols <- sapply(df, is.numeric)
num_cols[c("total", "cell.width", "cell.length")] <- FALSE
num_cols["shape"] <- TRUE

# delete all other columns
df_num <- df[, num_cols]

# create barplot to identify underrepresented bacterial shapes
ggplot(df, aes(x = shape)) +
  geom_bar(fill = "steelblue") +
  labs(x = "Bacterial Shape", y = "Number of species") +
  theme_minimal()

# exclude underrepresented shapes 
df_num <- df_num[!df_num$shape %in% c("pleomorphic", "filament", "helical"), ]

ggplot(df_num, aes(x = shape)) +
  geom_bar(fill = "steelblue") +
  labs(x = "Bacterial Shape", y = "Number of species") +
  theme_minimal()

# Remove columns with all zeros
zero_cols <- sapply(df_num, function(x) all(x == 0, na.rm = TRUE))
df_num <- df_num[, !zero_cols]

# Assign samples to nematode taxa
sample_names <- colnames(df_num)
sample_to_shape <- sapply(sample_names, function(x) {
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

library(dplyr)
library(tidyr)
library(tibble)  # for rownames_to_column

# bring dataset into long format
shape_long <- df_num %>%
  rownames_to_column(var = "Bacteria") %>%
  pivot_longer(
    cols = -c(Bacteria, shape),
    names_to = "Sample",
    values_to = "Count"
  ) %>%
  mutate(Art = sample_to_shape[Sample])

# Calculate composition of shapes
shape_counts <- shape_long %>%
  group_by(Sample, Art, shape) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")

# Ensure "PC" is first on x-axis
shape_counts$Art <- factor(
  shape_counts$Art,
  levels = c("PC", setdiff(unique(shape_counts$Art), "PC"))
)

# Ensure specific order of Art
shape_counts$Art <- factor(
  shape_counts$Art, 
  levels = c("PC", "Acrobeles sp.", "Stegelletina", "Cervidellus", "Acrobeloides", 
             "A. emmatus", "Mesorhabditis", "A. mariannae", 
              "Stegelleta")
)

library(ggplot2)

# Boxplot
ggplot(shape_counts, aes(x = Art, y = Count, fill = shape)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.5, position = position_dodge(width = 0.8)) +
  labs(
    x = NULL, 
    y = "Number of bacterial shapes per nematode taxa",
    fill = NULL
  ) +
  scale_y_log10() +
  scale_fill_manual(values = c(
    "coccus" = "#1f78b4",
    "rod"    = "#33a02c",
    "ovoid"  = "yellow"
  )) +
  theme_minimal(base_size = 22) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.25, "cm"),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        linetype = 0,
        shape = 22,
        color = "black"
      ),
      keywidth = 1.5,
      keyheight = 1.5
    )
  )



# stacked barplot
shape_sums <- shape_long %>%
  group_by(Art, shape) %>%
  summarise(TotalCount = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  group_by(Art) %>%
  mutate(Proportion = TotalCount / sum(TotalCount))

shape_sums$Art <- factor(shape_sums$Art, 
                         levels = c("PC", "Acrobeles sp.", "Stegelletina", "Cervidellus", "Acrobeloides", 
                                     "A. emmatus", "Mesorhabditis", "A. mariannae", "Stegelleta"))

ggplot(shape_sums, aes(x = Art, y = Proportion, fill = shape)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c(
    "coccus" = "#1f78b4",
    "rod"    = "#33a02c",
    "ovoid"  = "yellow"
  )) +
  labs(x = NULL, y = "Proportion of detected bacterial shapes") +
  theme_minimal(base_size = 22) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 1),
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 1),
    axis.line.x = element_blank(),
    axis.ticks.length = unit(0.25, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )




# ---- tests for all shapes simultaneously ----

library(dplyr)
library(tidyr)

run_tests_per_shape <- function(data, shape_name) {
  
  shape_data <- data %>% filter(shape == shape_name)
  
  shapiro_res <- shape_data %>%
    group_by(Art) %>%
    summarise(
      W = if(n() >= 3 & length(unique(Count)) > 1) shapiro.test(Count)$statistic else NA_real_,
      p_value = if(n() >= 3 & length(unique(Count)) > 1) shapiro.test(Count)$p.value else NA_real_,
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(shape = shape_name)
  
  valid_groups <- names(table(shape_data$Art)[table(shape_data$Art) >= 2])
  shape_valid <- shape_data %>% filter(Art %in% valid_groups)
  
  kruskal_res <- kruskal.test(Count ~ Art, data = shape_valid)
  kruskal_df <- kruskal_res$parameter  # Freiheitsgrade
  kruskal_p <- kruskal_res$p.value
  
  kruskal_out <- data.frame(
    shape = shape_name,
    df = kruskal_df,
    kruskal_p = kruskal_p
  )
  
  pairwise_res <- pairwise.wilcox.test(
    shape_valid$Count,
    shape_valid$Art,
    p.adjust.method = "BH"
  )
  
  pairwise_long <- as.data.frame(as.table(pairwise_res$p.value)) %>%
    rename(nematode1 = Var1, nematode2 = Var2, p_value = Freq) %>%
    mutate(
      signif = ifelse(!is.na(p_value) & p_value < 0.05, "*", ""),
      shape = shape_name
    )
  
  list(
    shapiro = shapiro_res,
    kruskal = kruskal_out,
    pairwise = pairwise_long
  )
}

shapes <- unique(shape_counts$shape)

all_results <- lapply(shapes, function(s) run_tests_per_shape(shape_counts, s))

# Shapiro
shapiro_all <- bind_rows(lapply(all_results, `[[`, "shapiro"))
write.csv(shapiro_all, "shapiro_results_all_shapes.csv", row.names = FALSE)

# Kruskal-Wallis
kruskal_all <- bind_rows(lapply(all_results, `[[`, "kruskal"))
write.csv(kruskal_all, "kruskal_results_all_shapes.csv", row.names = FALSE)

# Pairwise Wilcoxon
pairwise_all <- bind_rows(lapply(all_results, `[[`, "pairwise"))
write.csv(pairwise_all, "pairwise_wilcox_all_shapes.csv", row.names = FALSE)



# ---- Differences in Composition ---- 

library(vegan)

comm_mat <- shape_counts %>%
  pivot_wider(
    names_from = shape,
    values_from = Count,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample") %>%
  mutate(across(c("coccus", "ovoid", "rod"), as.numeric))

comm_mat_numeric <- comm_mat %>% select(-Art)

# Jaccard distance for community dissimilarity
dist_mat <- vegdist(comm_mat_numeric, method = "jaccard")

# PERMANOVA: tests for overall group differences in shape composition
group <- factor(comm_mat$Art)
adonis2(dist_mat ~ group, permutations = 999)

# betadisper test to see if dispersions are homogeneous
betadisp_res <- betadisper(dist_mat, group)
anova(betadisp_res)
# p>0.05 = dispersions within groups are statistically similar, whihc is good for permanova


library(vegan)
library(dplyr)

# Compute Jaccard distance
dist_mat <- vegdist(comm_mat_numeric, method = "jaccard")

# --- Pairwise PERMANOVA ---
group_levels <- levels(group)
group_pairs <- combn(group_levels, 2, simplify = FALSE)

pairwise_permanova <- list()

for (pair in group_pairs) {
  # Subset the distance matrix and group factor for these two groups
  subset_ids <- which(group %in% pair)
  dist_sub <- as.dist(as.matrix(dist_mat)[subset_ids, subset_ids])
  group_sub <- droplevels(group[subset_ids])
  
  # Run PERMANOVA on the subset
  adonis_res <- adonis2(dist_sub ~ group_sub, permutations = 999)
  
  # Store the p-value for this comparison
  pair_name <- paste(pair, collapse = "_vs_")
  pairwise_permanova[[pair_name]] <- adonis_res$`Pr(>F)`[1]
}

# Convert results to data.frame
pairwise_permanova_df <- data.frame(
  Comparison = names(pairwise_permanova),
  p_value = unlist(pairwise_permanova)
)

# Adjust for multiple testing (Benjamini-Hochberg)
pairwise_permanova_df$p_adj <- p.adjust(pairwise_permanova_df$p_value, method = "BH")

# View results
pairwise_permanova_df



# ---- NMDS: Similarity of shape composition----
library(ggforce)

nmds <- metaMDS(dist_mat, k = 2, trymax = 100)
nmds$stress 

# Extract coordinates
nmds_coords <- as.data.frame(scores(nmds))
nmds_coords$Sample <- rownames(nmds_coords)
nmds_coords$Nematode <- comm_mat$Art 

# Group size classification
group_sizes <- table(nmds_coords$Nematode)
large_groups <- names(group_sizes[group_sizes >= 3])
small_groups <- names(group_sizes[group_sizes < 3])

# Define custom colors
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

ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = Nematode, fill = Nematode)) +
  geom_point(size = 5) +
  stat_ellipse(
    data = subset(nmds_coords, Nematode %in% large_groups),
    aes(group = Nematode),
    type = "t",
    geom = "polygon",
    alpha = 0.2,
    linewidth = 1
  ) +
  geom_mark_ellipse(
    data = subset(nmds_coords, Nematode == "Acrobeloides"),
    alpha = 0.2,
    expand = 0.025,
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_mark_ellipse(
    data = subset(nmds_coords, Nematode %in% small_groups),
    aes(group = Nematode),
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
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 20, color = "black"),  
    axis.title.y = element_text(size = 20),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.line = element_line(color = "black", linewidth = 0.8)
  )






# ---- cluster dendrogramm: Similarity of shape composition

library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggdendro)

# Aggregate counts of each bacterial shape per nematode species
comm_mat_grouped <- shape_counts %>%
  select(Sample, Art, shape, Count) %>%
  pivot_wider(names_from = shape, values_from = Count, values_fill = 0) %>%
  group_by(Art) %>%
  summarise(across(c(coccus, rod, ovoid), \(x) mean(x, na.rm = TRUE))) %>%
  column_to_rownames("Art")


# Compute Jaccard distance for community dissimilarity
dist_mat_grouped <- vegdist(comm_mat_grouped, method = "jaccard")

# Perform hierarchical clustering using average linkage
hc_grouped <- hclust(dist_mat_grouped, method = "average")

# Prepare dendrogram data for ggplot
dend_data_grouped <- ggdendro::dendro_data(hc_grouped, type = "rectangle")

y_min <- min(dend_data_grouped$segments$y)

ggplot() +
  geom_segment(data = dend_data_grouped$segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               size = 1.2) +
  geom_text(data = dend_data_grouped$labels,
            aes(x = x, y = y_min - 0.9, label = label),
            angle = 90, hjust = 1, size = 5) +
  labs(x = NULL, y = "Jaccard distance") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 80, l = 10),
    axis.line.y = element_line(color = "black", size = 0.8),  # y-axis line
    axis.line.x = element_blank()
  ) +
  coord_cartesian(clip = "off", ylim = c(0, 1))


# ---- Ternary plot ----

# another ggplot2 version is needed which is compatible with ggtern. 
# Restart the R session, remove packages as seen below, then load code until shape_sums (line 112)
# then install and load packages below

remove.packages("ggplot2")
remove.packages("ggtern")

# old ggplot2
install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz",
                 repos = NULL, type = "source")

# ggtern
install.packages("https://cran.r-project.org/src/contrib/Archive/ggtern/ggtern_3.4.1.tar.gz",
                 repos = NULL, type = "source")

library(ggplot2)
library(ggtern)


shape_tern <- shape_sums %>%
  select(Art, shape, Proportion) %>%
  pivot_wider(
    names_from = shape,
    values_from = Proportion,
    values_fill = 0
  ) %>%
  rowwise() %>%
  mutate(total = coccus + rod + ovoid) %>%
  mutate(
    coccus = coccus / total,
    rod    = rod / total,
    ovoid  = ovoid / total
  ) %>%
  select(-total)

nem_colors <- c(
  "PC" = "#a1d99b",
  "Cervidellus"   = "#E69F00",
  "Acrobeloides"  = "#F0E442",
  "Stegelletina"  = "#999999",
  "Acrobeles sp." = "#56B4E9",
  "A. emmatus"    = "#0072B2",
  "A. mariannae"  = "#D55E00",
  "Mesorhabditis" = "#1b9e77",
  "Nothacrobeles" = "#666666",
  "Panagrolaimus" = "#1f78b4",
  "Stegelleta"    = "#CC79A7"
)

ggtern(data = shape_tern, aes(x = coccus, y = rod, z = ovoid, color = Art)) +
  geom_point(size = 10) +
  scale_color_manual(values = nem_colors) +
  theme_bw(base_size = 20) +
  labs(
    T = "Rod [%]",
    L = "Coccus [%]",
    R = "Ovoid [%]"
  ) +
  theme(
    tern.axis.arrow.show = FALSE,
    plot.title = element_blank(),
    legend.title = element_blank(), 
    legend.text = element_text(size = 20),
    legend.position = c(1, 1),
    legend.justification = c("right", "top")
  )

