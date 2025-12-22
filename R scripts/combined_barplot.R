library(dplyr)
library(tidyr)
library(ggplot2)

# bacterial species per nematode species

df <- data.frame(
  Nematode = c("Cervidellus", "Acrobeles sp.", "Mesorhabditis",
               "Acrobeloides", "A. emmatus", "A. mariannae",
               "Stegelleta", "Stegelletina", "PC"),
  `species.level` = c(27, 47, 86, 56, 108, 129, 118, 68, 2081),
  `shape.level`   = c(20, 31, 49, 37, 73, 86, 73, 49, 1049),
  `size.level`    = c(7, 21, 27, 13, 37, 44, 37, 24, 741)
)

# Long-Format
df_long <- df %>%
  pivot_longer(
    cols = c(`species.level`, `shape.level`, `size.level`),
    names_to = "Level",
    values_to = "Count"
  )

df_long$Level <- factor(df_long$Level, levels = c("size.level", "shape.level", "species.level"))

df_long$Nematode <- factor(df_long$Nematode, 
                           levels = rev(c("PC", setdiff(df_long$Nematode, "PC"))))

ggplot(df_long, aes(x = Nematode, y = Count, fill = Level)) +
  geom_bar(
    stat = "identity",
    position = position_dodge2(width = 0.8, padding = 0.25),
    color = "black"
  ) +
  scale_y_log10() +
  scale_fill_manual(
    values = c(
      "species.level" = "black",
      "shape.level"   = "grey40",
      "size.level"    = "grey80"
    ),
    labels = c(
      "species.level" = "Species (n=1373)",
      "shape.level"   = "Shape (n=741)",
      "size.level"    = "Size (n=460)"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of detected bacterial species",
    fill = "Level of analysis"
  ) +
  guides(
    fill = guide_legend(
      nrow = 3,
      title.position = "top",
      title.hjust = 0.5,
      reverse = TRUE
    )
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.y = element_text(hjust = 1),
    axis.ticks.y = element_line(color = "black"),    
    axis.line.y = element_line(color = "black"),     
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    legend.position = "right"
  ) +
  coord_flip()



# nematode individuals per species

data_new <- data.frame(
  Nema = c("PC","Cervidellus", "Acrobeles sp.", "Mesorhabditis",
           "Acrobeloides", "A. emmatus", "A. mariannae",
           "Stegelleta", "Stegelletina"),
  AfterExtraction = c(0, 10, 15, 25, 20, 30, 30, 30, 41),
  Species = c(0, 10, 15, 25, 15, 30, 30, 30, 41),
  Shape   = c(0, 10, 10, 25, 15, 30, 30, 30, 41),
  Size    = c(0, 10, 10, 20, 15, 30, 30, 30, 41)
)

data_long <- data_new %>%
  pivot_longer(
    cols = c(AfterExtraction, Species, Shape, Size),
    names_to = "Level",
    values_to = "Count"
  )

data_long$Level <- factor(data_long$Level, levels = c("AfterExtraction", "Species", "Shape", "Size"))

nematode_order <- rev(c("PC","Cervidellus", "Acrobeles sp.", "Mesorhabditis",
                        "Acrobeloides", "A. emmatus", "A. mariannae",
                        "Stegelleta", "Stegelletina"))
data_long$Nema <- factor(data_long$Nema, levels = nematode_order)

ggplot(data_long, aes(x = Nema, y = Count, fill = Level)) +
  geom_bar(
    stat = "identity",
    position = position_dodge2(width = 0.8, padding = 0.25, reverse = TRUE), # Balken umkehren
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "AfterExtraction" = "grey100",
      "Species" = "black",
      "Shape" = "grey40",
      "Size" = "grey80"
    ),
    labels = c(
      "AfterExtraction" = "After Extraction",
      "Species" = "Species",
      "Shape" = "Shape",
      "Size" = "Size"
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 4,         
      title.position = "top"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of nematode individuals per species",
    fill = "Level of analysis"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line.x = element_line(color = "black"),   
    axis.line.y = element_line(color = "black"),   
    axis.ticks.x = element_line(color = "black"),  
    axis.ticks.y = element_line(color = "black"),  
    legend.position = "right"
  ) +
  coord_flip()

