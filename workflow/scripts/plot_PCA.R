# redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# load required libraries
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(scales)

#### PCA plot ####
# Load PCA data
data <- read.delim(snakemake@input[[1]], header = TRUE, skip = 1)

# Load samples.csv
csv <- read.csv("config/samples.csv")
samples <- csv$sample
controls <- csv$control


product <- length(unique(csv$genotype)) * length(unique(csv$treatment))
colours <- brewer.pal(8, "Dark2")
shapes <- c(16, 15) # Square for IP, Circle for Control


# Keep only components 1 and 2, transpose and add sample information
df <- data[1:2, ] %>%
  dplyr::select(-c("Component", "Eigenvalue")) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sample = rownames(.)) %>%
  rename(PC1 = 1, PC2 = 2) %>%
  mutate(
    shape = ifelse(
      str_detect(sample, paste(controls, collapse = "|")),
      16,
      15
    )
  )
df$colour <- NA
for (i in seq_along(unique(samples))) {
  colour <- colours[i]
  sample <- unique(samples)[i]
  df$colour <- ifelse(str_detect(df$sample, sample), colour, df$colour)
}
df$shape <- factor(df$shape, levels = c(15, 16)) # Ensure correct order of shapes

# Calculate variance explained for each PC
PC1_var <- round((data$Eigenvalue[1] / sum(data$Eigenvalue)) * 100, 1)
PC2_var <- round((data$Eigenvalue[2] / sum(data$Eigenvalue)) * 100, 1)

# Create PCA plot
p <- ggplot(
  df,
  mapping = aes(x = PC1, y = PC2, colour = colour, shape = shape)
) +
  geom_point(size = 8) +
  geom_label_repel(
    data = df,
    aes(label = sample, fill = NULL),
    size = 5,
    max.overlaps = 30
  ) +
  scale_shape_manual(values = c(15, 16)) +
  theme_cowplot(18) +
  labs(
    x = paste0("PC1: ", PC1_var, "% variance"),
    y = paste0("PC2: ", PC2_var, "% variance")
  ) +
  theme(legend.position = "none")

# Save plot
ggsave(snakemake@output[["pca"]], p, height = 4, width = 6)

#### Scree plot ####
# Scale factor for utilising whole second y-axis range
# https://stackoverflow.com/questions/65559901/add-a-second-y-axis-to-ggplot
scalefactor <- max(data$Eigenvalue) / 100

# Prepare data for scree plot
df <- data %>%
  dplyr::select(c("Component", "Eigenvalue")) %>%
  mutate(Component = paste0("PC", Component)) %>%
  mutate(
    cumulative_variance = (cumsum(Eigenvalue) /
      sum(Eigenvalue) *
      100 *
      scalefactor)
  )

# Create scree plot
s <- ggplot(df, aes(Component, cumulative_variance)) +
  geom_bar(
    aes(Component, Eigenvalue),
    stat = "identity",
    colour = "black",
    fill = "#419179"
  ) +
  geom_line(
    mapping = aes(x = Component, y = cumulative_variance, group = 1),
    colour = "red",
    linewidth = 1
  ) +
  geom_point(
    mapping = aes(x = Component, y = cumulative_variance),
    colour = "red",
    fill = "white",
    shape = 21,
    size = 5,
    stroke = 1.5
  ) +
  theme_cowplot(15) +
  theme(
    axis.title.y.right = element_text(color = "red"),
    axis.text.y.right = element_text(color = "red")
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(
      transform = ~ .x / scalefactor,
      breaks = seq(0, 100, 25),
      name = "Cumulative variance explained (%)"
    ),
    expand = expansion(mult = c(0, .05))
  ) +
  labs(x = "Principal component", y = "Eigenvalue")

# Save plot
ggsave(snakemake@output[["scree"]], s, height = 4, width = 6)
