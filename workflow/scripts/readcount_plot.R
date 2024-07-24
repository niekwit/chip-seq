# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)
library(reshape2)

# Load read counts to df
df <- read.csv(snakemake@input[[1]], header = TRUE)
df <- df[order(df$sample), ]

# Create df for plotting
df_melt <- melt(df, value.name = "counts")

# Plot data
p <- ggplot(df_melt, aes(sample, counts)) +
  geom_bar(aes(fill = variable),
           position = "dodge",
           stat = "identity",
           color = "black") +
  theme_cowplot(16) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        legend.position = c(0.02, 0.9)) +
  scale_y_continuous(limits = c(0, max(df_melt$counts) * 1.15),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylab("Read count") +
  xlab(NULL) +
  guides(fill = guide_legend(NULL),
         y = guide_axis(minor.ticks = TRUE)) +
  scale_fill_manual(labels = c("Pre-deduplication",
                               "Post-deduplication"),
                    values = c("pre.dedup_counts" = "#419179",
                               "post.dedup_counts" = "#56A8CBFF"))

# Save plot
ggsave(snakemake@output[[1]],
       p,
       height = 6,
       width = length(df$sample) * 1.4)

# close log file
sink(log, type = "output")
sink(log, type = "message")