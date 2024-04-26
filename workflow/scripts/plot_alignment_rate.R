# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

# Get Bowtie2 alignment rates (in log files)
files <- snakemake@input

# Create df for storing alignment rates
df <- as.data.frame(matrix(ncol = 2, nrow = 0))
names(df) <- c("sample", "alignment.rate")

# Get sample name and mapping rates from log files
for (i in seq(files)){
  sample <- system(paste0("echo ", files[i], "| sed 's/.log//'"), intern = TRUE)
  sample <- basename(sample)

  rate <- system(paste0('grep "overall alignment rate" ',
                        files[i],
                        " | awk '{print $1}' | sed 's/%$//'"),
                        intern = TRUE)
  rate <- as.numeric(rate)

  # add to df
  df[i, "sample"] <- sample
  df[i, "mapping.rate"] <- rate
}

# Remove prepended X from samples names
# (only happens when they start with a number)
df$sample <- str_remove(df$sample, "^X")

# Create plot
p <- ggplot(df, aes(x = sample, y = mapping.rate)) +
  geom_bar(stat = "identity",
           fill = "aquamarine4",
           colour = "black") +
  theme_cowplot(16) +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylab("Overall alignment rate (%)") +
  xlab(NULL)

# Save plot
ggsave(snakemake@output[[1]],
       p,
       height = 6,
       width = length(df$sample) * 1.2)

# close log file
sink(log, type = "output")
sink(log, type = "message")
