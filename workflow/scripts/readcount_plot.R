# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type="output")
sink(log, type="message")

library(tidyverse)
library(cowplot)











