suppressMessages(library(GenomicRanges))
suppressMessages(library(AnnotationHub))
ah <- AnnotationHub()
query_data <- subset(ah, preparerclass == "excluderanges")

args <- commandArgs(trailingOnly = TRUE)
query <- args[1]

exclude <- query_data[[query]] %>% 
  sort() %>% 
  keepStandardChromosomes(pruning.mode = "tidy")

outfile <- args[2]
write.table(as.data.frame(exclude),
            file = outfile,
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)