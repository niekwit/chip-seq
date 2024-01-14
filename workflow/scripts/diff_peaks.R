library(DiffBind)
library(dplyr)
library(rtracklayer)
library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggupset)
library(ReactomePA)
library(stringr)
library(glue)

genome <- snakemake@params[1]

#load annotation databases
if (genome == "hg38"){
  library(EnsDb.Hsapiens.v86) ##CHANGE TO ENSEMBL
  edb <- EnsDb.Hsapiens.v86 ##CHANGE TO ENSEMBL
  annoDb <- "org.Hs.eg.db"
  organism <- "human"
} else if (genome == "hg19"){
  library(TxDb.Hsapiens.UCSC.hg19.knownGene) ##CHANGE TO ENSEMBL
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene ##CHANGE TO ENSEMBL
}

#get sample info
csv <- read.csv("samples.csv")

#create output dir
out.dir <- "diffpeaks"
dir.create(out.dir, showWarnings = FALSE)

#create empty sample sheet
sampleSheet <- as.data.frame(matrix(nrow = nrow(csv), ncol = 0))

###prepare sampleSheet contents
#create SampleIDs
SampleIDs <- paste0(csv$genotype,"_",csv$sample,"_",csv$factor,"_",csv$treatment)

#create Conditions
Condition <- csv$genotype

#create Treatment
Treatment <- csv$treatment

#create BAM file names
bamReads <- paste0("mapped/",csv$genotype,"_",csv$sample,"_",csv$factor,"_",csv$treatment,"_dedup.bam") 
bamControls <- paste0("mapped/",csv$genotype,"_",csv$sample,"_input_",csv$treatment,"_dedup.bam") 

#create ControlIDs
ControlID <- paste0(csv$genotype,"_",csv$sample,"_input_",csv$treatment)

#create peak file names
Peaks <- paste0("peaks/", 
                csv$genotype,"_",csv$sample,"_",csv$factor,"_",csv$treatment,"_vs_",
                csv$genotype,"_",csv$sample,"_input_",csv$treatment,"/",
                csv$genotype,"_",csv$sample,"_",csv$factor,"_",csv$treatment,"_peaks.xls") 

#add content to sampleSheet
sampleSheet <- sampleSheet %>%
  mutate(SampleID = SampleIDs) %>%
  mutate(Tissue = NA) %>%
  mutate(Factor = csv$factor) %>%
  mutate(Condition = paste0(csv$genotype,"_",csv$factor,"_",csv$treatment)) %>% #simplifies analysis
  mutate(Treatment = csv$treatment) %>%
  mutate(Replicate = csv$sample) %>%
  mutate(bamReads = bamReads) %>%
  mutate(ControlID = ControlID) %>%
  mutate(bamControl = bamControls) %>%
  mutate(Peaks = Peaks) %>%
  mutate(PeakCaller = "macs") 
  
#write sampleSheet to file
write.csv(sampleSheet,
          file=file.path("diffpeaks","diffbind_samples.csv"),
          row.names = FALSE)


#initiate analysis
#diffbind <- dba.analyze(sampleSheet,
#                        bGreylist=FALSE,
#                        design="~Condition")

diffbind <- dba(sampleSheet=sampleSheet) %>% #skipped blacklist step as these regions are already absent
            dba.count() %>%
            dba.normalize() %>%
            dba.contrast(minMembers = 2) %>%
            dba.analyze(design="~Condition",
                        bGreylist=FALSE)


#save diffbind object to file 
save(diffbind, file = "diffpeaks/diffbind.RData")

#plot correlation heatmap
pdf("diffpeaks/sample_correlation_heatmap.pdf")
plot(diffbind)
dev.off()

#plot PCA
pdf("diffpeaks/PCA.pdf")
dba.plotPCA(diffbind,DBA_CONDITION,label=DBA_ID)
dev.off()

#generate binding affinity heatmaps
hmap <- colorRampPalette(c("red", "black", "forestgreen"))(n = 13)
pdf("diffpeaks/affinity_heatmap.pdf")
dba.plotHeatmap(diffbind, correlations=FALSE,scale="row", colScheme = hmap)
dev.off()

#Plot profiles
tryCatch({
  pdf("diffpeaks/plotProfile.pdf")
  dba.plotProfile(diffbind, samples=diffbind$masks$All, doPlot=TRUE)
  dev.off()
}, error = function(e){
  warning("WARNING: failed to create profile plot")
})

#for each condition perform differential binding analysis as control
conditions <- unique(sampleSheet$Condition)
for (condition in conditions){
  #set control condition
  diffbind <- dba.contrast(diffbind, reorderMeta=list(Condition=condition))
  
  #get number of contrasts
  contrasts <- length(diffbind$contrasts)
  
  #iterate over contrasts and generate plots when DBs > 1
  for(i in 1:contrasts){
    dbs <- dba.report(diffbind, contrast=i)
    if(sum(dbs$Fold>0) > 1 | sum(dbs$Fold<0) > 1){
      contrast <- diffbind$contrasts[[i]][[5]]
      condition <- paste0(contrast[2],"_vs_",contrast[3])
      print(paste0("Generating plots for ", contrast[2]," vs ",contrast[3]))
      
      #create dir for output
      dir.out <- file.path("diffpeaks", condition)
      if(dir.exists(dir.out) == FALSE){
        dir.create(dir.out)
      }
      
      #Generate MA plots
      print("MA plots")
      file <- file.path(dir.out,paste0("MA-plot-",condition,".pdf"))
      pdf(file)
      dba.plotMA(diffbind, contrast = i)
      dev.off()
      
      #Generate volcano plots
      print("Volcano plots")
      file <- file.path(dir.out,paste0("volcano-",condition,".pdf"))
      pdf(file)
      dba.plotVolcano(diffbind, contrast = i)
      dev.off()
      
      #Plot venn diagrams
      tryCatch({
        print("Creating Venn diagram")
        file <- file.path(dir.out,paste0("venn-",condition,".pdf"))
        pdf(file)
        dba.plotVenn(diffbind, contrast=i, bDB=TRUE,bGain=TRUE,bLoss=TRUE,bAll=FALSE)
        dev.off()
      }, error = function(e){
        print("ERROR: failed to create Venn diagram")
      }, warning = function(w){
        print("WARNING: failed to create Venn diagram")
      })
      
      #export DBs to bed file
      print(paste0("Exporting differential binding sites to a bed file for ", contrast[2]," vs ",contrast[3]))
      bed.file <- file.path(dir.out,paste0("DB-",condition,".bed"))
      export.bed(dbs,con=bed.file)
      
      #prepend "chr" to chromosome names in bed file
      #system(glue("sed -i 's/^/chr/' {bed.file}"))
      
      #coverage plot
      peak <- readPeakFile(bed.file)
      file <- file.path(dir.out,paste0("coverage-",condition,".pdf"))
      pdf(file)
      covplot(peak, weightCol="V4")
      dev.off()
      
      peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=edb, annoDb=annoDb)
      df <- as.data.frame(peakAnno)
      write.csv(df, 
                file = file.path(dir.out,glue("DB-{condition}_annotated.csv")),
                row.names = FALSE)
      
      file <- file.path(dir.out,paste0("piechart-",condition,".pdf"))
      pdf(file)
      plotAnnoPie(peakAnno)
      dev.off()
      
      file <- file.path(dir.out,paste0("upsetplot-",condition,".pdf"))
      pdf(file)
      upsetplot(peakAnno, vennpie=TRUE)
      dev.off()
      
      ###Functional enrichment analysis with ReactomePA
      
      pathways_all <- enrichPathway(df$ENTREZID,
                                    organism=organism)
      
      tryCatch({
        file <- file.path(dir.out,paste0("dotplot_all_peaks-",condition,".pdf"))
        pdf(file)
        dotplot(pathways_all)
        dev.off()
      }, error = function(e){
        print("ERROR: failed to create dotplot_all_peaks")
      }, warning = function(w){
        print("WARNING: failed to create dotplot_all_peaks")
      })
      
      #add gene symbols to go.all
      go.all <- pathways_all@result
      go.all$geneSymbol <- NA
      suppressMessages(for(j in 1:nrow(go.all)){
        entrez.genes <- unlist(strsplit(go.all[j,8], "/"))
        entrez.genes <- mapIds(org.Hs.eg.db, entrez.genes,'SYMBOL','ENTREZID')
        entrez.genes <- paste(entrez.genes, collapse = "/")
        go.all[j,10] <- entrez.genes
      })
      write.csv(go.all, 
                file = file.path(dir.out,paste0("GO_all-peaks_",condition,".csv")),
                row.names = FALSE)
      
    }
  }
  
}















