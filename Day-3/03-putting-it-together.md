# Putting it together: the complete workflow

### Learning objectives:
- Review the major steps required to perform a complete DE analysis using DESeq2


### A complete workflow

The entire DESeq2 workflow essentially boils down to just a few functions run sequentially. In this lesson we will review them to consolidate our knowledge oh how to perform a complete DE analysis with DESeq2. 

Read in the count data:
```r
#Putting it all together

library(DESeq2)



#Set directory, files, and contrast
WORKING_DIRECTORY = "C:/Users/Tim Sullivan/Downloads/RNA-seq-Differential-Expression-workshop-June-2021-master/data"
COUNTS_FILE = "all_counts.txt"
METADATA_FILE = "sample_metadata.csv"
TREATMENTS=c("untreated", "Dex", "Alb", "Alb_Dex")
CONTRAST_NAME = "tx.group"
CONTRAST_BASE = "untreated"
CONTRAST_TEST = "Dex"


cts <- as.matrix(read.table(paste(WORKING_DIRECTORY,COUNTS_FILE, sep="/"),
                            sep="\t",
                            header = TRUE,
                            row.names=1,
                            stringsAsFactors = F))

colData <- read.csv(paste(WORKING_DIRECTORY,METADATA_FILE, sep="/"), row.names=1)
head(colData)

colData$tx.group <- factor(colData$tx.group, levels=TREATMENTS)


dds_matrix <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design =  as.formula(paste(" ", CONTRAST_NAME, sep="~")))

#DeSeq function includes estimateSizeFactors(), estimateDispersions(), and nbinomWaldTest()
#?DESeq
dds <- DESeq(dds_matrix)

#rld <- rlog(dds, blind = FALSE)

png(paste(CONTRAST_TEST, "vs", CONTRAST_BASE, "disp_est.png", sep="_"))
plotDispEsts(dds)
dev.off()

res <- results(dds,
               contrast = c(CONTRAST_NAME, CONTRAST_BASE, CONTRAST_TEST),
               alpha = 0.05,
               lfcThreshold = 0)


png(paste(CONTRAST_TEST, "vs", CONTRAST_BASE, "preshrink_MA.png", sep="_"))
plotMA(dds)
dev.off()

res_shrink <- lfcShrink(dds,contrast = c(CONTRAST_NAME, CONTRAST_BASE, CONTRAST_TEST), type="normal")

png(paste(CONTRAST_TEST, "vs", CONTRAST_BASE, "shrunk_MA.png", sep="_"))
plotMA(res_shrink)
dev.off()

# order results
res_shrink_ord <- res_shrink[order(res$padj),]

write.csv(as.data.frame(res_shrink_ord), file=paste(WORKING_DIRECTORY,"dex_vs_untreated_deseq.csv", sep="/"), row.names=T, quote=F )

```

While additional setups exist, for example exploratory data analysis (e.g. PCA, hierarchical clustering) or visualization of DE results (volcano plots, MA plots) the core steps of the DE analysis can be run with just these few functions detailed above (meaning it is important to understand what they are doing under the hood).
