# Putting it together: the complete workflow

### Learning objectives:
- Review the major steps required to perform a complete DE analysis using DESeq2
- Create a baseline for a DESeq2 workflow

### A complete workflow

The entire DESeq2 workflow essentially boils down to just a few functions run sequentially. In this lesson we will review them to consolidate our knowledge oh how to perform a complete DE analysis with DESeq2. 

```r
#Putting it all together

#Import packages
library(DESeq2)

#Set directory, file names, and contrasts
WORKING_DIRECTORY = "C:/Users/Tim Sullivan/Downloads/RNA-seq-Differential-Expression-workshop-June-2021-master/data"
setwd(WORKING_DIRECTORY)
COUNTS_FILE = "all_counts.txt"
METADATA_FILE = "sample_metadata.csv"
TREATMENTS=c("untreated", "Dex", "Alb", "Alb_Dex")
CONTRAST_NAME = "tx.group"
CONTRAST_BASE = "untreated"
CONTRAST_TEST = "Dex"

#Load Counts and Metadata
cts <- as.matrix(read.table(paste(WORKING_DIRECTORY,COUNTS_FILE, sep="/"), sep="\t", header = TRUE, row.names=1, stringsAsFactors = F))
colData <- read.csv(paste(WORKING_DIRECTORY,METADATA_FILE, sep="/"), row.names=1)

#Build DESeq objects
colData$tx.group <- factor(colData$tx.group, levels=TREATMENTS)
dds_matrix <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design =  as.formula(paste(" ", CONTRAST_NAME, sep="~")))

#DeSeq function includes estimateSizeFactors(), estimateDispersions(), and nbinomWaldTest(). ?DESeq for more information.
dds <- DESeq(dds_matrix)

#Exploratory Plotting
rld <- rlog(dds, blind = FALSE)
var_feature_n <- 500
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[1:500]
rld_sub <- assay(rld)[select, ]
rld_sub <- t(rld_sub)
pca <- prcomp(rld_sub)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- percentVar[1:5]
names(percentVar) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
png(paste(CONTRAST_TEST, "vs", CONTRAST_BASE, "PCA_Variance_Explained.png", sep="_"))
barplot(percentVar, col = "indianred", las = 1, ylab = "% Variance", cex.lab = 1.2)
dev.off()
pca_df <- as.data.frame(pca$x)
pca_df$tx.group <- dds@colData$tx.group
pca_df$sample_ids <- colnames(dds)
pca_df$col <- NA
for(i in 1:length(levels(pca_df$tx.group))){
  ind1 <- which(pca_df$tx.group == levels(pca_df$tx.group)[i])
  pca_df$col[ind1] <- i
}

png(paste(CONTRAST_TEST, "vs", CONTRAST_BASE, "PCA_1_2.png", sep="_"))
plot(pca_df[, 1], pca_df[, 2],
     xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"),
     ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
     main=paste0("PC1 vs PC2 for ", var_feature_n, " most variable genes"),
     pch=16, cex=1.35, cex.lab=1.3, cex.axis = 1.15, las=1,
     panel.first = grid(),
     col=pca_df$col)
text(pca_df[, 1], pca_df[, 2], labels = pca_df$tx.group, cex=0.6, font=2, pos=4)
dev.off()

#Plot Mean-Variance relationship -- Look at this!
png(paste(CONTRAST_TEST, "vs", CONTRAST_BASE, "disp_est.png", sep="_"))
plotDispEsts(dds)
dev.off()

#Gather significant results
res <- results(dds, contrast = c(CONTRAST_NAME, CONTRAST_BASE, CONTRAST_TEST), alpha = 0.05, lfcThreshold = 0)

#Pre-shrinkage MA Plot
png(paste(CONTRAST_TEST, "vs", CONTRAST_BASE, "preshrink_MA.png", sep="_"))
plotMA(dds)
dev.off()

#Shrinkage
res_shrink <- lfcShrink(dds,contrast = c(CONTRAST_NAME, CONTRAST_BASE, CONTRAST_TEST), type="normal")

#Post-shrinkage MA Plot
png(paste(CONTRAST_TEST, "vs", CONTRAST_BASE, "shrunk_MA.png", sep="_"))
plotMA(res_shrink)
dev.off()

#Order results
res_shrink_ord <- res_shrink[order(res$padj),]

#Annotate results with gene name
annotation_table <- read.delim("GRCh38.p12_ensembl-97.txt", stringsAsFactors = T, header = T)
annotation_matrix <- match(rownames(res_shrink_ord), annotation_table$Gene.stable.ID)
res_shrink_ord$gene <- as.character(annotation_table$Gene.name[annotation_matrix])

#Write output CSV
write.csv(as.data.frame(res_shrink_ord), file=paste(WORKING_DIRECTORY,"dex_vs_untreated_deseq.csv", sep="/"), row.names=T, quote=F )

```

Despite DESeq2 being able to be run in an automated way, it is important to understand the functions and calculations within.  If something doesn't look right, you now have the skills and knowledge to be able to pull apart the counts matrix, dds, and results objects.  Additionally, plots of heatmaps based on counts, subsetted by variance, or subsetted by the top differentially expressed genes can be useful.  The above lines of code should be used as a starting point for a fully explored and annotated differential expression analysis.
