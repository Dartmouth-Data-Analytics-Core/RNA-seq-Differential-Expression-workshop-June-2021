# Day3 - 03 putting it all together

library(DESeq2)

# Read in the count data:
cts <- as.matrix(read.table("data/all_counts.txt",
                            sep="\t",
                            header = TRUE,
                            row.names=1,
                            stringsAsFactors = F))

# Read in the metadata:
colData <- read.csv("data/sample_metadata.csv", row.names=1)

# convert tx.group into factor class variable
colData$tx.group <- factor(colData$tx.group, levels=c("untreated", "Dex", "Alb", "Alb_Dex"))

# Construct a DESeq2 dataset from the raw counts, the metadata, and the desired design variable to be tested for differential expression.
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ tx.group)

# Apply the DESeq2 analysis pipeline:
dds <- DESeq(dds)

# Perform regularized log transformation:
rld <- rlog(dds, blind = FALSE)

# Check the dispersion estimates to evaluate model fit:
plotDispEsts(dds)

# extract results
res <- results(dds,
               name = "tx.group_Dex_vs_untreated",
               alpha = 0.05,
               lfcThreshold = 0)

# order by adj Pval
res_ord <- res[order(res$padj),]

# add gene annotation to results
anno <- read.delim("data/GRCh38.p12_ensembl-97.txt",
                   stringsAsFactors = T,
                   header = T)
anno <- anno[order(anno$Chromosome.scaffold.name),]
mat1 <- match(rownames(res_ord), anno$Gene.stable.ID)
res_ord$gene <- as.character(anno$Gene.name[mat1])

# subset results for only genes with adjusted P-values < 0.05
res_order_FDR_05 <- res_ord[res_ord$padj<0.05,]


#Perform Empirical Bayes shrinkage of raw fold-change estimates:
res_shrink <- lfcShrink(dds,
                        coef=paste0(resultsNames(dds)[which(resultsNames(dds)=="group_Dex_vs_untreated")]),
                        type="apeglm")
