# Putting it together: the complete workflow

### Learning objectives:
- Review the major steps required to perform a complete DE analysis using DESeq2

### Set-up

Load required R-packages:
```r
library(DESeq2)
```
ß
Load the DESeq2 dataset we already generated:
```r
dds <- readRDS("DESeq2.rdata")
```

### The complete workflow

The entire DESeq2 workflow essentially boils down to just a few functions run sequentially. In this lesson we will review them to consolidate our knowledge oh how to perform a complete DE analysis with DESeq2. 

Read in the count data:
```r
cts <- as.matrix(read.table("Day-2/all_counts.txt",
                            sep="\t",
                            header = TRUE,
                            row.names=1,
                            stringsAsFactors = F))
```

Read in the metadata:
```r
colData <- read.csv("data/sample_metadata.csv", row.names=1)
```

Construct a DESeq2 dataset from the raw counts, the metadata, and the desired design variable to be tested for differential expression.
```r
dds <- DESeqDataSetFromMatrix(countData = cts,
  colData = colData,
  design = ~ group)
```

Apply the DESeq2 analysis pipeline:
```r
dds <- DESeq(dds)
```

Perform regularized log transformation:
```r
rld <- rlog(dds, blind = FALSE)
```

Check the dispersion estimates to evaluate model fit:
```r
plotDispEsts(dds)
```

Extract, order, annotate, and subset the results from the DESeq2 object:
```r
res <- results(dds,
  name = "group_Dex_vs_untreated",
  alpha = 0.05,
  lfcThreshold = 0)

# order by adj Pval
res_ord <- res[order(res$padj),]

# add gene annotation to results
anno <- read.delim("Day-2/GRCh38.p12_ensembl-97.txt",
                   stringsAsFactors = T,
                  header = T)
anno <- anno[order(anno$Chromosome.scaffold.name),]
mat1 <- match(rownames(res_ord), anno$Gene.stable.ID)
res_ord$gene <- as.character(anno$Gene.name[mat1])

# subset results for only genes with adjusted P-values < 0.05
res_order_FDR_05 <- res_ord[res_ord$padj<0.05,]
```

Perform Empirical Bayes shrinkage of raw fold-change estimates:
```r
res_shrink <- lfcShrink(dds,
                        coef=paste0(resultsNames(dds)[which(resultsNames(dds)=="group_Dex_vs_untreated")]),
                        type="apeglm")
```

While additional setups exist, for example exploratory data analysis (e.g. PCA, hie) or visualization of DE results (volcano plots, MA plots) the core steps of the DE analysis can be run with just these few functions detailed above (meaning it is important to understand what they are doing under the hood).
