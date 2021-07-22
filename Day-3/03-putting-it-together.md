Recap of the full workflow
--------------------------

In the above analysis we got into the details of how the statistics
behind differential expreission work and ran some extra code to
demonstrate the utility of these statistics. However, the entire DESeq2
workflows boils down to just a few functions run sequentially. Lets do a
quick recap of these functions to help consolidate what we have learnt.

Read in the data:

```r
    cts <- as.matrix(read.table("Day-2/all_counts.txt",
                                sep="\t", header = TRUE, row.names=1,
                                stringsAsFactors = F))
```

Read in the metadata:

```r
    sra_res <- read.csv("Day-2/sra_result.csv", row.names=1)
    sra_res$Sample <- sra_res$Sample.Accession
    sra_run <- read.csv("Day-2/SraRunInfo.csv", row.names=1)
```

Construct a DESeq2 dataset from the raw counts, the metadata, and the
desired design variable to be tested for differential expression.

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

Use rlog to perform exploratory analyses: - Principal components
analysis (PCA) - Unsupervised hierachical clustering

Check the disperion estimates to evalute model fit:

```r
    plotDispEsts(dds)
```

Extract, order, annotate, and subset the results from the DESeq2 object

```r
    res <- results(dds,
      name = "group_Dex_vs_untreated",
      alpha = 0.05,
      lfcThreshold = 0)

    # order by adj Pval
    res_ord <- res[order(res$padj),]

    # add gene annotation to results
    anno <- read.delim("Day-2/GRCh38.p12_ensembl-97.txt", stringsAsFactors = T, header = T)
    anno <- anno[order(anno$Chromosome.scaffold.name),]
    mat1 <- match(rownames(res_ord), anno$Gene.stable.ID)
    res_ord$gene <- as.character(anno$Gene.name[mat1])

    # subset results for only genes with adjusted P-values < 0.05
    res_order_FDR_05 <- res_ord[res_ord$padj<0.05,]
```

Perform empirical bayes shrinkage of raw fold-change estimates:

```r
    res_shrink <- lfcShrink(dds,
                        coef=paste0(resultsNames(dds)[which(resultsNames(dds)=="group_Dex_vs_untreated")]),
                        type="apeglm")
```

Generate visualizations: - MA plots (raw vs shrunken fold-changes) -
Volcano plots - Heatmaps

Session Information
-------------------
```r
    sessionInfo()
```
