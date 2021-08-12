To do:
- add learning objectives
- confirm code works with current directory structure
- adapt text so that this markdown reads appropriately as a standalone section
- reduce overall text where possible
- change paths to figures folder so that figures are shown in markdown


# Results annotation & visualization


```r
library(vsn)
library(dplyr)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(readr)
library(circlize)
library(EnhancedVolcano)
library(apeglm)
library(xtable)
library(kableExtra)
```

We want to also add the annotation data for each gene (symbol, genome coordinates, etc.) to the results. Since we used Ensembl version 97 to annotate these data, we need to use the Ensembl 97 annotation data to annotate these results. We can obtain this for our species of interest in a flat file format using the [BioMart on the Ensembl website](http://uswest.ensembl.org/biomart/martview/b0399bb192186dea3aedf87d82a4580c).

```r
    # read in the flat file we downloaded and have a look at it
    anno <- read.delim("data/GRCh38.p12_ensembl-97.txt", stringsAsFactors = T, header = T)
    anno <- anno[order(anno$Chromosome.scaffold.name),]
    dim(anno)


    # have a look at the first few rows
    head(anno)
```

Lets have a look at the Chromosome distribution of features

```r
    tab1 <- table(anno$Chromosome.scaffold.name)
    tab1[1:22]
```

Lets also quickly check that nothing is duplicated in the ENSG ID column of our annotation, as this would cause problems when merging with our results.

```r
    any(duplicated(anno$Gene.stable.ID))
```

Now lets add the annotation for each gene name directly to the results.

```r
    # use match() to find corresponding indicies (rows) for each ENSG ID
    mat1 <- match(rownames(res_ord), anno$Gene.stable.ID)
    table(is.na(mat1))

    # add gene names to results as a new column
    res_ord$gene <- as.character(anno$Gene.name[mat1])
    head(res_ord, 20)

```

Lets also add some other columns that might be of interest to us when reviewing the results.

```r
    res_ord$chr <- as.character(anno$Chromosome.scaffold.name[mat1])
    res_ord$start <- as.character(anno$Gene.start..bp.[mat1])
    res_ord$end <- as.character(anno$Gene.end..bp.[mat1])
    res_ord$strand <- as.character(anno$Strand[mat1])
```

------------------------------------------------------------------------

### Visualization of Differential Expression

#### Volcano plot

Volcano plots are a useful visualization for exploring your results, the **log2 fold change** (x-axis) is plotted against the **-log10 P-value**. Since the -log10() of a really small number is a very large value, any gene that has a very small P-value and was significantly differentially expressed, will appear higher up along the y-axis. In contrast, the -log10 of 1 (`-log10(1)`) is equal to `0`, therefore genes with low statistical significance (P-values approaching 1) will appear lower down on the y-axis.

Similarly, genes with larger fold changes will appear further along the x-axis, in both directions. Genes with a positive fold change represent genes whose expression was greater than the group of the experimental design variable used as baseline, while genes with a negative fold change represent genes whose expression was lower than in the baseline group.

The fold-change value of genes with non-significant fold changes is not meaningful, as there is not enough statistical confidence in these fold changes.

```r
    plot(res$log2FoldChange, -log10(res$pvalue),
         main = "Volcano plot",
         las = 1, col = "indianred",
         ylab = "- log10 P-value", xlab = "log2 Fold change")

    # add horizontal lines to help guide interpretation
    abline(h=-log10(0.05/nrow(res)), lty = 2, col = "black") # Bonferonni
    abline(h=-log10(0.05), lty = 2, col = "black") # nominal P-value
```

Here we can clearly see that there are quite a few genes above our significance threshold in both the up and downregulation directions (negative and positive fold changes), that also have absolute log2 fold change values of at least 2 or more. Of particular interest, there seem to be a few genes with very large fold change values & -log10 P-values, making them especially interesting as their effect size is large AND our confidence in this fold change is good.

It is a little hard to make specific inferences from this plot at the individual gene level, so some labels for interesting data points ( and some colors) would definitely improve this volcano plot, and make it more informative. We will use the **ggpolot2** R package to do this, and we will color each point based on a combination of fold change and P-value, as these determine which genes are of most interest to us.

```r
    # save a dataframe from the results() output
    res_tmp <- as.data.frame(res_ord)

    # add a column that will be used to save the colors we want to plot
    res_tmp$cols <- c()

    # set the significance cut off (alpha) and fold change threshold to be used for coloring of genes
    alpha <- 0.05/nrow(res)
    fc_cutoff <- 2

    # loop through our dataframe and add values to the color column based on magnitude of alpha and LFCs
    res_tmp$cols <- NA
    for(i in 1:nrow(res_tmp)){
        if(is.na(res_tmp$pvalue[i])){
          res_tmp$cols[i] <- NA
        }
        else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i] > fc_cutoff){
          res_tmp$cols[i] <- "indianred"
        }
        else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i] < -fc_cutoff){
          res_tmp$cols[i] <- "indianred"
        }
        else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i]>-fc_cutoff & res_tmp$log2FoldChange[i]<fc_cutoff){
          res_tmp$cols[i] <- "cornflowerblue"
        }
        else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] > fc_cutoff){
          res_tmp$cols[i] <- "gray47"
        }
        else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] < -fc_cutoff){
          res_tmp$cols[i] <- "gray47"
        }
        else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] < fc_cutoff){
          res_tmp$cols[i] <- "gray10"
        }
    }

    res_tmp$ENSG <- rownames(res_tmp)

    # generate the splot
    p = ggplot(res_tmp, aes(log2FoldChange, -log10(pvalue))) +
        geom_point(aes(col=col), alpha = 0.5, size =2.5, colour = res_tmp$cols, fill = res_tmp$cols)  +
        xlab("Log2 fold change") + ylab("-log10 Q-value") +
        ylim(0, 9) +
        xlim(-5, 11) +
        geom_hline(yintercept = -log10(alpha), color = "black", linetype = "dashed", size = 0.4) +
        theme(legend.key = element_blank()) +
        ggtitle("Control vs Dex")

    # print the plot
    print(p)

```

This is nice, but some labels for potentially interesting genes would be useful. Lets add some using the **ggrepel** package.

```r
    p2 <- p +
      # add labels to genes w/ LFC > 2 and above alpha threshold
      geom_label_repel(data = subset(res_tmp, log2FoldChange > 2 & pvalue < alpha), aes(label = gene),
                         box.padding   = 0.35,
                         nudge_x = 0.1,
                         nudge_y = 0.1,
                         point.padding = 1,
                         label.size = 0.1,
                         segment.size = 0.3,
                         segment.color = 'grey50', size = 3) +
      # add labels to genes w/ LFC < -2 and above alpha threshold
      geom_label_repel(data = subset(res_tmp, log2FoldChange < -2 & pvalue < alpha), aes(label = gene),
                         box.padding   = 0.35,
                         nudge_x = -0.1,
                         nudge_y = 0.1,
                         point.padding = 1,
                         label.size = 0.1,
                         segment.size = 0.3,
                         segment.color = 'grey50', size = 3) +
      # add vertical fold change lines
      geom_vline(xintercept = fc_cutoff, colour = "black", linetype="dotted") +
      geom_vline(xintercept = -fc_cutoff, colour = "black", linetype="dotted")

    # print the plot
    print(p2)
```

This looks a lot better, and gives us a lot more information than the first, very basic plot we generated.

Food for thought: detecting truly differentially expressed genes is dependent on the technical variance between your replicates. If the technical variance is high, you generally need a large fold-change to achieve statistical significance. The more replicates you have, the more you are able to reduce this technical variance, which increases your statistical power, and enables you to confidently detect differential expression of smaller fold changes. For example, for an experiment where there are 300 truly differentially expressed genes between your conditions, you may detect 200 of these with 3 replicates, while you may detect 250 with 5 replicates.

**Save our results to .csv files**

```r
    # subset @ 5% adjusted pval sig. level
    res_order_FDR_05 <- res_ord[res_ord$padj<0.05,]
    nrow(res_order_FDR_05)

    # write both to csv files
    write.csv(as.data.frame(res_ord), file= "DE_results.csv")
    write.csv(as.data.frame(res_order_FDR_05), file="DE_results.FDR.0.05.csv")
```
-----------------------

### Other visualizations - MA plots

MA plots are also useful ways to visualize results from a DE analysis of
RNA-seq data. These involve plotting the log2 fold-change (the so called
M-value, representing the *M* in *MA-plot*) against the average
expression level of a gene (the *A* in *MA-plot*).

The MA-plot allows us to inspect the **full range of expression values
over which we detected significant DEGs, and what the magnitude of these
fold-changes is**. In a typical experiment, we expect to see DEGs across
most of the range of expression values. To help identify genes that were
significantly DE, any gene with an adjusted P-value of &lt; 0.05 (or
whatever threshold is set) is colored in red.

```r
    plotMA(res_ord, ylim=c(-6,6), main = "Raw Log2 Fold change")
```

The **log2 fold-change** plotted above is the raw LFC value estimated by
the negative binomial GLM that we used in modeling. However, as we
discussed above, the individual estimates of variance or dispersion for
a single gene are often unreliable, and this holds true
`log2 fold change` also.

**To obtain more useful LFC estimates,** `DESeq2` performs a statsitical
procedure that involves **shrinking the raw fold change estimates toward
zero** for genes that are less likely to contain reliable or highly
important information.

This is done in a very similar way to the shrinkage using empirical
bayes that we discussed for the **dispersion estimates**.

**For shrinking LFC values, LFCs are penalized for properties such as:
**  
- low count values  
- high dispersion (& thus reduced confidence in expression levels)

DESeq2 provides a function `lfcShrink()` that must be implemented
separately of the standard workflow implemented using `DESeq2()`.

```r
    # calculate shrunken fold change estimate
    res_shrink <- lfcShrink(dds,
                        coef=paste0(resultsNames(dds)[which(resultsNames(dds)=="group_Dex_vs_untreated")]),
                        type="apeglm")
```

After performing the shrinkage procedure, we compare the raw and
shrunken LFCs to assess the impact of shrinkage.

**Raw estimates of log2 FC:**

```r
    plotMA(res_ord, ylim=c(-6,6), main = "Raw Log2 Fold change")
```

**Shrunken estimates of log2 FC:**

```r
    plotMA(res_shrink, ylim=c(-6,6), main = "Shrunken Log2 Fold change")
```

We can see that **significantly DE genes are detected across the full
range of expression values** (x-axis), which is a good sign that our
differential expression modeling has worked well. We can also see that
we have a handful of genes with larger expression values (&gt; LFC 2)
which potentially represent the most important individual genes, while
the majority of our DEGs have a LFC &lt; 1.5 (ish).

Comparing to the raw LFCs, we can also see that the **majority of genes
with lower expression values have have their LFCs shrunk toward zero**.
This is important as genes with low counts may simply end up with a
large LFC since this is easy to do at small count values, but these are
unlikely to be accurate fold-changes, so we don’t want to prioritize
their importance by giving them a large LFC.

It’s always good to look at the shrunken estimates, to confirm that you
don’t have a lot of DEGs at very small count values. If you do, you may
want to look at the expression levels for those genes to investigate
these findings in more detail.

**As the mean or counts increase, it is evident that the level of
shrinkage is less**, although may still be high for genes with greater
dispersion estimates. As we move toward the more highly expessed genes,
you can see how more genes at lower fold change values are able to be
identified as significant, which is due to the fact that there is more
information avaiable for these genes, so we can be more confident during
hypothesis tetsing of these genes.

**Note:** This shrinkage does not really change the hypothesis testing,
therefore is performed independently, as is for use in prioritizing your
results further for visual inspection or some sort of functional
analysis (e.g. pathway analysis).

------------------------------------------------------------------------

#### Hierachical clustering on the DEGs

A final visualization that is useful to generate is a heatmap based on
unsupervised hierachical clustering of the DEGs identified. We can do
this by limiting the matrix of rlog values to only those for the DEGs,
and then performing the clustering specifically on these data.

```r
    rld <- rlog(dds, blind = FALSE)
    ind_to_keep <- c(which(colData(rld)$group=="untreated"), which(colData(rld)$group=="Dex"))

    # set up gene expression matrix
    mat1 <- assay(rld)[rownames(res_order_FDR_05), ind_to_keep]

    # scale matrix by each col. values
    mat_scaled = t(apply(mat1, 1, scale))

    # set up colors for heatmap
    col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
    cols1 <- brewer.pal(11, "Paired")
    cols2 <- brewer.pal(9, "Greens")

    # subset coldata for samples in untx and ex groups
    colData_sub <- colData(dds)[ind_to_keep, ]

    # set up annotation bar for samples
    ha1 = HeatmapAnnotation(Group = colData_sub$group,
                            col = list(Group = c("untreated" = cols1[1], "Dex" = cols1[2])),
                                       show_legend = TRUE)

    # se up column annotation labels (samples)
    ha = columnAnnotation(x = anno_text(colData_sub$SRR,
                                        which="column", rot = 45,
                                        gp = gpar(fontsize = 10)))

    # generate heatmap object
    ht1 = Heatmap(mat_scaled, name = "Expression", col = col,
                  top_annotation = c(ha1),
                  bottom_annotation = c(ha),
                  show_row_names = FALSE)

    # plot the heatmap
    draw(ht1, row_title = "Genes", column_title = "Hierachical clustering of DEGs (padj<0.05)")
```
