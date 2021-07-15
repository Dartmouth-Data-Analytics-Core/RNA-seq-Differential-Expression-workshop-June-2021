
# TO DO:
- incoprorate some features of bioinfo workshop statistics II section to help reduce text here 


### Exploratory data analysis & quality control

Before we run the differential expression analysis, we should explore
our dataset to learn a little more about it. Importantly, we should
determine how samples are related to each other based on their gene
expression profiles, ensure replicates cluster together, assess the
potential for any batch effects (technical variation between samples run
on different machines or by different technicians) that may exist, and
more generally build expectations for the DE analysis.

For these exploratory analyses, we use several unsupervised methods to
assess how our samples cluster relative to one another. These
unsupervised approached include **principal components analysis (PCA)**,
and **unsupervised Hierarchical clustering** and are commonly referred
to as *data reduction* appraoches.

#### Principal components analysis (PCA)

PCA is a mathematical procedure that calculates vectors that explain
varition in the dataset (in this case, variation in gene expression),
and orders samples along these vectors. We would expect samples that are
more similar to each other, e.g. replicates, to be very close to each
other along these axes of varaition, while we might expect samples in
different treatment groups to be further away.

Each vector that is calculated is called a principal component (PC), and
each principal component explains less varaition in our dataset than the
last. e.g, PC1 explains more variation in gene expression differences
between the samples than PC2. If we plot PC1 against PC2, samples will
‘cluster’ with other samples that have similar gene expression profiles,
and be further away from samples with more distant expression profiles.

<center>
![](../figures/pca_example.png)
</center>

Again, StatQuest has an excellent
[video](https://www.youtube.com/watch?v=_UVHneBUBW0) that explains the
fundamental concepts of PCA, and provoides more details how to the PCs
themselves are calculated.

To perform mathematical procedures such as PCA, it is better to work
with a transformed version of the counts, rather than the original
untransformed values. DESeq2 actually utilizes its own transformation
procedure, called the **regularized logatrithm (rlog)** implemented in
the `rlog()` function. The rlog is similar in principle to a standard
log transformation of the data, but is able to more appropriately
transform the counts for genes with low expression values.

DESeq2 is also capable of implementing a **variance stabilising
transformation (VST)** for count data, this is generally recommended for
larger datasets due to increased speed.

For this analysis, we will use the rlog, which produces values on the
log2 scale that are also normalized for library size during the rlog
procedure. Lets perform the rlog transformation on our data.

```r
    rld <- rlog(dds, blind = FALSE)
    head(assay(rld))
```

We can illustrate the benefit of using the rlog over standard log
transformation (+ a pseudo-count for genes with 0 counts where the log
of 0 is infinity) by comparing the transformed values for two samples
against each other.

```r
    par(mfrow=c(1,2))
    plot(log2(cts[,1]+1), log2(cts[,2]+1), col = "cornflowerblue", xlab = "Sample 1", ylab = "Sample 2", main = "Log2 + 1")
    plot(assay(rld)[,1], assay(rld)[,2], col = "indianred", xlab = "Sample 1", ylab = "Sample 2", main = "rlog")
```

We can use these transformed values to investigate how many features
(genes) in our dataset exhibit variability across samples. This is
useful to know as we only want to use variable features for PCA. Genes
that don’t explain any variation in the dataset aren’t useful for
helping us explore differences between the samples.

```r
    # calculate gene expression level variance between samples
    var <- rev(rowVars(assay(rld))[order(rowVars(assay(rld)))])
    # plot variance for genes accross samples
    plot(var, las = 1, main="Sample gene expression variance", xlab = "Gene", ylab = "Variance")
    abline(v=1000, col="red") ; abline(v=500, col="green") ; abline(v=250, col="blue")
```

At around 500 the variance starts to spike upwards, so this is the
number of variable features (genes) to use. Lets restrict the dataset to
500 genes for purposes of the PCA. Now lets extract and visualize the
variance explained by each PC to determine which are most informative.

```r
    # modify variable feature number to be used in PCA and hierachical clutering based on no. of most variable features
    var_feature_n <- 500
    # perform PCA and order by variance
    rv <- rowVars(assay(rld))
    select <- order(rv, decreasing = TRUE)[seq_len(min(var_feature_n, length(rv)))]
    pca <- prcomp(t(assay(rld)[select, ]))
    # extract the varioance explained by each PC
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    names(percentVar)[1:5] <- c("PC1", "PC2", "PC3", "PC4", "PC5")
    percentVar <- percentVar[1:5]
    # plot variance for top 10 PCs
    barplot(percentVar[1:5], col = "indianred", las = 1, ylab = "% Variance", cex.lab = 1.2)
```

We can see that the majority of variance is explained by the first few
PCs, therefore visualizing where samples fall along these PCs will be
the most informative way to identify major differences between them,
based on their gene expression profiles. Lets generate a PCA plot for
PC1 vs PC2.

```r
    # construct data frame w/ PC loadings and add sample labels
    pca_df <- as.data.frame(pca$x)
    pca_df$tx.group <- dds@colData$tx.group
    pca_df$sample_ids <- colnames(dds)
    # add colors for plotting to df
    pca_df$col <- NA
    for(i in 1:length(levels(pca_df$tx.group))){
      ind1 <- which(pca_df$tx.group == levels(pca_df$tx.group)[i])
      pca_df$col[ind1] <- i
    }
    # plot PC1 vs PC2
    plot(pca_df[, 1], pca_df[, 2],
         xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"),
         ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
         main=paste0("PC1 vs PC2 for ", var_feature_n, " most variable genes"),
         pch=16, cex=1.35, cex.lab=1.3, cex.axis = 1.15, las=1,
         panel.first = grid(),
         col=pca_df$col)
    text((pca_df[, 2])~(pca_df[, 1]), labels = pca_df$tx.group, cex=0.6, font=2, pos=4)
```

Looking at the plot, we can see some clear clustering by treatment group
for some of the samples. For example, **untreated** samples appear to
have much lower PC1 values than all the **Dex** samples, suggesting that
some of the largest variability in gene expression differences between
samples in this dataset explains differences between **untreated** and
**Dex** treated samples. Therefore we expect the most substantial
differential expression to be found for the **untreated** vs **Dex**
analysis.

The **Alb** treated samples, and the co-treated samples **(Alb\_Dex)**
do not seem to consistently cluster along PC1 and PC2. One explanation
for this could be that treatment with **Alb** or co-treatment with **Alb
+ Dex** have inconsistent effects on gene expression across study
subjects. This is often the case with *in vivo* data sets, rather than
those collected in *in vitro* systems such as cell culture. This could
indicate that a greater number of replicates is required to accurately
assess the impact of such treatment on gene expression, as there are
usually a greater number of uncontrollable variables in *in vivo*
datasets (e.g. genotype).

Another explanation could be there is another variable, unrelated to
sample groupings by treatment, that explains a lot of the variability in
the dataset, such as a batch effect, that is driving the clustering
observed in the PCA plot. For example, if the two Alb treated samples
with high values for PC1 and PC2 were known to have been processed in a
different batch to the Alb treated sample with lower values for these
PCs, we could infer that there is a batch effect that is driving the
differences between these samples. Unfortunately, no batch information
was provided for this dataset, so we cannot confidently attribute the
lack of clustering to batch effects.

It is interesting that 4 of the samples, 1 from each treatment group,
have much higher PC2 values than all other samples. This is a pattern we
might expect to see if these 4 samples were processed in a separate
batch, and the gene expression variation captured by PC2 is related to
the batch effect, rather than our biological factor of interest
(treatment group). For purposes of this example, lets create a fake
variable for sample batch, and include this on the PCA plot.

```r
    pca_df$batch <- NA
    pca_df$batch[pca_df$PC2 < 10] <- "Batch 1"
    pca_df$batch[pca_df$PC2 > 10] <- "Batch 2"
    pca_df$batch <- factor(pca_df$batch, levels = c("Batch 1", "Batch 2"))
    # replot PC1 vs PC2 with the shape set to batch
    plot(pca_df[pca_df$batch=="Batch 1", 1], pca_df[pca_df$batch=="Batch 1", 2],
         xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"),
         ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
         main=paste0("PC1 vs PC2 for ", var_feature_n, " most variable genes"),
         pch=16, cex=1.35, cex.lab=1.3, cex.axis = 1.15, las=1,
         panel.first = grid(),
         col=pca_df$col,
         ylim=c(-14,20))
    points(pca_df[pca_df$batch=="Batch 2", 1], pca_df[pca_df$batch=="Batch 2", 2],
           col=pca_df$col,
           pch=2, cex=1.35)
    legend(9.5, 10.5, levels(pca_df$batch), pch = c(16, 2))
    legend(1.5, 11.5, levels(pca_df$tx.group), pch = 16, col = pca_df$col)
    text((pca_df[, 2])~(pca_df[, 1]), labels = pca_df$tx.group, cex=0.6, font=2, pos=4)
```

If this was a real batch variable, this would clearly indicate a batch
effect. If batch effects are not accounted for in the differential
expression analysis, these samples will increase the per gene variance
(and therefore dispersion paracmeter in the negative-binomial model) and
prevent the detection of genes that are truly differentially expressed
between our experimental conditions of interest, effectively reducing
statistical power. This could lead to both false-negatives and
false-positives.

If you detect a batch effect in your data, you can try to: \* use a
statistical procedure to remove the batch effect from your data \*
adjust for batch as a term in your statistical model when running
differential expression \* remove the samples driving the batch effect,
if you have a particular reason to suspect these specific samples

How you handle a batch effect is a complicated issue, and is largely
dependent on the extent of the batch effect. If the batch effect is very
large, it may be too diffcult to effectively remove it statistically, or
regress out variation attributable to it in the DE analysis. This is
where practicing your protcol and confirming you can get consistent
results across replciates and batches comes in. If you do this work
ahead of time, you reduce the risk of having to deal with this
complicated batch effect in the analysis.

**Take home message on batch effects:** If your experiment includes
multiple batches, you should always include them in your unsupervised
analyses to check for a batch effect.

#### Hierarchical clustering

Hierarchical clustering is another complimentary approach to explore the
relationships between your samples. While supervised clustering
approaches exist, we will perform an unsupervised analysis so that we do
not impose any restrictions on the clustering of the samples.

Hierachical clustering is often associated with **heatmaps**, as it is a
useful way to explore the results of hierachical clustering. Here we
represent genes are rows, and individual samples as columns. The
**denrograms ** on the rows and the columns represent the *‘distances’*
calculated between each of the genes/samples. Presenting the data in
this way is useful as it allows us to identify samples whose patterns of
gene expression are similar to each other, but also modules of genes
that change in a similar way across our samples, and may share some
common function of interest.

<center>
![](../figures/heatmaps.png)
</center>

The first step in a hierachical clustering analaysis is to *scale your
data*. This means that expression levels are all transformed onto the
same scale before clustering. This is important to do as we can only
visualize so many colors at once, and a very highly expressed gene would
mean that all the other genes would essentially invisible on this scale.
Scaling for clustering in this way is typically performed by calculating
Z-scores, where the mean for each gene across all the samples is
subtracted from each of the individual expression values for each gene,
this centers the expression values around 0. We then divide these values
by the standard deviation to make sure the data is more tightly grouped,
and we can represent lots of genes in the same scale.

Although we will not go into full detail here on how the actual
clustering algorithm works to group samples and genes, once more
**StatQuest** has an [excellent
video](https://www.youtube.com/watch?v=oMtDyOn2TCc) on this topic.

Similarly to the PCA, we perform the clustering using the **rlog
transformed data** and the **500 most variable features**, as features
that do not vary across samples are not informative for dimension
reduction appraoches.

```r
    # select top X no. of variable genes
    topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), var_feature_n)
    # set up gene expression matrix
    mat1 <- assay(rld)[topVarGenes,]
    # scale matrix by each col. values
    mat_scaled = t(apply(mat1, 1, scale))
    # set up colors for heatmap
    col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
    cols1 <- brewer.pal(11, "Paired")
    cols2 <- brewer.pal(9, "Greens")
    # set up annotation bar for samples
    ha1 = HeatmapAnnotation(Group = colData(dds)$tx.group,
                            col = list(Group = c("untreated" = cols1[1], "Dex" = cols1[2],
                                                 "Alb" = cols1[5], "Alb_Dex" = cols1[6])),
                                       show_legend = TRUE)
    # se up column annotation labels (samples)
    ha = columnAnnotation(x = anno_text(colData(dds)$SRR,
                                        which="column", rot = 45,
                                        gp = gpar(fontsize = 10)))
    # generate heatmap object
    ht1 = Heatmap(mat_scaled,
                  name = "Expression",
                  col = col,
                  top_annotation = c(ha1),
                  bottom_annotation = c(ha),
                  show_row_names = FALSE)
    # plot the heatmap
    draw(ht1, row_title = "Genes", column_title = "Top 500 most variable genes")
```

As we saw in the PCA, the Alb and co-treated samples do not form any
clear clusters. We may want to remove them and perform the clustering
again so that we can compare the untreated and Dex samples more easily.

```r
    ind_to_keep <- c(which(colData(rld)$group=="untreated"), which(colData(rld)$group=="Dex"))
    topVarGenes <- head(order(rowVars(assay(rld)[,ind_to_keep]), decreasing=TRUE), var_feature_n)
    # set up gene expression matrix
    mat1 <- assay(rld)[topVarGenes, ind_to_keep]
    # scale matrix by each col. values
    mat_scaled = t(apply(mat1, 1, scale))
    # set up colors for heatmap
    col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
    cols1 <- brewer.pal(11, "Paired")
    cols2 <- brewer.pal(9, "Greens")
    # subset coldata for samples in untx and ex groups
    colData_sub <- colData(dds)[ind_to_keep, ]
    # set up annotation bar for samples
    ha1 = HeatmapAnnotation(Group = colData_sub$tx.group,
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
    draw(ht1, row_title = "Genes", column_title = "Top 500 most variable genes")
```

There does indeed seem to be relatively good clustering between
untreated and Dex samples, suggesting there are unique gene expression
programs defining the Dex samples from the untreated. However, one of
these samples in the Dex group seems to be clustered further away from
the other Dex samples. This could be the Dex treated sample that
clustered away from the other Dex treated samples on the PCA. We can add
an annotation bar for the fake batch effect we created earlier to this
plot to confirm this.

```r
    # which samples had values > 10 for PC2
    pca_df$sample_ids[pca_df$PC2 > 10 & pca_df$tx.group=="untreated"]

    pca_df$sample_ids[pca_df$PC2 > 10 & pca_df$tx.group=="Dex"]

    # set the batch variable for these samples as batch 2
    colData_sub$batch <- "Batch 1"
    colData_sub$batch[colData_sub$SRR=="SRR1039516"] <- "Batch 2"
    colData_sub$batch[colData_sub$SRR=="SRR1039517"] <- "Batch 2"

    # set up annotation bar for samples
    ha1 = HeatmapAnnotation(group = c(as.character(colData_sub$tx.group)),
                            batch = c(as.character(colData_sub$batch)),
                            col = list(group = c("untreated" = cols1[1], "Dex" = cols1[2]),
                                       batch = c("Batch 1" = cols1[5], "Batch 2" = cols1[6])),
                                       show_legend = TRUE)
    # generate heatmap object
    ht2 = Heatmap(mat_scaled, name = "Expression", col = col,
                  top_annotation = c(ha1),
                  bottom_annotation = c(ha),
                  show_row_names = FALSE)
    # plot the heatmap
    draw(ht2, row_title = "Genes", column_title = "Top 500 most variable genes")
```

Based on our newly labeled plot it does seem that these 2 samples are
outliers based on the hierachical clustering, which would support the
presence of a batch effect if these data were infact collected in
multiple batches. Again without the metadata to indicate this is a batch
effect using a method to correct for a batch effect is inappropriate. In
this case I would reach out to the seqeuncing center or authors of the
paper where the data was published and ask for additional metadata to
confirm our suspicion about batch effects in this data.

------------------------------------------------------------------------

That wraps up exploratory analysis. In the next R markdown, we will
perform the differential expression analysis.

Session Information
-------------------
```r
    sessionInfo()
```
