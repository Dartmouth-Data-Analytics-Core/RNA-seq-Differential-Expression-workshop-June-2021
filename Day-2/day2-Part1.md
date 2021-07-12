*Bulk RNA-seq data analysis data analysis workshop, July 2020*

Exploratory data analysis in R
------------------------------

### Introduction and set-up

Several popular R packges designed for exploration and statistical
analysis of bulk RNA-seq data exist, including
[*EdgeR*](https://www.bioconductor.org/packages/release/bioc/html/edgeR.html),
[*limma-voom*](http://bioconductor.org/packages/release/bioc/html/limma.html),
[*DESeq2*](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
For the purposes of this workshop, we will use DESeq2 to perform the
parts of the analysis, including reading in the data, normalization of
read counts, and fitting statistical models to test differential
expression. [Detailed
tutorials](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
for using DESeq2 can be found on its Bioconductor page.

DESeq2 is a very well organized package that applies robust algorithms
to perform several aspects of RNA-seq data analysis. If you plan to use
DESeq2 for your work, you should read both the tutorials made available
on their Bioconductor page, and the original manuscript for DESeq2, in
[Love *et al*,
2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
to develop an understanding of the theory behind DESeq2 and the
proccesses implemented by its functions.

Despite DESeq2’s extensive functionality, a different package may be
appropriate for the analysis of RNA-seq data with unique experimental
designs, for example using linear-mixed effects models to perform
differential expression analysis of clustered datasets.

**Note:** You must change the below line, and all other lines loading
images, to the directory on your computer!!

<center>
![Overview](figures/overview.png)
</center>
Set the root directory for the whole markdown. THIS MUST BE SET TO THE
LOCATION OF THE FOLDER YOU DOWNLOADED!

```r
    knitr::opts_knit$set(root.dir = '/Users/shannon/Documents/GitHub/RNA-seq-Differential-Expression-workshop-June-2021/')
```

Lets start by loading the packages we will need:

```r
    library(dplyr)
    library(ggplot2)
    library(tximport)
    library(DESeq2)
    library(biomaRt)
    library(vsn)
    library(pheatmap)
    library(gplots)
    library(RColorBrewer)
    library(ComplexHeatmap)
    library(circlize)
    library(xtable)
    library(kableExtra)
```
------------------------------------------------------------------------

#### The dataset

The dataset that we are using comes from [this
paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625).
This data was collected from human airway smooth muscle cells to test
gene pathways effected by exposure to Glucocorticoid, which have been
historically used for their anti-inflammatory effects to treat
asthmatics. Four cell lines were treated with either a control vehicle
(untreated), **dexamethasone (dex)**, **albuterol (alb)**, or both
**dexamethasone and albuterol (co-treated)** for 18 hours before
transcriptomes were extracted.

### Read in raw count data

Now we can read in our data. How you read your data into DESeq2 depends
on what format your raw reads counts are in (individual files for each
sample, or a gene expression matrix) and how your read counts were
quantified (e.g. at the gene or transcript level). `DESeq2` provides a
specific function `DESeqDataSetFromHTSeqCount` to read in gene-level
read count abundances from *htseq-count*.

```r
    # read in the matrix we generated using htseq-count 
    cts <- as.matrix(read.table("Day-2/all_counts.txt", 
                                sep="\t", header = TRUE, row.names=1, 
                                stringsAsFactors = F))
    # quick look at the matrix 
    head(cts)
    tail(cts)

    # filter out these last 5 rows 
    cts <- cts[1:(nrow(cts)-5),]
    tail(cts)
```

If you estimated **transcript-level counts** (rather than **gene-level
counts** produced by htseq-count) using a method like *RSEM*, *Salmon*,
or *kallisto*, using the tximport() function from the [tximport
package](https://f1000research.com/articles/4-1521/v1). You may have
estimated transcript-level counts if you used a library-preparation
protocol that captures full length transcript information.

Even if you only plan to do a differential expression analysis at the
gene-level, it has been shown that [transcript-level estimates can
improve gene-level
inferences](https://f1000research.com/articles/4-1521/v1), therefore if
you are able to estimate counts at the transcript-level for your data,
it is beneficial to do so. Briefly, this method works by collapsing
transcript-level estimates into gene-level estimates, while an offset
matrix is calculated based on the average transcript length, that is
used in the differential expression analysis to correct for biases that
may be introduced by transcript-length differences between samples. You
can read more about how to do this in the [documnetation for
tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html).

If you collected **3’-end data**, e.g. with the **Lexogen QuantSeq
assay**, you should not do this correction for length, as there is no
length bias in your data. Doing this correction would introduce bias
into your data and likely distort your differential expression results.
For 3’-end data, it is best to read in the raw count matrix directly
using (`DESeqDataSetFromHTSeqCount`) or simply (`read.table()`).

------------------------------------------------------------------------

### Read in sample metadata

We also need to read in the sample annotation (metadata) that we
downloaded from the SRA, which contains sample labels, experimental
labels, and sequencing run information, etc.

```r
    # read in the file from the SRA metadata that has sample/experimental labels 
    colData <- read.csv("Day-2/sample_metadata.csv", row.names=1)
    head(colData)

    # order by SRA run accession 
    colData <- colData[order(colData$SRR),]
    # quick look 
    head(colData)
```

Lets have a look at our experimental design variable (drug treatment:)

```r
    # now make this a factor as it will be the variable we will use define groups for the differential expression analysis 
    colData$tx.group

    ##  [1] untreated Dex       Alb       Alb_Dex   untreated Dex       Alb      
    ##  [8] Alb_Dex   untreated Dex       Alb       Alb_Dex   untreated Dex      
    ## [15] Alb       Alb_Dex  
    ## Levels: Alb Alb_Dex Dex untreated
```

It is important that we make this variable a (`factor`) class variable,
with the reference group set as the variable we want to be considered
baseline expression. This has already been done for us, however if it
had not, you can create an ordered factor variable from a character
string in R using:

```r
    colData$tx.group <- factor(colData$tx.group, levels=c("untreated", "Dex", "Alb", "Alb_Dex"))
```
------------------------------------------------------------------------

### Construct the *DESeq2* data set & explore the characteristics of the data

*DESeq2* uses an object class called the `DESeqDataSet` that stores the
read counts, metadata, experimental design, and all the intermediate
values calculated during the analysis. `DESeqDataSet` extends the
`SummarizedExperiment` class object from the `SummarizedExperiment`
R/Bioconductor package that is commonly used to store data from
expression studies and other genomics assays in R.

Three elements are required to generate the `DESeqDataSet`:  
- matrix of raw counts  
- sample metadata (colData)  
- a design formula

Lets create the `DESeqDataSet` object.

```r
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = colData,
                                  design = ~ tx.group)
```

We could have also done this using the `DESeqDataSetFromHTSeqCount()`
function by specifying a `SampleTable` that includes the path to the
htseq-count files, however since we compiled the read counts into one
file, we can just load the dataset directly.

Before moving on, lets explore our DESeq2 class object a bit to get to
familar with its contents.

```r
    # have a quick look at the object 
    dds
    # print structure 
    str(dds)
    # several accessor functions exist to access specific data 'slots'
    head(counts(dds))
    head(colData(dds))
    # specific slots can also be accessed using the '@'
    dds@colData
```

Lets drop genes that have less than 10 reads across all samples, as
there just isn’t enough information for these genes to fit robust
statistical models to.

```r
    # drop genes with low counts 
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    dim(dds)

    ## [1] 24419    16
```

Lets also save the DESeq object at this point (so that we don’t have to
do the above everytime we want to work with our data).

```r
    save(dds, file = "DESeq2.rdata")
```

------------------------------------------------------------------------

### Normalization of raw counts

Before comparing expression levels of specific genes between samples,
and performing any differential expression analyses, it is important
that we normalize the data to account for variation in expression that
is not related to true differential expression. There are two major
sources variation that we need to adjust for in the normalization
process for RNA-seq data when we wish to compare expression levels
**between** samples:

#### Library size/sequencing depth

Although we generally try to pool samples together (each sample is
tagged with a barcode and samples are combined) at similar
concentrations in a sequencing run, some samples will end up being
sequenced more than others, leading to slight differences in how many
reads are produced for that sample, and therefore sequencing depth and
size. Furthermore, if samples are sequenced on separate runs, their
sequencing depths may be very different. If we don’t account for this
variation in sequencing depth, we might conclude some genes are
expressed at greater levels in a sample that has simply been sequenced
to a higher depth.

![](../figures/library_size.png)

#### Library composition

The presence of truly differentially expressed genes (in particular,
DEGs with very large fold changes) between samples will cause the number
of reads for other genes in those samples to be skewed. For example, in
the below example, gene C is differentially expressed between the two
samples, with much higher expression in sample 1. This high number of
reads causes fewer reads to be detected for other genes in this sample,
making it appear that these other genes are expressed at lower levels
than in sample 2, however this is simply an artifact of library
composition differences between the samples.

![](../figures/library_composition.png)

To correct for **library size** AND **library composition**, DESeq2 uses
a algorithm referred to as the **median-of-ratios** method. Although we
won’t go over how the algorithm works in detail, a brief summary of the
steps is:

1.  Take the log of all values in raw count matrix  
2.  Average each row (genes)
3.  Filter out genes with Infinity values
4.  Subtract average log count value from log of count for each cell
    (due to the laws of working with logarithms, this is essentially
    calculating the ratio of the counts for gene X in 1 sample to the
    average counts for gene X across all samples)
5.  Calculate the median of the ratios in each sample (column)
6.  Take exponents of medians to get the **size factors** for each
    sample/library.
7.  Divide the count for each gene in each sample by the size factor
    calculated for that sample.

This procedure will generate a matrix of read counts that are corrected
for both **library size** and **library composition**, and are stored in
our (`DESeqDataset`) object. DESeq2 uses the function
(`estimateSizeFactors()`) to perform this algorithm and calculate size
factors for each sample. Lets do this for our (`DESeqDataset`).

```r
    dds <- estimateSizeFactors(dds)
```

Note: [This video](https://www.youtube.com/watch?v=UFB993xufUU) from
StatQuest provides an excellent summary of the steps performed by
(`estimateSizeFactors()`) in order to calculate these size factors.

Once we have calculated the size factors, it can be helpful to look at
their distribution to get a feel for how they vary and how much
normalization between the samples is required.

```r
    sizeFactors(dds)

    hist(sizeFactors(dds), 
         breaks=6, col = "cornflowerblue",
         xlab="Size factors", ylab="No. of samples", 
         main= "Size factor distribution over samples")
```

After we have calculated the size factors, we can use the `counts()`
function, with `normalized` set to `TRUE`), to return the matrix of
counts where each column (each library/sample) have been divided by the
size factors calculated by the `estimateSizeFactors()` function.

```r
    counts_norm <- counts(dds, normalized=TRUE)
    head(counts_norm)
```
   
Comparing the normalized to the raw counts, we can clearly see that they
are different.

```r
    head(counts(dds, normalized=FALSE))
```
We can use this table of normalized read counts to compare values for
individual genes across samples. We might want to use this to (sanity)
check the expression of a few genes of interest, before we actually do
any statistical modelling. The abstract of the paper describes *DUSP1*,
a phosphatase with dual specificity for tyrosine and threonine, as a
well-known glucocorticoid-responsive gene.

```r
    # lets make a function to generate a quick plot of the normalized counts 
    gene_plot <- function(ENSG, gene_symbol){
      # save the normalized counts in a dataframe 
      cnts <- counts(dds, normalized=TRUE)
      colnames(cnts) <- colData(dds)$SRR
      # extract the counts for specified ENSG ID and add sample group data 
      df1 <- data.frame(log2(cnts[ENSG,]), colData(dds)$tx.group)
      colnames(df1) <- c(paste0("log2_gene"), "sample_group")
      # use ggplot2 to make a plot of counts vs sample group 
      p1<- ggplot(df1, aes(sample_group, log2_gene)) + 
        geom_jitter(aes(color = sample_group)) + 
        ggtitle(paste0(gene_symbol), " - Log2 Normalized counts")
      # print the plot 
      print(p1)
    }
    # now apply the function to print a plot for a specified gene 
    gene_plot(ENSG = "ENSG00000120129", gene_symbol = "DUSP1")
```

DUSP1 expression is consistently higher in the DEX samples than the
untreated, suggesting this gene is differentially expressed after DEX
treatment, validating prior knowledge and giving us confidence that our
experiment worked, sample labels are all correct, and we are well
positioned to make new discoveries with these data.

**Important note:** the normalized count matrix is normalized for
**library size and composition**, which means we can compare expression
levels of individual genes across samples. The read counts are NOT
normalized for gene length, so we cannot use this matrix to compare
expression levels between genes within the same sample. This is
important because some genes may simply pick up more reads than others
because they are larger, making them appear more highly expressed than a
smaller gene, which may not be the case.

For such comparisons between genes, we need to use measures such as:  
- *Transcripts per million (TPM)*  
- *Fragments per kilobase million (FPKM)*  
- *Reads per kilobase million (RPKM)*

[This
video](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)
provides an excellent explanation of *RPKM*, *FPKM*, & *TPM*, and
explains why it is better to use TPM if you need to correct for
**library size** AND **gene length**.

<center>
![](../figures/gene_length.png)
</center>

------------------------------------------------------------------------

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
