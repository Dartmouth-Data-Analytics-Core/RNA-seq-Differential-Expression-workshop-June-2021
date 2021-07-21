TO DO: 
- reduce text and try add more images to improve flow 
- remove/reduce content covered in statistical inference section (linear models, GLM principles etc.)


# Differential expression analysis in R
-------------------------------------

### Learning objectives:
-
------------------------------

### Introduction to DEG in R

After exploring our data, we are ready to run the differential
expression analysis. As our exploratory analysis showed the Alb treated
and co-treated samples did not cluster together, so going forward we
will focus on the comparison between untreated samples and Dex treated
samples. We will continue using *DESeq2* to perform the analysis.

**NOTE:** You must change the below line, and all other lines loading
images, to the directory on your computer!!

</center>
![Overview](/Users/shannon/Documents/GitHub/RNA-seq-Differential-Expression-workshop-June-2021/figures/overview.png)
</center>
Set the root directory for the whole markdown

```r
    knitr::opts_knit$set(root.dir = '/Users/shannon/Documents/GitHub/RNA-seq-Differential-Expression-workshop-June-2021/')
```

Lets start by loading the required libraries again.

```r
    library(tximport)
    library(DESeq2)
    library(biomaRt)
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

Read in the DESeq2 dataset we created in PART-1, which contains the raw
counts, normalization factors, and sample metadata.

```r
    load("Day-2/DESeq2.rdata")
```
------------------------------------------------------------------------

### Apply the DESeq2 procedure to the data

Now we apply the `DEseq()` function to the dataset to perform the
analysis. This function is the main powerhouse of the DESeq2 package and
does **a lot** under the hood. It is important you understand the
general principles of this analysis of before running your own analysis.

**At a high level, the major function performed by `DESeq2` are:**  
- Estimation of size factors (`estimateSizeFactors()`)  
- Estimation of dispersion (`estimateDispersions`)  
- Fitting of the negative binomial generalized linear model (GLM) and
wald statistics for differential expression testing (`nbinomWaldTest`)

Lets run `DESeq2` on our dataset:

```r
    # run the DEseq2 analysis
    dds <- DESeq(dds)
```

Before running the differential expression analysis, lets have a look at
some of the standard characteristics of RNA-seq data. The first and most
obvious thing to do is look at how the distribution of the raw counts.

```r
    hist(counts(dds, normalized=FALSE)[,5], breaks = 500, col="blue",
         xlab="Raw expression counts", ylab="Number of genes",
         main = "Count distribution for sample X")
```

Perhaps the most obvious feature of this distribution is the large
number of genes with very low count values. This occurs as there are
many genes expressed at low levels relative to the highly expressed
genes, which are fewer in number. This causes the distribution to have a
long right tail, ultimately caising the dynamic range of RNA-seq data to
be very large.

These features of how RNA-seq data is distributed are important in
selecting the statistical model used to test differential expression.
Importantly, we can see from the histogram that the data is **not**
normally distributed, therefore any statistical model based on the
normal distribution is not appropriate for this dataset. By looking
again at the matrix of raw counts, it is actually clear that RNA-seq is
integer count data, therefore we should use a statistical model for
count-based data.

```r
    head(counts(dds, normalized=FALSE))
```

At this point it might be useful to define a few terms that are really
important to know in order to understand as we fit statistical models to
RNA-seq data.

**mean **- the average count of a gene across samples **variance** - the
spread of count values across samples for a gene **dispersion** - the
amount that the variance deviates from the mean

One commonly used distribution for count data is **Poisson
distribution**, however, there is a feature of RNA-seq data that makes
the Poisson distribution a little to simplistic for such data, called
**overdispersion**.

**Overdispersion** describes the situation where the varaince for a set
of observations generally exceeds the mean of those observations. We can
visualize overdispersion in RNA-seq data by plotting the mean-variance
relationship for a group of replicates in our data.

```r
    # calculate mean and varaince for group of replicates
    mean_counts <- apply(counts(dds, normalized=FALSE)[,1:3], 1, mean)
    variance_counts <- apply(counts(dds, normalized=FALSE)[,1:3], 1, var)

    # plot the mean variance trend
    plot(log10(mean_counts), log10(variance_counts),
         ylim=c(0,9), xlim=c(0,9),
         ylab = "log10 (mean counts)", xlab = "log10 (varaince)",
         main = "Mean-variance trend", las = 1)

    # add line for x=y
    abline(0,1,lwd=2,col="red")
```

**We can clearly see a few features of the mean variance trend from this
plot:** 1. The data does not fall along the x = y line, as it would if
the mean = varaince. Instead, the varaince is generally greater than the
mean, making the varaince overdispersed. 2. There is more difference in
the varaince between low count genes than there is amongst higher count
genes, therefore the varaince is unequal across the range of count
values (non-constant variance is sometimes referred to as
**heteroscadicity**).

To account for this **overdispersion**, we use a generalization of the
*Poisson distribution* called the **negative-binomial (NB)
distribution**. The NB dist. includes a **dispersion parameter** that
accounts for the amount the variance exceeds the mean (the *Poisson
variance*). It is clearly important that we do this, because **the
varaince changes dramatically depending on the expression level of the
gene you are observing**.

We can plot a few different NB distributions to examine how the
dispersion parameter affects the spread of the data.

```r
    # generate a random varaible using the negative binomial distribution
    ### dispersion = 10
    par(mfrow=c(3,1))
    hist(rnbinom(n = 10000, mu = 100, size = 1/0.001),
         xlim = c(0, 300), xlab = "", breaks = 500,
         main = " Dispersion 0.001")
    ### dispersion = 10
    hist(rnbinom(n = 10000, mu = 100, size = 1/0.01),
         xlim = c(0, 300), xlab = "", breaks = 500,
         main = " Dispersion 0.01")
    ### dispersion = 10
    hist(rnbinom(n = 10000, mu = 100, size = 1/0.1),
         xlim = c(0, 300), xlab = "", breaks = 500,
         main = " Dispersion 0.1")
```

Note: The above example for plotting NB distributions at various
disperions was adapted from the *Data Analysis for the Life Sciences
series* available on edX and at
[rafalab](https://rafalab.github.io/pages/harvardx.html), and is an
excellent resource to learn more about how we model RNA-seq data for
differential expression analysis.

It is clear that as the disperion increases, the varaition around the
mean also increases. The mean, variance, and dispersion are linked by
the equation:

variance = mean + dispersion x 2 mean-squared ( var = mu + disp. \* mu^2
)

In order to accurately model differential expression for the genes in
our dataset, `DESeq2` uses this equation to obtain estimates for the
dispersion of each gene within each sample group (e.g. Control and Dex
separately).

**However,** for the small number of replicates avaiable in RNA-seq
data, these estimates of dispersion at the gene-level are often
inaccurate (yet another reason to use more replicates..).

To improve these gene-level estimates of dispersion, `DESeq2` uses
another statistical model called **empirical bayes** to *‘shrink’* these
inital dispersion estimates toward a *‘prior’* mean, which is calculated
by fitting a curve to the inital dispersion estimates.

This procedure produces **more accuarate estimates of disperion** as it
shares information across genes with similar expression levels to
predict a more approriate dispersion for those genes. This is rational
as the formula linking the mean, variance, and dispersion tells us that
the variance is the only thing affecting the magnitude of the dispersion
for genes with the similar mean expression.

**The major factors affecting how much a gene’s dispersion is shrunk
toward the prior mean are:**  
1. the number of samples in the group under consideration (use more
replicates!)  
2. how far the inital dispersion is from the prior mean

![DEseq2 dispersion estimation](../figures/dispersion_estimation.png)


This Figure taken frm the `DESeq2` paper demonstrates the process of
*shrinkage*, where the inital dispersion esimates for each gene
(estimated by maximum-likelihood) are shrunken towards the *prior mean*
(based on the fitted curve in red) to a final MAP estimate. For
dispersion estimates further away from the line, you can see that their
estimates are shrunken more than those are are originally closer to the
line.

**Why do we need to know about the dispersions:**  
This curve also shows the expected trend for dispersion estimates over a
range of expression levels. Importantly, the dispersion tends to
decrease as the mean increases. Through inspecting disperion estimates
for our own data, we can determine if the **NB model** is a good fit for
our data, and therefore if it can be used to accurately test DE.

**We can plot the dispersion estimates for our own data using:**

```r
    plotDispEsts(dds)
```

This is an example of a well calibrated set of dispersion estimates due
to these two features: the final MAP estimates are well scattered around
the fitted line, and the dispersion trend decreases with increasing mean
expression.

If the MAP estimates were more structured in these plots, we would be
concerned that the model is not estimating dispersions well for our
data, indicating something may be wrong with the dataset, e.g. outlier
samples, a batch effect, low quality samples/data, potential
contamination etc.

**It is important to confirm your dispersion estimates are well
calibrated before performing your differential expression analysis, as
accurate estimation of dispersion is critical in controlling the
false-positive rate in experiments with smaller sample sizes (most
RNA-seq experiments)**.

------------------------------------------------------------------------

### Differential expression analysis - Hypothesis testing

Now that we understand how the dispersions are estimated, we are ready
to fit the data and test each gene for differential expression!

We fit the data using a **generalized linear model (GLM)**. GLM’s are a
family of statistical models that generalize standard linear regression
in two ways:  
- use of probability distributions other than the normal distribution -
the use of a *link-function* that links the expression values in the
linear model to the experimental groups, in a way that these other
distributions (such as the NB) can be used.

Since we are need to model our counts using the negative-binomial
distribution, the GLM we will fit is of the NB family of GLMs.

**The DESeq2 model:**

![](../figures/neg-binom.png)

In order to fit the GLM, we need the **mean count of each gene** across
the samples in each experimental group, and the **dispersion of that
gene** in those groups. The mean count is a combination of the expected
expression level and the size factor, so that our model is corrected for
**library size and composition**.

The process of fitting the model to the expression and dispersion values
for each gene results in final set of **model coefficients** for each
sample group, which can be interpreted as the **log2 fold-change** in
expression for that gene between the baseline group and each comparison
group.

Each of the model coefficients has an associated **standard error**
associated with it, which we can use to calculate a **P-value** and
perform a process called **hypothesis testing**. Through hypothesis
testing we test the *null hypothesis* that the log2 fold-change between
experimental groups for an individual gene is not significnatly
different from 0 (no change in expression).

**The default test used by `DESeq2` for hypothesis testing is the
*Wald-test*, which is implemented as follows: **  
1. The *coefficient (log 2 fold-change)* is divided by the *standard
error* (measure of statistical accuracy of the measurement).  
2. The resulting *Z-statistic* is compared to a standard normal
distribution (mean = 0, sd = 1) in order to compute a P-value.  
3. If the P-value is less than our pre-determined threshold for
significance, we reject the null hypothesis and accept the alternative,
that the gene is significantly DE.

**Note:** `DESeq2` can also implement a *likelihood ratio test* (LRT),
which is used to compare expression accross more than two groups. For
example, if you collected samples over a range of time points and you
wanted to test if gene expression changed significantly over these time
points, you could use the LRT instead of the wald-test.

`DESeq2` already performed all of the steps for hypothesis testing using
the wald-test for us when we ran the `DESeq2()` function. All we have to
do is tell DESeq2 which results we want to look at, which can be done
using the `results()` function, and specifying the coefficients that we
want by using the `names` agument.

```r
    # quickly check the available coefficients we could extract
    resultsNames(dds)

    ## [1] "Intercept"                  "group_Dex_vs_untreated"    
    ## [3] "group_Alb_vs_untreated"     "group_Alb_Dex_vs_untreated"

    # get results for DEG analysis (and order by Pval) by specifying design
    res <- results(dds,
      name = "group_Dex_vs_untreated",
      alpha = 0.05,
      lfcThreshold = 0)
```

**A couple of things to note here:**

-   `alpha` is set to 0.05 (5%) to correct the P-values for multiple
    hypothesis testing (example with more detail on this coming up
    below). By default, the “BH” method is used (Benjamini & Hochberg)
    which controls the false discovery rate (FDR). Corrected P-values
    are found in the `padj` column of the `results()` output, while the
    uncorrected P-values are found in the `pvalue` column. Other methods
    to control for multiple hypothesis testing can be specified using
    the `pAdjustMethod` argument in the `results()` function, such as
    the more conservative **Bonferonni** method.

-   `lfcThreshold` is set to 0, and is the default value. This tests the
    hypothesis that the log2 fold change values between our experimental
    conditions are equal to 0. Different fold change values can be
    specified, which can be useful if you observe a large number of
    significantly differentially expressed genes with small fold
    changes, and you want to restrict the test to the genes with the
    largest differences (fold changes) between your conditions (we could
    also achieve this by restricting the results to genes with
    significant P-values AND have an absolute fold change &gt; a
    specific threshold, however when we do this, the P-values loose some
    of their meaning).

------------------------------------------------------------------------

#### Note on P-values:

P-value thresholds **do not need to be set at 0.05** for every
experiment. You can be more or less stringent than this dependningh on
the nature of your experiment: if you want to be very conservative and
restrict your results to few results that are likely to be true
positives, you may wish to restrict the results to a more stringent
threshold. If your experiment is very preliminary and you care less
about capturing some false positives than missing true positives, you
may wish to relax your threshold.

**Additional note:** To extract the results, we could also use the
`contrast` argument in a similar way to how we used the `names` agument.
The first group specified to `contrast` is used as the numerator in
calculating the fold change, and the second group is used as the
denominator, therefore the second group is used as the baseline for the
comparison.

```r
    res <- results(dds, alpha = 0.05,
      contrast = c("group", "Dex", "untreated"),
      lfcThreshold = 0)
```

This is useful when we have multiple levels in the experimental design
variable and we wish to extract coefficients for the results from
testing specific levels against one another. It is generally the same as
using the `names` argument to extract coefficients, with some exceptions
that are discussed in the DESeq2 documentation.

**Lets have a quick look at the results and how many genes were
statistically significant at an adjusted P-value threshold of 0.05. **

```r
    # order by adj Pval
    res_ord <- res[order(res$padj),]

    # quick check for how many DEGs with significance @ 5% level in either FC direction
    sum(res$padj < 0.05, na.rm=TRUE)

    sum(res$padj < 0.05 & res$log2FoldChange>2, na.rm=TRUE)

    sum(res$padj < 0.05 & res$log2FoldChange < -2, na.rm=TRUE)

```  

You may have noticed I am using `na.rm=TRUE` in the `sum()` function
above. Why might this be?

```r
    table(is.na(res$padj))
```

This is not a mistake, but rather part of a deliberate filtering process
conducted by `DESeq2`, in order to flag genes that have little or no
change of being differentially expressed.

This is of value as it means we can correct for fewer total tests and
increase our statistical power to identify true positives.The three ways
which `DESeq2` filters results are:  
- Genes with counts = 0 in all samples - Genes with extreme outliers
(determined using Cook’s distance) - *Independent filtering*
(identifying genes with low counts)

*Independent filtering*, DESeq2 carries out an iterative process where
it maximizes the value of the number of rejections over the quantiles of
the mean normalized counts. Once the maximum number of rejections is
identified, DESeq2 will select the quantile of the normalized counts
that is 1 standard deviation below this maximum, and filter any results
with mean counts below this threshold. It is essentially a fancy (and
cool) way of reducing the number of tests we need to run.

We can plot the number of rejections of the null hypotesis against mean
counts, along with a vertical line, to help us understand at which mean
count value DESeq2 chose to filter results for. Any genes with a mean
expression value below this line will have their `padj` values set to
NA, and discarded during multiple testing correction.

```r
    plot(metadata(res_ord)$filterNumRej,
         type="b", ylab="number of rejections",
         xlab="quantiles of filter (mean norm. counts)")
    lines(metadata(res_ord)$lo.fit, col="red")
    abline(v=metadata(res_ord)$filterTheta)
```

Its worth removing these results with NAs before moving forward to make
our lives a little easier when handling the adjusted P-values.

```r
    res_ord <- res_ord[!is.na(res_ord$padj),]
```
