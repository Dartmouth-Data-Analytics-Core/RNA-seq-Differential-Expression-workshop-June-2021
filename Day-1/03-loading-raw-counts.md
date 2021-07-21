
## TO DO

incorproate section and implementations of norm approaches from rnaseq I workshop


# Count data & normalization for RNA-seq
------------------------------

### Learning objectives:
- Learn how to read raw count data from an RNA-seq experiment into R
- Understand why read counts must be normalized in RNA-seq data
- Learn the principles behind the major normalization strategies in RNA-seq and when to apply them
- Learn how to perform these normalization strategies

------------------------------

### Importing count data into R  

Several popular R-packges designed for exploration and statistical
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


![](../figures/overview.png)



```r
setwd('~/Documents/GitHub/RNA-seq-Differential-Expression-workshop-June-2021/')
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
cts <- as.matrix(read.table("data/all_counts.txt",
                            sep="\t",
                            header = TRUE,
                            row.names=1,
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
