## Day3-01 Differential expression analysis in R

library(DESeq2)
setwd('~/Documents/GitHub/RNA-seq-Differential-Expression-workshop-June-2021/')

dds <- readRDS("DESeq2.rds")

# run the DEseq2 analysis
dds <- DESeq(dds)


hist(counts(dds, normalized=FALSE)[,5],
     breaks = 500, col="blue",
     xlab="Raw expression counts",
     ylab="Number of genes",
     main = "Count distribution for sample X")


head(counts(dds, normalized=FALSE))

## plot mean v var

# calculate mean and variance for group of replicates
mean_counts <- apply(counts(dds, normalized=FALSE)[,1:3], 1, mean)
variance_counts <- apply(counts(dds, normalized=FALSE)[,1:3], 1, var)

# plot the mean variance trend
plot(log10(mean_counts), log10(variance_counts),
     ylim=c(0,9), xlim=c(0,9),
     ylab = "log10 (mean counts)", xlab = "log10 (variance)",
     main = "Mean-variance trend", las = 1)

# add line for x=y
abline(0,1,lwd=2,col="red")


## neg binomial distributions

# set the plotting window to 3 rows and 1 column
par(mfrow=c(3,1))

### dispersion = 0.001
hist(rnbinom(n = 10000, mu = 100, size = 1/0.001),
     xlim = c(0, 300), xlab = "", breaks = 500,
     main = " Dispersion 0.001")

### dispersion = 0.01
hist(rnbinom(n = 10000, mu = 100, size = 1/0.01),
     xlim = c(0, 300), xlab = "", breaks = 500,
     main = " Dispersion 0.01")

### dispersion = 0.1
hist(rnbinom(n = 10000, mu = 100, size = 1/0.1),
     xlim = c(0, 300), xlab = "", breaks = 500,
     main = " Dispersion 0.1")


## plot disp est for dds
par(mfrow=c(1,1))
plotDispEsts(dds)

## check results of GLM - wald test

# quickly check the available coefficients we could extract
resultsNames(dds)

# get results for DEG analysis (and order by Pval) by specifying design
res <- results(dds,
               name = "tx.group_Dex_vs_untreated",
               alpha = 0.05)

# check dminesions on tables
dim(head)

# print top of results table
head(res)


# order by adj Pval
res_ord <- res[order(res$padj),]

# print top of results table
head(res_ord)


# quick check for how many DEGs with significance @ 5% level in either FC direction
sum(res$padj < 0.05, na.rm=TRUE)
sum(res$padj < 0.05 & res$log2FoldChange > 2, na.rm=TRUE)
sum(res$padj < 0.05 & res$log2FoldChange < -2, na.rm=TRUE)


res <- results(dds,
               name = "tx.group_Dex_vs_untreated",
               alpha = 0.05,
               pAdjustMethod = "bonferroni")

# how many significant DEGs
sum(res$padj < 0.05, na.rm=TRUE)

# count na values
table(is.na(res$padj))

#threshold 1
0.05/10000

# threshold 2
0.05/20000

# independent filtering
# look at independent filtering results table
metadata(res_ord)$filterNumRej

# plot number of rejections at specific quantiles
plot(metadata(res_ord)$filterNumRej,
     type="b", ylab="number of rejections",
     xlab="quantiles of filter (mean norm. counts)")

# add line connecting points together
lines(metadata(res_ord)$lo.fit, col="red")

# add vertical line at selected independent filtering threshold
abline(v=metadata(res_ord)$filterTheta)

res_ord <- res_ord[!is.na(res_ord$padj),]

write.csv(as.data.frame(res_ord), file = "DE_results.csv")

# save RDS object
saveRDS(dds, file = "DESeq2.rds")
