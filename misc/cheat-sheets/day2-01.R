# Day2 - 01 Exploratory analysis

setwd('~/Documents/GitHub/RNA-seq-Differential-Expression-workshop-June-2021/')

# read in the RDS object
dds <- readRDS("DESeq2.rds")

library(ggplot2)
library(DESeq2)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# drop genes with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# check the new dimensions
dim(dds)

rld <- rlog(dds, blind = FALSE)

# check first few rows
head(assay(rld))


# set plotting window to 1 row vs 2 columns
par(mfrow=c(1,2))

# plot standard log counts
cts <- counts(dds, normalized=FALSE)
plot(log2(cts[,1]+1), log2(cts[,2]+1), col = "cornflowerblue", xlab = "Sample 1", ylab = "Sample 2", main = "Log2 + 1")

# plot rlog counts
plot(assay(rld)[,1], assay(rld)[,2], col = "indianred", xlab = "Sample 1", ylab = "Sample 2", main = "rlog")


## look at variance distribution of genes
# calculate gene expression level variance between samples
var <- rev(rowVars(assay(rld))[order(rowVars(assay(rld)))])

# reset plotting window to 1 row and 1 column
par(mfrow=c(1,1))
# plot variance for genes accross samples
plot(var,
     las = 1,
     main="Sample gene expression variance",
     xlab = "Gene", ylab = "Variance")

# add vertical lines at specific gene number indexes
abline(v=1000, col="red")
abline(v=250, col="blue")
abline(v=500, col="green")


## check variance explained by PCs
# modify variable feature number to be used in PCA and hierachical clutering based on no. of most variable features
var_feature_n <- 500

# calculate the row variance
rv <- rowVars(assay(rld))

# order variance by size and select top 500 with most variance
select <- order(rv, decreasing = TRUE)[1:500]

# subset rlog values for genes with top varaince ranks
rld_sub <- assay(rld)[select, ]

# transpose the matrix (rows to columns and columns to rows)
rld_sub <- t(rld_sub)

# run principal components analysis
pca <- prcomp(rld_sub)

# extract the variance explained by each PC
percentVar <- pca$sdev^2/sum(pca$sdev^2)

# subset for first 5 elemets
percentVar <- percentVar[1:5]

# give the string names
names(percentVar) <- c("PC1", "PC2", "PC3", "PC4", "PC5")

# plot variance for top 10 PCs
barplot(percentVar, col = "indianred", las = 1, ylab = "% Variance", cex.lab = 1.2)

## make PCA plot

# construct data frame w/ PC loadings and add sample labels
pca_df <- as.data.frame(pca$x)

# add a column contaiing tx group
pca_df$tx.group <- dds@colData$tx.group

# add column containing sample IDs
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

# add sample names to data points
text((pca_df[, 2])~(pca_df[, 1]), labels = pca_df$tx.group, cex=0.6, font=2, pos=4)


## batch effect

# add variable to PCA dataframe for batch
pca_df$batch <- NA

# (artifically) select samples with a PC2 value below 10 to batch 1
pca_df$batch[pca_df$PC2 < 10] <- "Batch 1"
# (artifically) select samples with a PC2 value above 10 to batch 2
pca_df$batch[pca_df$PC2 > 10] <- "Batch 2"

# convert string to factor
pca_df$batch <- factor(pca_df$batch, levels = c("Batch 1", "Batch 2"))

# plot PC1 vs PC2 but only for batch 1 samples
plot(pca_df[pca_df$batch=="Batch 1", 1], pca_df[pca_df$batch=="Batch 1", 2],
     xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"),
     ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
     main=paste0("PC1 vs PC2 for ", var_feature_n, " most variable genes"),
     pch=16, cex=1.35, cex.lab=1.3, cex.axis = 1.15, las=1,
     panel.first = grid(),
     col=pca_df$col,
     ylim=c(-14,20))

# add points for batch 2 samples but use different shape to denote batch
points(pca_df[pca_df$batch=="Batch 2", 1], pca_df[pca_df$batch=="Batch 2", 2],
       col=pca_df$col,
       pch=2, cex=1.35)

# add legend
legend(9.5, 10.5, levels(pca_df$batch), pch = c(16, 2))
legend(1.5, 11.5, levels(pca_df$tx.group), pch = 16, col = pca_df$col)

# add sample names as text to points
text((pca_df[, 2])~(pca_df[, 1]), labels = pca_df$tx.group, cex=0.6, font=2, pos=4)

######### part 2 hierarchical clustering

# select top X no. of variable genes
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), var_feature_n)

# set up gene expression matrix
mat1 <- assay(rld)[topVarGenes,]

# set up colors for heatmap
col = colorRamp2(c(0, 9, 18), c("blue", "white", "red"))
cols1 <- brewer.pal(11, "Paired")

# set up annotation bar for samples
ha1 = HeatmapAnnotation(Group = colData(dds)$tx.group,
                        col = list(Group = c("untreated" = cols1[1],
                                             "Dex" = cols1[2],
                                             "Alb" = cols1[5],
                                             "Alb_Dex" = cols1[6])),
                        show_legend = TRUE)

# se up column annotation labels (samples)
ha = columnAnnotation(x = anno_text(colData(dds)$SRR,
                                    which="column", rot = 45,
                                    gp = gpar(fontsize = 10)))
# generate heatmap object
ht1 = Heatmap(mat1,
              name = "Expression",
              col = col,
              top_annotation = c(ha1),
              bottom_annotation = c(ha),
              show_row_names = FALSE,
              show_column_names = FALSE)

# plot the heatmap
draw(ht1, row_title = "Genes", column_title = "Top 500 most variable genes")


## reduced heatmap

# select sample groups to keep
ind_to_keep <- c(which(colData(rld)$group=="untreated"), which(colData(rld)$group=="Dex"))

# select top variable features
topVarGenes <- head(order(rowVars(assay(rld)[,ind_to_keep]), decreasing=TRUE), var_feature_n)

# set up gene expression matrix
mat1 <- assay(rld)[topVarGenes, ind_to_keep]

# set up colors for heatmap
col = colorRamp2(c(0, 9, 18), c("blue", "white", "red"))
cols1 <- brewer.pal(11, "Paired")

# subset coldata for samples in untx and ex groups
colData_sub <- colData(dds)[ind_to_keep, ]

# set up annotation bar for samples
ha1 = HeatmapAnnotation(Group = colData_sub$tx.group,
                        col = list(Group = c("untreated" = cols1[1],
                                             "Dex" = cols1[2])),
                        show_legend = TRUE)

# se up column annotation labels (samples)
ha = columnAnnotation(x = anno_text(colData_sub$SRR,
                                    which="column", rot = 45,
                                    gp = gpar(fontsize = 10)))
# generate heatmap object
ht1 = Heatmap(mat1, name = "Expression", col = col,
              top_annotation = c(ha1),
              bottom_annotation = c(ha),
              show_row_names = FALSE,
              show_column_names = FALSE)

# plot the heatmap
draw(ht1, row_title = "Genes", column_title = "Top 500 most variable genes")


### reduced heatmap with batch annotation

# which samples had values > 10 for PC2
pca_df$sample_ids[pca_df$PC2 > 10 & pca_df$tx.group=="untreated"]
# or <10 for PC2
pca_df$sample_ids[pca_df$PC2 > 10 & pca_df$tx.group=="Dex"]

# set the batch variable for these samples as batch 2
colData_sub$batch <- "Batch 1"
colData_sub$batch[colData_sub$SRR=="SRR1039516"] <- "Batch 2"
colData_sub$batch[colData_sub$SRR=="SRR1039517"] <- "Batch 2"

# add a new color palate for batch
cols2 <- brewer.pal(3, "Greens")

# set up annotation bar for samples
ha1 = HeatmapAnnotation(group = c(as.character(colData_sub$tx.group)),
                        batch = c(as.character(colData_sub$batch)),
                        col = list(group = c("untreated" = cols1[1], "Dex" = cols1[2]),
                                   batch = c("Batch 1" = cols2[2], "Batch 2" = cols2[3])),
                        show_legend = TRUE)

# generate heatmap object
ht2 = Heatmap(mat1, name = "Expression", col = col,
              top_annotation = c(ha1),
              bottom_annotation = c(ha),
              show_row_names = FALSE,
              show_column_names = FALSE)

# plot the heatmap
draw(ht2, row_title = "Genes", column_title = "Top 500 most variable genes")
