# Day 3 -02 annotation and visualization

library(vsn)
library(dplyr)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(readr)
library(circlize)
library(apeglm)
library(ggplot2)


# read in the flat file we downloaded and have a look at it
anno <- read.delim("data/GRCh38.p12_ensembl-97.txt", stringsAsFactors = T, header = T)
anno <- anno[order(anno$Chromosome.scaffold.name),]
dim(anno)


# have a look at the first few rows
head(anno)

tab1 <- table(anno$Chromosome.scaffold.name)
tab1[1:22]

any(duplicated(anno$Gene.stable.ID))


# use match() to find corresponding indicies (rows) for each ENSG ID
mat1 <- match(rownames(res_ord), anno$Gene.stable.ID)
table(is.na(mat1))

# add gene names to results as a new column
res_ord$gene <- as.character(anno$Gene.name[mat1])
head(res_ord, 20)

res_ord$chr <- as.character(anno$Chromosome.scaffold.name[mat1])
res_ord$start <- as.character(anno$Gene.start..bp.[mat1])
res_ord$end <- as.character(anno$Gene.end..bp.[mat1])
res_ord$strand <- as.character(anno$Strand[mat1])


### volcano plot
plot(res$log2FoldChange, -log10(res$pvalue),
     main = "Volcano plot",
     las = 1, col = "indianred",
     ylab = "- log10 P-value", xlab = "log2 Fold change")

# add horizontal lines to help guide interpretation
abline(h=-log10(0.05/nrow(res)), lty = 1, col = "black") # Bonferonni
abline(h=-log10(0.05), lty = 2, col = "black") # nominal P-value

## beautified volcano plot
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
  # if pvalue is NA don't assign a color - no plotting
  if(is.na(res_tmp$pvalue[i])){
    res_tmp$cols[i] <- NA
  }
  # if pvalue is < alpha AND LFC is > FC cutoff color red
  else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i] > fc_cutoff){
    res_tmp$cols[i] <- "indianred"
  }
  # if pvalue is < alpha AND LFC is < -FC cutoff color red
  else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i] < -fc_cutoff){
    res_tmp$cols[i] <- "indianred"
  }
  # if pvalue is < alpha AND LFC is not within cut off value color blue
  else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i]>-fc_cutoff & res_tmp$log2FoldChange[i]<fc_cutoff){
    res_tmp$cols[i] <- "cornflowerblue"
  }
  # if pvalue is > alpha AND LFC is > cut off value color gray
  else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] > fc_cutoff){
    res_tmp$cols[i] <- "gray47"
  }
  # if pvalue is > alpha and LFC is < -cut off value color gray
  else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] < -fc_cutoff){
    res_tmp$cols[i] <- "gray47"
  }
  # if pvalue is > alpha and LFC is not within cutoff values color light gray
  else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] < fc_cutoff){
    res_tmp$cols[i] <- "gray10"
  }
}

res_tmp$ENSG <- rownames(res_tmp)

# generate the plot
p = ggplot(res_tmp, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=col), alpha = 0.5, size =2.5, colour = res_tmp$cols, fill = res_tmp$cols)  +
  xlab("Log2 fold change") + ylab("-log10 Q-value") +
  ylim(0, 9) +
  xlim(-5, 11) +
  geom_hline(yintercept = -log10(alpha), color = "black", linetype = "dashed", size = 0.4) +
  # add vertical fold change lines
  geom_vline(xintercept = fc_cutoff, colour = "black", linetype="dotted") +
  geom_vline(xintercept = -fc_cutoff, colour = "black", linetype="dotted")
  theme(legend.key = element_blank()) +
  ggtitle("Control vs Dex")

# print the plot
print(p)


## add labels to interesing genes

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

# print the plot
print(p2)


## save results to csv file

# subset @ 5% adjusted pval sig. level
res_order_FDR_05 <- res_ord[res_ord$padj<0.05,]
nrow(res_order_FDR_05)

# write both to csv files
write.csv(as.data.frame(res_ord), file= "DE_results.csv")
write.csv(as.data.frame(res_order_FDR_05), file="DE_results.FDR.0.05.csv")

plotMA(res_ord, ylim=c(-6,6), main = "Raw Log2 Fold change")

# calculate shrunken fold change estimate
res_shrink <- lfcShrink(dds,
                        coef=paste0(resultsNames(dds)[which(resultsNames(dds)=="tx.group_Dex_vs_untreated")]),
                        type="apeglm")

par(mfrow=c(2,1))
plotMA(res_ord, ylim=c(-6,6), main = "Raw Log2 Fold change")
plotMA(res_shrink, ylim=c(-6,6), main = "Shrunken Log2 Fold change")

### hierarchical clustering

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
              show_row_names = FALSE,
              show_column_names = FALSE)

# plot the heatmap
par(mfrow=c(1,1))
draw(ht1, row_title = "Genes", column_title = "Hierachical clustering of DEGs (padj<0.05)")
