# Day1-03 Normalization

# read in the RDS object
dds <- readRDS("DESeq2.rds")

cts <- counts(dds, normalized=FALSE)

library(ggplot2)
library(tximport)
library(DESeq2)
library(biomaRt)
library(vsn)
library(pheatmap)


## CPM
# look at the counts object
head(cts)

# write a function that will calculate TPM
cpm <- function(counts) {
  cpm <- c()
  for(i in 1:length(counts)){
    cpm[i] <- counts[i] / sum(counts) * 1e6
  }
  cpm
}

# apply function to the columns of raw counts data
# we start at the third column because the first two columns have the ensemble IDs and gene names
cts_cpm <- apply(cts[, 3:5], 2, cpm)
## NOTE: we are calculating tpm for first 3 samples only to save time..
# add gene info columns back in
cts_cpm <- cbind(cts[, c(1,2)], cts_cpm)

# write to file
write.csv(cts_cpm, file="cts_CPM.csv")


## TPM
# read in gene lengths matrix (pre made for you)
gene_lengths <- read.table("data/gene-lengths-grch38.tsv", sep="\t", stringsAsFactors=FALSE, header=TRUE)

# look at the lengths object
head(gene_lengths)

# write a function that will calculate TPM
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  tpm <- c()
  for(i in 1:length(counts)){
    tpm[i] <- rate[i] / sum(rate) * 1e6
  }
  tpm
}

# apply function to the columns of raw counts data
cts_tpm <- apply(cts[, 3:5], 2, tpm, gene_lengths$length)
## NOTE: we are calculating tpm for first 3 samples only to save time..

# add gene info columns back in
cts_tpm <- cbind(cts[, c(1,2)], cts_tpm)

# write to file
write.csv(cts_tpm, file="cts_TPM.csv")


# read in file containing all TPM counts (pre-made for you)
cts_tpm_full <- read.csv("data/all_counts_TPM-full.csv")

# get expression values for DUSP1 row
DUSP1_tpm <- cts_tpm_full[cts_tpm_full$gene_name=="DUSP1",]

# remove gene info columns
DUSP1_tpm <- DUSP1_tpm[ ,c(4:ncol(DUSP1_tpm))]

# convert to a numeric vector
DUSP1 <- as.numeric(DUSP1_tpm[1,])

# generate barplot of gene expression across samples
ppi=300
png("DUSP1_tpm.png")
barplot(DUSP1,
        col="lightblue", ylab="TPM", xlab="sample",
        main = "DUSP1 expression", las = 1)
dev.off()


## FPKM

# write a function that will calculate TPM
fpkm <- function(counts, lengths) {
  rate <- counts / lengths
  fpkm <- c()
  for(i in 1:length(counts)){
    fpkm[i] <- rate[i] / sum(counts) * 1e9
  }
  fpkm
}

# apply function to the columns of raw counts data
cts_fpkm <- apply(cts[, 3:5], 2, fpkm, gene_lengths$length)
## NOTE: we are calculating fpkm for first 3 samples only to save time..

# add gene info columns back in
cts_fpkm <- cbind(cts[, c(1,2)], cts_fpkm)

# write to file
write.csv(cts_fpkm, file="cts_FPKM.csv")

## DESeq2 normalization

dds <- estimateSizeFactors(dds)
# extract size factors
sizeFactors(dds)

# plot histogram
hist(sizeFactors(dds),
     breaks=6, col = "cornflowerblue",
     xlab="Size factors", ylab="No. of samples",
     main= "Size factor distribution over samples")

# calculate normalized counts
counts_norm <- counts(dds, normalized=TRUE)

# print top rows
head(counts_norm)
head(counts(dds, normalized=FALSE))


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
