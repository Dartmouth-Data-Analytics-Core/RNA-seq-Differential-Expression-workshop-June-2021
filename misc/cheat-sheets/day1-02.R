# Day1-02 Data Management and Setup

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

setwd('~/Documents/GitHub/RNA-seq-Differential-Expression-workshop-June-2021/')

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


# read in the file from the SRA metadata that has sample/experimental labels
colData <- read.csv("data/sample_metadata.csv", row.names=1)
head(colData)

# order by SRA run accession
colData <- colData[order(colData$SRR),]

# quick look
head(colData)

colData$tx.group
colData$tx.group <- factor(colData$tx.group, levels=c("untreated", "Dex", "Alb", "Alb_Dex"))
colData$tx.group


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ tx.group)

# have a quick look at the object
dds

# print structure
str(dds)

# several accessor functions exist to access specific data 'slots'
head(counts(dds))
head(colData(dds))

# specific slots can also be accessed using the '@'
dds@colData

saveRDS(dds, file = "DESeq2.rds")
