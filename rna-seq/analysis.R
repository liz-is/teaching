# Import & pre-process ----------------------------------------------------

# Import data from featureCounts
## Filename with output from featureCounts
countfile <- "scratch/counts.txt"
## Read in the data
countdata <- read.table(countfile, header=TRUE, row.names=1)
## Take a look at the first few lines
head(countdata)
colnames(countdata)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,-(1:5)] 
head(countdata)
colnames(countdata)

# Rename columns
## Manually
c("ctl1", "ctl2", "ctl3", "ctl4", "uvb1", "uvb2", "uvb3", "uvb4", "uvb5")
## Using paste and rep
?paste
paste("ctl", 1:4)
paste("ctl", 1:4, sep="")
c(paste("ctl", 1:4, sep=""), paste("uvb", 1:5, sep=""))
## Using gsub
?gsub
gsub(pattern=".fastq_tophat.accepted_hits.bam", replacement="", x=colnames(countdata))
colnames(countdata) <- gsub(pattern=".fastq_tophat.accepted_hits.bam", replacement="", x=colnames(countdata))
head(countdata)

# Convert to matrix
class(countdata)
countdata <- as.matrix(countdata)
class(countdata)
head(countdata)

# Create condition variable
condition <- gsub("\\d", "", colnames(countdata))
condition

# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet
coldata <- data.frame(row.names=colnames(countdata), condition)
coldata
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)
head(res)
table(res$padj<0.05)
# Order by adjusted p-value
res <- res[order(res$padj), ]
# Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
# Write results
write.csv(resdata, file="foff-ipsc-run1-diffexpr-results.csv")


# Data Visualization ------------------------------------------------------

# Plot dispersions
plotDispEsts(dds, main="Dispersion plot")

## Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))

plotPCA(rld)

# Sample distance heatmap
head(assay(rld))
assay(rld)[1:5,1:5]
t(assay(rld))[1:5,1:5]
dist(t(assay(rld)))
as.matrix(dist(t(assay(rld))))
sampleDists <- as.matrix(dist(t(assay(rld))))
heatmap(sampleDists)
## better heatmap with gplots
library(gplots)
heatmap.2(sampleDists)
heatmap.2(sampleDists, col=colorpanel(64, "steelblue", "white"), key=FALSE, trace="none")
heatmap.2(sampleDists, col=colorpanel(64, "black", "white"), key=FALSE, trace="none")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

# MA Plot
par(pch=16)
with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x"))
with(subset(res, padj<.05), points(baseMean, log2FoldChange, col="red", pch=16))
## add points
library(calibrate)
?textxy
res$Gene <- rownames(res)
with(subset(res, padj<.05), textxy(baseMean, log2FoldChange, labs=Gene, cex=1, col=2))

# Volcano plot
## Set point character
par(pch=16)
with(res, plot(log2FoldChange, -log10(pvalue), main="Volcano plot"))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), col="green"))
## Add legend
legend("topleft", legend=c("FDR<0.05", "|LFC|>1", "both"), pch=16, col=c("red","orange","green"))
## Label points
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=1))
