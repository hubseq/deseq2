source('/.Rprofile')
print("Inside R script!")
print(args)
library('DESeq2')
library('ggplot2')
library('dplyr')

args <- commandArgs(TRUE)
print(args[1])

# output_dir
odir <- args[4]
group1 <- args[6]
group2 <- args[8]

# get count and group data into R objects
countData <- read.csv(args[10], header = TRUE, sep = ",")
metaData <- read.csv(args[2], header = TRUE, sep = ",")
print(metaData)

# will create for loop later
print(paste('"^', group1, '$"', sep=''))
print(paste('"^', group2, '$"', sep=''))
metaData_group1_vs_2 <- filter(metaData, grepl(paste('^', group1, '$', sep=''), metaData$group) | grepl(paste('^', group2, '$', sep=''), metaData$group))
metaData_group_other <- filter(metaData, group != group1 & group != group2)
countData_group1_vs_2 <- select(countData, -contains(gsub("-",".",metaData_group_other$sample_id)))

print("metaData_group1_vs_2")
print(metaData_group1_vs_2)
print("metaData_group_other")
print(metaData_group_other)
print("countData_group1_vs_2")
head(countData_group1_vs_2)

# set up DESeq object and run DESeq
dds <- DESeqDataSetFromMatrix(countData=countData_group1_vs_2, 
                              colData=metaData_group1_vs_2, 
                              design=~group, tidy = TRUE)
dds <- DESeq(dds)

# create results table
res <- results(dds)
print("finished first DESeq")

# make gene_id a proper row
res <- cbind(Row.Names = rownames(res), res)
rownames(res) <- NULL
print(colnames(res))
colnames(res) <- c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")

# output results as CSV
write.csv(res, file=paste(odir, "deseq2.group1-vs-group2-results.csv", sep=''))

# volcano plots
jpeg(paste(odir,'deseq2.group1-vs-group2-volcano-plot.jpg', sep=''))
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05) - change to pvalue for now
with(subset(res, pvalue<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()
