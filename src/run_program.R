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

# get count and group data into R objects
countData <- read.csv(args[6], header = TRUE, sep = ",")
metaData <- read.csv(args[2], header = TRUE, sep = ",")
metaData

# will create for loop later
metaData_group1_vs_2 <- filter(metaData, grepl("^group1$", metaData$group) | grepl("^group2$", metaData$group))
metaData_group3_vs_4 <- filter(metaData, group == "group3" | group == "group4")
countData_group1_vs_2 <- select(countData, -contains(gsub("-",".",metaData_group3_vs_4$sample_id)))
countData_group3_vs_4 <- select(countData, -(starts_with(gsub("-",".",metaData_group1_vs_2$sample_id)) & ends_with(gsub("-",".",metaData_group1_vs_2$sample_id))))

metaData_group1_vs_2
metaData_group3_vs_4
head(countData)
head(countData_group1_vs_2)
head(countData_group3_vs_4)

# set up DESeq object and run DESeq
dds <- DESeqDataSetFromMatrix(countData=countData_group1_vs_2, 
                              colData=metaData_group1_vs_2, 
                              design=~group, tidy = TRUE)
dds <- DESeq(dds)

# create results table
res <- results(dds)
print("finished first DESeq")

# output results as CSV
write.csv(res, file=paste(odir, "deseq2.group1-vs-group2-results.csv", sep=''))

dds2 <- DESeqDataSetFromMatrix(countData=countData_group3_vs_4,
                              colData=metaData_group3_vs_4,
                              design=~group, tidy = TRUE)
dds2 <- DESeq(dds2)
res2 <- results(dds2)
write.csv(res2, file=paste(odir, "deseq2.group3-vs-group4-results.csv", sep=''))
