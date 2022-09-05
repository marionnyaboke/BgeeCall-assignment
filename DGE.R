#!/usr/bin/Rscript

#PART ONE ----

#subset gene id and gene counts from gene_level_abundance+calls.tsv for FBdv:00007079 
{
calls_tsv <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX109278/gene_level_abundance+calls.tsv", header = TRUE)
head(calls_tsv)

calls_tsv2 <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX109279/gene_level_abundance+calls.tsv", header = TRUE)
head(calls_tsv2)

subset1 <- calls_tsv[ , c("id", "counts")] # all rows, two columns
head(subset1)

subset2 <- calls_tsv2[ , c("id", "counts")] 
}

#subset gene id and gene counts from gene_level_abundance+calls.tsv for UBERON:0000066 
{
  calls_tsv3 <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX493950/gene_level_abundance+calls.tsv", header = TRUE)
  head(calls_tsv)
  
  calls_tsv4 <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX493999/gene_level_abundance+calls.tsv", header = TRUE)
  head(calls_tsv2)
  
  calls_tsv5 <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX1720957/gene_level_abundance+calls.tsv", header = TRUE)
  head(calls_tsv)
  
  calls_tsv6 <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX1720958/gene_level_abundance+calls.tsv", header = TRUE)
  head(calls_tsv2)
  
  subset3 <- calls_tsv3[ , c("id", "counts")] # all rows, two columns

  subset4 <- calls_tsv4[ , c("id", "counts")] 
  
  subset5 <- calls_tsv5[ , c("id", "counts")] # all rows, two columns

  subset6 <- calls_tsv6[ , c("id", "counts")] 
}

# writing a txt file
{
write_delim(subset1, file = "/home/nyamarim/Marion/SRX493950-counts.txt")

write_delim(subset2, file = "/home/nyamarim/Marion/SRX493999-counts.txt")

write_delim(subset3, file = "/home/nyamarim/Marion/SRX493950-counts.txt")

write_delim(subset4, file = "/home/nyamarim/Marion/SRX493999-counts.txt")

write_delim(subset5, file = "/home/nyamarim/Marion/SRX1720957-counts.txt")

write_delim(subset6, file = "/home/nyamarim/Marion/SRX1720958-counts.txt")
}

# Take 'all' counts and combine into one big dataframe

# required packages
library(tibble)

# # where are we?
cntdir <- here::here("/home/nyamarim/Marion/")
pat <- "-counts.txt"
file.all <- list.files(path = cntdir,
                         pattern = pat,
                         all.files = TRUE,
                         recursive = FALSE,
                         ignore.case = FALSE,
                         include.dirs = FALSE)

# we choose the 'all' series
myfiles <- file.all
DT <- list()

# read each file as array element of DT and rename the last 2 cols
# we created a list of single sample tables
for (i in 1:length(myfiles) ) {
  infile = paste(cntdir, myfiles[i], sep = "/")
  DT[[myfiles[i]]] <- read.table(infile, header = T, stringsAsFactors = FALSE)
  cnts <- gsub("(.*)-counts.txt", "\\1", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on first ID columns
reads_count <- DT[[myfiles[1]]]

# inspect
head(reads_count)

# we now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(reads_count, y, by = c("ID"))
  reads_count <- z
}

# ID column becomes rownames
rownames(reads_count) <- reads_count$ID
reads_count <- reads_count[,-1]

# add total counts per sample
reads_count <- rbind(reads_count,
                     tot.counts=colSums(reads_count))

# transpose table for readability

reads_count_summary <- reads_count[grep("^__", rownames(reads_count), 
                                        perl=TRUE, invert=FALSE), ]


# transpose table
t(reads_count_summary)

# write summary to file
#write.csv(reads_count_summary,
#          file = here::here("/home/nyamarim/Marion/"reads_count_summary.csv"),
#          row.names = TRUE)

# take all data rows to a new table
reads_count <- reads_count[grep("__", rownames(reads_count), perl=TRUE, invert=TRUE), ]

# inspect final merged table
head(reads_count, 3)

# colnames(reads_count) <- gsub("(_.*$)", "", colnames(reads_count))

# reads_count <- rownames_to_column(reads_count,"transcript_id")
write.csv(reads_count, 
          file = here::here("/home/nyamarim/Marion/read_counts.csv"))

# cleanup intermediate objects
rm(y, z, i, DT, reads_count)

#PART TWO ----

# Load packages for differential gene expression analysis
library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)
library(ggplot2) #Best plots
library(ggrepel)
library(dplyr)

# Install missing package
#BiocManager::install("HTSFilter")

#library(HTSFilter)

# set the directory from which files are imported
directory <- "~/Marion/"
dir(directory)

# create metadata file
sample_names <- c('SRX109278', 'SRX109279', 'SRX493950', 'SRX493999', 'SRX1720957', 'SRX1720958')
devstage <- c("FBdv:00007079", "FBdv:00007079", "UBERON:0000066", "UBERON:0000066", "UBERON:0000066", "UBERON:0000066")
data <- data.frame(sample_names, devstage)

write.csv(data, 
          file = here::here("/home/nyamarim/Marion/Metadata.csv"), row.names = FALSE)


# read files
rawCountTable <- read.csv(paste0(directory,"read_counts.csv"), header=TRUE,
                          row.names=1)
sampleInfo <- read.csv(paste0(directory,"Metadata.csv"), header=TRUE,
                       row.names=1)

# Inspect
head(rawCountTable)
nrow(rawCountTable)

rawCountTable <- rawCountTable[,match(rownames(sampleInfo),
                                      colnames(rawCountTable))]
colnames(rawCountTable)
rownames(sampleInfo)

# create the ‘condition’ column
condition <- c("FBdv:00007079", "FBdv:00007079", "UBERON:0000066", "UBERON:0000066", "UBERON:0000066", "UBERON:0000066")


sampleInfo = data.frame(sampleInfo, condition)


# create a DGEList data object
dgeFull <- DGEList(rawCountTable, group=sampleInfo$condition.2)


# add the sample information object in the DGEList data
dgeFull$sampleInfo <- sampleInfo

# check the number of genes with no expression in all samples
table(rowSums(dgeFull$counts==0)==15)

FALSE 
14297 

# Filtering non-expressed and lowly-expressed genes.
keep.exprs <- filterByExpr(dgeFull, group=sampleInfo$condition)
filtered.counts <- dgeFull[keep.exprs,, keep.lib.sizes=FALSE]

# preparing the data object for the analysis 
# select the subset paired-end samples from dgeFull
dge <- DGEList(filtered.counts$counts)
dge$sampleInfo <- dgeFull$sampleInfo[dgeFull$sampleInfo$type=="paired-end",]

table(rowSums(filtered.counts$counts==0)==15)
FALSE 
10613 

# preparing the data object for the analysis 
# select the subset paired-end samples from degFull
#dge <- DGEList(dgeFull$counts)
#dge$sampleInfo <- dgeFull$sampleInfo[dgeFull$sampleInfo$type=="paired-end",]

# data exploration and quality assessment
# extract pseudo-counts (ie \(\log_2(K+1)\))
pseudoCounts <- log2(dge$counts+1)
head(pseudoCounts)

# histogram for pseudo-counts
hist(pseudoCounts[,"SRR10729843"])

# boxplot for pseudo-counts
boxplot(pseudoCounts, col="cyan")

# MDS for pseudo-counts (using limma package)
plotMDS(pseudoCounts)

# PCA plot
plotMDS(pseudoCounts, gene.selection="common")

# Differential expression analysis
# remove genes with zero counts for all samples

dge <- DGEList(dge$counts[apply(dge$counts, 1, sum) != 0, ],
               group=sampleInfo$condition)
dge$sampleInfo <- dge$sampleInfo

# estimate the normalization factors
dge <- calcNormFactors(dge, method="TMM")

# estimate common and tagwise dispersion
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# perform an exact test for the difference in expression between the conditions
dgeTest <- exactTest(dge)

# Independant filtering
#  remove low expressed genes

filtData <- HTSFilter(dge)$filteredData

dgeTestFilt <- exactTest(filtData)

# Diagnostic plot for multiple testing
# plot a histogram of unadjusted p-values
hist(dgeTest$table[,"PValue"], breaks=50)

# plot a histogram of unadjusted p-values after filtering
hist(dgeTestFilt$table[,"PValue"], breaks=50)

# Inspecting the results
# extract a summary of the differential expression statistics
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
head(resNoFilt)

resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))
head(resFilt)

# compare the number of differentially expressed genes with and without filtering
# before independent filtering
sum(resNoFilt$table$FDR < 0.01)

# after independent filtering
sum(resFilt$table$FDR < 0.01)

# extract and sort differentially expressed genes
sigDownReg <- resFilt$table[resFilt$table$FDR<0.01,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
head(sigDownReg)

sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
head(sigUpReg)

# write the results in csv files
write.csv(sigDownReg, file="/home/nyamarim/Marion/sigDownReg-filtered.csv")
write.csv(sigUpReg, file="/home/nyamarim/Marion/sigUpReg-filtered.csv")

# Interpreting the DE analysis results
# create a MA plot with 10% differentially expressed genes

plotSmear(dgeTestFilt,
          de.tags = rownames(resFilt$table)[which(resFilt$table$FDR<0.01)])

# create a Volcano plot

# create a column with thresholds
significance = -log10(resFilt$table$FDR) < 1
resFilt$table$significance <- significance
dim(resFilt$table)

## setting the values
resFilt$table$diffexpressed <- "No"

# if logFC > 1 and FDR < 0.01, set as "UP"
resFilt$table$diffexpressed[resFilt$table$logFC > 1 & resFilt$table$FDR < 0.01] <- "Up"

# if log2Foldchange < -1 and pvalue < 0.01, set as "DOWN"
resFilt$table$diffexpressed[resFilt$table$logFC < -1 & resFilt$table$FDR < 0.01] <- "Down"

# set different colors
mycolors <- c("cornflowerblue", "black", "firebrick")
names(mycolors) <- c("Down", "No", "Up")

dev.off()


input <- cbind(gene=rownames(resFilt$table), resFilt$table) #convert the rownames to a column
volc = ggplot(input, aes(logFC, -log10(PValue))) + #volcanoplot with logFC versus pvalue
  geom_point() +
  theme_classic()+
  geom_vline(xintercept=c(-1.0, 0.8), col="orange", linetype="dashed") +
  geom_hline(yintercept=-log10(0.01), col="orange" , linetype="dashed")+
  geom_point(aes(col=diffexpressed)) + #add points colored by significance
  scale_color_manual(values= mycolors) + 
  ggtitle("FBdv:00007079 vs UBERON:0000066")
volc <- volc+geom_text_repel(data=head(input, 500), aes(label=gene), size = 3) #adding text for the top 20 genes
volc
#traceback()

# transform the normalized counts in log-counts-per-million
y <- cpm(dge, log=TRUE, prior.count = 1)

# select 1% differentially expressed genes and produce a heatmap
selY <- y[rownames(resFilt$table)[resFilt$table$FDR<0.01 & 
                                    abs(resFilt$table$logFC)>1.5],]

cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]

dev.off() 

cim(selY, margins = c(10, 10))

