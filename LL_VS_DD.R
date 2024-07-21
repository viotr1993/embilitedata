#Read featureCount files and assign conditions
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("WB", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E6", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "wb_e6_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "wb_e6_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "wb_e6_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"wb_e6_dd_ll.csv")
###################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("WB", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E10", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "wb_e10_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "wb_e10_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "wb_e10_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"wb_e10_dd_ll.csv")
##################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("WB", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E17", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "wb_e17_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "wb_e17_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "wb_e17_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"wb_e17_dd_ll.csv")
##################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("WB", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E20", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "wb_e20_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "wb_e20_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "wb_e20_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"wb_e20_dd_ll.csv")
#############################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("CER", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E10", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "CER_e10_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "CER_e10_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "CER_e10_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"CER_e10_dd_ll.csv")
##################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("CER", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E17", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "CER_e17_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "CER_e17_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "CER_e17_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"CER_e17_dd_ll.csv")
##################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("CER", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E20", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "CER_e20_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "CER_e20_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "CER_e20_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"CER_e20_dd_ll.csv")
#####################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("HYP", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E10", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "HYP_e10_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "HYP_e10_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "HYP_e10_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"HYP_e10_dd_ll.csv")
##################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("HYP", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E17", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "HYP_e17_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "HYP_e17_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "HYP_e17_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"HYP_e17_dd_ll.csv")
##################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("HYP", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E20", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "HYP_e20_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "HYP_e20_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "HYP_e20_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"HYP_e20_dd_ll.csv")
#################################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("PI", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E17", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "PI_e17_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "PI_e17_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PI_e17_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"PI_e17_dd_ll.csv")
##################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("PI", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E20", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "PI_e20_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "PI_e20_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PI_e20_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"PI_e20_dd_ll.csv")
####################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("EYE", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E6", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "eye_e6_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "eye_e6_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "eye_e6_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"eye_e6_dd_ll.csv")
#############
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("EYE", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E10", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "eye_e10_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "eye_e10_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "eye_e10_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"eye_e10_dd_ll.csv")
####################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("EYE", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E17", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "eye_e17_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "eye_e17_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "eye_e17_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"eye_e17_dd_ll.csv")
#########################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("EYE", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("E20", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "eye_e20_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "eye_e20_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "eye_e20_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"eye_e20_dd_ll.csv")
######################
##################
countdata<-read.csv('broiler_count_3_noutliers.csv', header= T, row.names=1)
countdata[1:5]<- NULL
countdata_WHE20 <- countdata[,grepl("PI", colnames(countdata))]
countdata_WHE20 <- countdata_WHE20[,grepl("LL|DD", colnames(countdata_WHE20))]
countdata_WHE20 <- countdata_WHE20[,grepl("ADULT", colnames(countdata_WHE20))]

#set sample conditions
#set sample conditions
sample_id<-as.data.frame(names(countdata_WHE20))
library(stringr)
{
  sample_id$ST<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,1]
  sample_id$TS<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,2]
  sample_id$Treat<-str_split(sample_id$`names(countdata_WHE20)`,'_',simplify = T)[,3]
  sample_id$Treat<-gsub('1','',sample_id$Treat)
  sample_id$Treat<-gsub('2','',sample_id$Treat)
  sample_id$Treat<-gsub('3','',sample_id$Treat)
  sample_id$Treat<-gsub('4','',sample_id$Treat)
  sample_id$Treat<-gsub('A','',sample_id$Treat)
}
condition<-sample_id$Treat
stage<-sample_id$ST
tissue<-sample_id$TS


# Analysis with DESeq2, first based on tissue type 

library(DESeq2)
library (ggplot2)

#colate the data into a dataframe
coldata <- data.frame(row.names=colnames(countdata_WHE20),stage,condition)
#start DGE 

dds <- DESeqDataSetFromMatrix(countData = countdata_WHE20,
                              colData = coldata,
                              design= ~condition)
##
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
vsd <- vst(dds)
mat <- assay(vsd)
assay(vsd)<-mat

#Plot PCA, but with shapes and forms so it is clearer

data <- plotPCA(vsd, intgroup = c( "stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=stage, )) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Tissue and incubation condition")
###############
######Testing for light vs dark 
#start DGE 
##

outputPrefix <- "PI_adult_dd_ll"
dds <- estimateSizeFactors(dds)
gene_counts<-counts(dds,normalized=T)
#write.csv(gene_counts,file="gene_counts_chicken.csv")
idx <- rowSums( counts(dds, normalized=F) >= 5 ) >= 3
dds <- dds[idx,]
# Run the DESeq pipeline
# shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="tissue_Fissure_vs_Dorsal", type="apeglm")
#find and order genes with signifiant p value
dds <- DESeq(dds)
resultsNames(dds)
res_DDvsLD <- results(dds, contrast=c("condition","LL","DD"), alpha = 0.1)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
#shrink
res_DDvsLD <- lfcShrink(dds, coef="condition_LL_vs_DD", type="apeglm", res=res_DDvsLD)
#plots
plotMA(res_DDvsLD, ylim=c(-2,2), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)
hist(res_DDvsLD$pvalue, col = "blue", 
     border = "white", main = "After skrinkage", xlab = "p-values")
#subset only p value lesss or equal 0.05
res_DDvsLD<- subset(res_DDvsLD, padj<0.05)
res_DDvsLD <-res_DDvsLD[order(res_DDvsLD$padj),]
# save data results and normalized reads to csv
resdata_DDvsLD <- merge(as.data.frame(res_DDvsLD), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata_DDvsLD)[1] <- 'gene'
write.csv(resdata_DDvsLD, file = paste0(outputPrefix, "-results-with-normalized.csv"))
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')
# produce DataFrame of results of statistical tests
mcols(res_DDvsLD, use.names = T)
write.csv(as.data.frame(mcols(res_DDvsLD, use.name = T)),file = paste0(outputPrefix, "-test-ddld.csv"))
# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean_DDvsLD <- replaceOutliersWithTrimmedMean(dds)
ddsClean_DDvsLD <- DESeq(ddsClean_DDvsLD)
tab_DDvsLD <- table(initial = results(dds)$padj <= 0.05,
                    cleaned = results(ddsClean_DDvsLD)$padj <= 0.05)
addmargins(tab_DDvsLD)
write.csv(as.data.frame(tab_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean_DDvsLD = subset(res_DDvsLD, padj<0.05)
res_DDvsLD <- resClean_DDvsLD[order(resClean_DDvsLD$padj),]
write.csv(as.data.frame(resClean_DDvsLD),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Convert final results .csv file into .txt file
results_csv <- "PI_adult_dd_ll-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "PI_adult_dd_ll-replaceoutliers-results.txt"
# IMPORT GTF
gtf <- rtracklayer::import('Gallus_gallus_gca000002315v5.GRCg6a.109.gtf')
gtf_df=as.data.frame(gtf)
gtf_df[1:9]<-NULL
gtf_df[2]<-NULL
gtf_df[3]<-NULL
gtf_df[4:15]<-NULL
library(tidyverse)
gtf_df<-gtf_df %>% distinct(gene_id, .keep_all = TRUE)

a <- read.table(results_txt, head=TRUE)

names(a) <- c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj")
m <- merge(a, gtf_df, by="gene_id")
write.csv(m,"PI_adult_dd_ll.csv")
################
