#Misbah Sehar 16/03/2023
#Updated: 11/07/2023
#Updated: 03/06/2024

#Aim to carry out differential expression analysis comparing well-watered and antitranspirant-treated plant anther samples with unsprayed samples.

library("DESeq2")
library("dplyr")
library("tidyr")
library(ggplot2)
library(ggrepel)


setwd("~/Desktop/RNASeq_data_analysis_files/kallisto_results")
abundances <- read.table("Anther_lines_tpm.tsv")
head(abundances)

#Remove low confidence (LC) genes - could add additional filtering here
abundances <- abundances[!grepl("LC",row.names(abundances)),]
head(abundances)
tpm_data <- abundances[,1:12]

#Importing count data and subsetting in the same way as the tpm data
counts <- read.table("Anther_lines_count.tsv")
head(counts)
counts <- counts[rownames(counts) %in% rownames(tpm_data),]

#Providing metadata to the count matrices
samples <- read.csv("sample_index.csv")
rownames(samples) <- samples$SampleName
samples <- samples[,c("Treatment","Type")]
samples$Treatment <- factor(samples$Treatment)
samples$Type <- factor(samples$Type)

#Check if metadata is in the correct order - these commands should both give output 'TRUE'
all(rownames(samples) == colnames(counts))

##DESeq2 Analysis according to Treatment factor
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=samples, design=~Treatment) #for design=, put the factor based on which you want to make comparisons
dds$Treatment <- relevel(dds$Treatment,"Unsprayed") #set unsprayed as reference level
dds <- DESeq(dds)

#Normalize and transform data
dds_norm <- vst(dds) 
head(assay(dds_norm))

#PCA plot of samples
plotPCA(dds_norm, intgroup = "Treatment")
pca_results <- plotPCA(dds_norm, intgroup = "Treatment", returnData = TRUE) # save pca results and values

##separate tpm data to do different contrasts and do DESeq2 separately for each contrast as
#there is so much variation in all samples which can be seen from PCA plots therefore it is better
#to do DESeq2 separately by selecting only data for that contrast samples each time

#well-watered (WW) (or also named as benchmark (BM) in raw data) and unsprayed (US) comparison (contrast 1)
tpm_data_con1 <- subset(tpm_data[,4:9],)

#add maxtmp columns in the end of tpm_data table from each treatment to show maximum tpm from the each treatment reps
tpm_data_con1$maxtpm_WW <- apply(tpm_data_con1[,1:3],1,max)
tpm_data_con1$maxtpm_US <- apply(tpm_data_con1[,4:6],1,max)

#now select only genes rows with > 0.5 maxtpm in WW 
genes_0.5tpm_WW <- tpm_data_con1[tpm_data_con1$maxtpm_WW > 0.5,]
head (genes_0.5tpm_WW)
dim(genes_0.5tpm_WW)

#now select only genes rows with > 0.5 maxtpm in US 
genes_0.5tpm_US <- tpm_data_con1[tpm_data_con1$maxtpm_US > 0.5,]
head (genes_0.5tpm_US)
dim(genes_0.5tpm_US)

#Merge rows of > 0.5 WW and US
genes_0.5tpm_WWandUS <- transform(merge(genes_0.5tpm_WW, genes_0.5tpm_US, all=TRUE, by='row.names'), row.names=Row.names, Row.names=NULL)
head (genes_0.5tpm_WWandUS)
dim(genes_0.5tpm_WWandUS)

setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results")

#Write a csv file for genes with >0.5tpm merged from WW and US samples 
write.csv(genes_0.5tpm_WWandUS, file="genes_0.5tpm_WWandUS.csv")

setwd("~/Desktop/RNASeq_data_analysis_files/kallisto_results")

#Importing count data and subsetting in the same way as the tpm data
counts <- read.table("Anther_lines_count.tsv")
head(counts)
counts_con1 <- subset(counts[,4:9],)

getwd()

counts <- counts_con1[rownames(counts_con1) %in% rownames(genes_0.5tpm_WWandUS),]

#Providing metadata to the count matrices
samples <- read.csv("sample_index.csv")
rownames(samples) <- samples$SampleName
samples <- subset (samples[4:9,])
samples <- samples[,c("Treatment","Type")]
samples$Treatment <- factor(samples$Treatment)
samples$Type <- factor(samples$Type)

#Check if metadata is in the correct order - these commands should both give output 'TRUE'
all(rownames(samples) == colnames(counts))

##DESeq2 Analysis according to Treatment factor
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=samples, design=~Treatment) 
dds$Treatment <- relevel(dds$Treatment,"Unsprayed")# setting unsprayed as reference level
dds <- DESeq(dds)
head(counts(dds))

#res_WW_v_US <-results(dds)
res_WW_v_US <- results(dds, contrast=c("Treatment","Well-watered","Unsprayed")) ##comparisons of well-watered and unsprayed samples (contrast 1)
ordered_res_WW_v_US <- res_WW_v_US[order(res_WW_v_US$padj),]
head (ordered_res_WW_v_US)
tail(ordered_res_WW_v_US)
dim(ordered_res_WW_v_US)

#remove gene with NA value
ordered_res_na.rm_WWvUS <- na.omit(ordered_res_WW_v_US)
head(ordered_res_na.rm_WWvUS)
tail(ordered_res_na.rm_WWvUS)
dim(ordered_res_na.rm_WWvUS)

WW_vs_US_results_na.rm <- (ordered_res_na.rm_WWvUS[ordered_res_na.rm_WWvUS$padj < 0.05,])
dim(WW_vs_US_results_na.rm)

setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results/DESeq2_results_Final_2_newfilter")

#Write a csv file with adjusted p value < 0.05 
write.csv(WW_vs_US_results_na.rm, file="WW_vs_US_results_0.05.csv")

#filter to keep only upregulated genes that are significantly upregulated (padj<0.05)
DEgenes_WWvUS_up_0.05 <- subset(WW_vs_US_results_na.rm, log2FoldChange > 1)
head(DEgenes_WWvUS_up_0.05)
dim(DEgenes_WWvUS_up_0.05)
write.csv(DEgenes_WWvUS_up_0.05, file="DEgenes_WWvsUS_upreg2fold_0.05.csv")

#filter to keep only downregulated genes that are significantly downregulated (padj<0.05)
DEgenes_WWvUS_down_0.05 <- subset(WW_vs_US_results_na.rm, log2FoldChange < -1)
head(DEgenes_WWvUS_down_0.05)
dim(DEgenes_WWvUS_down_0.05)
write.csv(DEgenes_WWvUS_down_0.05, file="DEgenes_WWvsUS_downreg2fold_0.05.csv")


##ABA and unsprayed samples comparison (contrast 2)
tpm_data_con2 <- subset(tpm_data[,c(1:3,7:9)])
tpm_data_con2$maxtpm_ABA <- apply(tpm_data_con2[,1:3],1,max)
tpm_data_con2$maxtpm_US <- apply(tpm_data_con2[,4:6],1,max)

#now select only genes rows with >0.5 maxtpm in ABA
genes_0.5tpm_ABA <- tpm_data_con2[tpm_data_con2$maxtpm_ABA > 0.5,]
head (genes_0.5tpm_ABA)

#now select only genes rows with > 0.5 maxtpm in US 
genes_0.5tpm_US <- tpm_data_con2[tpm_data_con2$maxtpm_US > 0.5,]
head (genes_0.5tpm_US)
dim(genes_0.5tpm_US)

#Merge rows of >0.5 ABA and US
genes_0.5tpm_ABAandUS <- transform(merge(genes_0.5tpm_ABA, genes_0.5tpm_US, all=TRUE, by='row.names'), row.names=Row.names, Row.names=NULL)
head (genes_0.5tpm_ABAandUS)
dim(genes_0.5tpm_ABAandUS)

setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results")

#Write a csv file for genes with 0.5tpm merged from ABA and US samples 
write.csv(genes_0.5tpm_ABAandUS, file="genes_0.5tpm_ABAandUS.csv")

#Importing count data and subsetting in the same way as the tpm data
setwd("~/Desktop/RNASeq_data_analysis_files/kallisto_results")
counts <- read.table("Anther_lines_count.tsv")
head(counts)
counts_con2 <- subset(counts[,c(1:3,7:9)])

counts <- counts_con2[rownames(counts_con2) %in% rownames(genes_0.5tpm_ABAandUS),]

#Providing metadata to the count matrices
samples <- read.csv("sample_index.csv")
rownames(samples) <- samples$SampleName
samples <- subset (samples[c(1:3,7:9),])
samples <- samples[,c("Treatment","Type")]
samples$Treatment <- factor(samples$Treatment)
samples$Type <- factor(samples$Type)

#Check if metadata is in the correct order - these commands should both give output 'TRUE'
all(rownames(samples) == colnames(counts))

##DESeq2 Analysis according to Treatment factor
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=samples, design=~Treatment) #for design=, put the factor based on which you want to make comparisons
dds$Treatment <- relevel(dds$Treatment,"Unsprayed") #set unsprayed as reference level
dds <- DESeq(dds)
head(counts(dds))

res_ABA_v_US <- results(dds, contrast=c("Treatment","ABA","Unsprayed")) #comparisons of ABA and unsprayed samples (contrast 2)
ordered_res_ABA_v_US <- res_ABA_v_US[order(res_ABA_v_US$padj),]
head (ordered_res_ABA_v_US)
tail(ordered_res_ABA_v_US)
dim(ordered_res_ABA_v_US)

#remove gene with NA value
ordered_res_na.rm_ABAvUS <- na.omit(ordered_res_ABA_v_US)
head(ordered_res_na.rm_ABAvUS)
tail(ordered_res_na.rm_ABAvUS)
dim(ordered_res_na.rm_ABAvUS)

ABA_vs_US_results_na.rm <- (ordered_res_na.rm_ABAvUS[ordered_res_na.rm_ABAvUS$padj < 0.05,])
dim(ABA_vs_US_results_na.rm)

setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results/DESeq2_results_Final_2_newfilter/Contrast2_ABAvsUS")

#Write a csv file with adjusted p value < 0.05 
write.csv(ABA_vs_US_results_na.rm, file="ABA_vs_US_results_0.05.csv")

#filter to keep only upregulated genes that are significantly upregulated (padj<0.05)
DEgenes_ABAvUS_up_0.05 <- subset(ABA_vs_US_results_na.rm, log2FoldChange > 1)
head(DEgenes_ABAvUS_up_0.05)
dim(DEgenes_ABAvUS_up_0.05)
write.csv(DEgenes_ABAvUS_up_0.05, file="DEgenes_ABAvsUS_upreg2fold_0.05.csv")

#filter to keep only downregulated genes that are significantly downregulated (padj<0.05)
DEgenes_ABAvUS_down_0.05 <- subset(ABA_vs_US_results_na.rm, log2FoldChange < -1)
head(DEgenes_ABAvUS_down_0.05)
dim(DEgenes_ABAvUS_down_0.05)
write.csv(DEgenes_ABAvUS_down_0.05, file="DEgenes_ABAvsUS_downreg2fold_0.05.csv")


##comparison of Vapor Gard (VG) and unsprayed samples (contrast 3)
tpm_data_con3 <- subset(tpm_data[,7:12])
tpm_data_con3$maxtpm_US <- apply(tpm_data_con3[,1:3],1,max)
tpm_data_con3$maxtpm_VG <- apply(tpm_data_con3[,4:6],1,max)

#now select only genes rows with >0.5 maxtpm in VG
genes_0.5tpm_VG <- tpm_data_con3[tpm_data_con3$maxtpm_VG > 0.5,]
head (genes_0.5tpm_VG)
dim(genes_0.5tpm_VG)
#now select only genes rows with > 0.5 maxtpm in US 
genes_0.5tpm_US <- tpm_data_con3[tpm_data_con3$maxtpm_US > 0.5,]
head (genes_0.5tpm_US)
dim(genes_0.5tpm_US)

#Merge rows of >0.5 VG and US
genes_0.5tpm_VGandUS <- transform(merge(genes_0.5tpm_VG, genes_0.5tpm_US, all=TRUE, by='row.names'), row.names=Row.names, Row.names=NULL)
head (genes_0.5tpm_VGandUS)
dim(genes_0.5tpm_VGandUS)

setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results/DESeq2_results_Final_2_newfilter")

#Write a csv file for genes with 0.5tpm merged from VG and US samples 
write.csv(genes_0.5tpm_VGandUS, file="genes_0.5tpm_VGandUS.csv")

#Importing count data and subsetting in the same way as the tpm data
setwd("~/Desktop/RNASeq_data_analysis_files/kallisto_results")
counts <- read.table("Anther_lines_count.tsv")
head(counts)
counts_con3 <- subset(counts[,7:12])

counts <- counts_con3[rownames(counts_con3) %in% rownames(genes_0.5tpm_VGandUS),]

#Providing metadata to the count matrices
samples <- read.csv("sample_index.csv")
rownames(samples) <- samples$SampleName
samples <- subset (samples[7:12,])
samples <- samples[,c("Treatment","Type")]
samples$Treatment <- factor(samples$Treatment)
samples$Type <- factor(samples$Type)

#Check if metadata is in the correct order - these commands should both give output 'TRUE'
all(rownames(samples) == colnames(counts))

##DESeq2 Analysis according to Treatment factor
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=samples, design=~Treatment) #for design=, put the factor based on which you want to make comparisons
dds$Treatment <- relevel(dds$Treatment,"Unsprayed") #set unsprayed as reference level
dds <- DESeq(dds)
head(counts(dds))

res_VG_v_US <- results(dds, contrast=c("Treatment","VG","Unsprayed"))
ordered_res_VG_v_US <- res_VG_v_US[order(res_VG_v_US$padj),]
head (ordered_res_VG_v_US)
tail(ordered_res_VG_v_US)
dim(ordered_res_VG_v_US)

#remove gene with NA value
ordered_res_na.rm_VGvUS <- na.omit(ordered_res_VG_v_US)
head(ordered_res_na.rm_VGvUS)
tail(ordered_res_na.rm_VGvUS)
dim(ordered_res_na.rm_VGvUS)

VG_vs_US_results_na.rm <- (ordered_res_na.rm_VGvUS[ordered_res_na.rm_VGvUS$padj < 0.05,])
dim(VG_vs_US_results_na.rm)

setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results/DESeq2_results_Final_2_newfilter/Contrast3_VGvsUS/")

#Write a csv file with adjusted p value < 0.05 
write.csv(VG_vs_US_results_na.rm, file="VG_vs_US_results_0.05.csv")

#filter to keep only upregulated genes that are significantly upregulated (padj<0.05)
DEgenes_VGvUS_up_0.05 <- subset(VG_vs_US_results_na.rm, log2FoldChange > 1)
head(DEgenes_VGvUS_up_0.05)
dim(DEgenes_VGvUS_up_0.05)
write.csv(DEgenes_VGvUS_up_0.05, file="DEgenes_VGvsUS_upreg2fold_0.05.csv")

#filter to keep only downregulated genes that are significantly downregulated (padj<0.05)
DEgenes_VGvUS_down_0.05 <- subset(VG_vs_US_results_na.rm, log2FoldChange < -1)
head(DEgenes_VGvUS_down_0.05)
dim(DEgenes_VGvUS_down_0.05)
write.csv(DEgenes_VGvUS_down_0.05, file="DEgenes_VGvsUS_downreg2fold_0.05.csv")
