## Aim to do GO enrichment analysis on DE genes found in four different contrasts of meiosis anther data
## Philippa's script edited by Misbah Sehar
## 10th May 2023
## Updated 12th July 2023 

#### read in information about lengths and GO terms #########
# read in GO terms
all_go <- read.table("~/Desktop/RNASeq_data_analysis_files/scripts/IWGSC_stress_GO.txt", sep=",")
head(all_go)
all_go <- all_go[,c(1,2)]
colnames(all_go) <- c("Gene", "GO_term")
head(all_go)
dim(all_go)

# convert from v1.0 to v1.1 
head(gsub("01G", "02G", all_go$Gene))

all_go$Gene <- (gsub("01G", "02G", all_go$Gene))
head(all_go)
dim(all_go)

all_go_HC <- all_go[!grepl("LC", all_go$Gene),]
head(all_go_HC)
dim(all_go_HC)

length(unique(all_go_HC$Gene)) # number of HC genes with go terms before removing ones which don't match v1.0 to v1.1

# only keep genes which were >99 % ID > 90% coverage from v1.0 to v1.1 
genes_to_transfer <- read.table("~/Desktop/RNASeq_data_analysis_files/scripts/genes_to_transfer_qcov90_pident99_same_ID.txt", sep=",")
head(genes_to_transfer)
colnames(genes_to_transfer) <- c("gene_v1.1", "gene_v1.0")

all_go <- all_go[all_go$Gene %in% genes_to_transfer$gene_v1.1,]
head(all_go)
dim(all_go)

length(unique(all_go$Gene)) # number of genes with go terms

# select only genes which were used for DESeq2 in contrast 1
setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results/DESeq2_results_Final_2_newfilter/")
genes_0.5tpm_WWandUS <- read.csv("genes_0.5tpm_WWandUS.csv")

head(genes_0.5tpm_WWandUS)
dim(genes_0.5tpm_WWandUS)

all_go <- subset(all_go, Gene %in% genes_0.5tpm_WWandUS$X)
dim(all_go)

length(unique(genes_0.5tpm_WWandUS$X)) # number of genes expressed

length(unique(all_go$Gene)) #  number of genes with go terms which were expressed

#create vector for gene_lengths

# need to get lengths of genes not of transcripts
setwd("~/Desktop/RNASeq_data_analysis_files/kallisto_results")
lengths <- read.csv("Anther_gene_length.csv", header=T)
head(lengths)
colnames(lengths) <- c("gene", "length")
head(lengths)

t1 <- subset(lengths, gene %in% genes_0.5tpm_WWandUS$X)
head(t1)
dim(t1)

# turn into a vector called gene.lens to use with GOSeq
gene.lens <- as.numeric(t1$length)
names(gene.lens) = t1$gene
head(gene.lens)
length(gene.lens)


####### Do GO term enrichment ####

assayed.genes <- as.vector(t1$gene)
length(assayed.genes)

library(goseq)
library(BiasedUrn)
library(geneLenDataBase)

# do the GO term enrichment for first contrast (WW vs US)
setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results/DESeq2_results_Final_2_newfilter/Contrast1_WWvsUS/")

my_data_con1 <- read.csv("WW_vs_US_results_0.05.csv")
head(my_data_con1)
dim(my_data_con1)

upreg2fold_0.05 <- (my_data_con1[my_data_con1$padj < 0.05 & my_data_con1$log2FoldChange > 1,])
downreg2fold_0.05 <- (my_data_con1[my_data_con1$padj < 0.05 & my_data_con1$log2FoldChange < -1,])

# Select separate table each time (for up and downreg. genes) to do GO analysis on each and save separate files  
genes_for_GO <- upreg2fold_0.05

genes_for_GO <- downreg2fold_0.05

head(genes_for_GO)
dim(genes_for_GO)

#now do GO stats analysis on the genes expressed in each pattern compared to all genes expressed
#create a named binary vector for genes where one means differentially expressed and 0 means not differentially expressed
de.genes <- genes_for_GO$X
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
head(gene.vector)

#now carry out the GOseq analysis
pwf = nullp(gene.vector, bias.data = gene.lens, plot.fit = TRUE)
GO.wall = goseq(pwf, gene2cat = all_go) #this gave table with p-values

# add new column with over represented GO terms padj
GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BH")

GO_enriched_contrast <- GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",]
head(GO_enriched_contrast)
dim(GO_enriched_contrast)

setwd("~/Desktop/RNASeq_data_analysis_files/GO_analysis_files_newfilter/GO_analysis_con1/")

### create individuals files for 0.05 up & down regulated genes in this contrast 
write.table(GO.wall[GO.wall$over_rep_padj <0.05,], file = paste0("con1_upreg2fold0.05_GO.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
write.table(GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",], file = paste0("con1_upreg2fold0.05_GO_BP.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)

write.table(GO.wall[GO.wall$over_rep_padj <0.05,], file = paste0("con1_downreg2fold0.05_GO.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
write.table(GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",], file = paste0("con1_downreg2fold0.05_GO_BP.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)


## Start GO analysis again from the start for contrast 2 (ABA vs US)

#### read in information about lengths and GO terms #########
# read in GO terms
all_go <- read.table("~/Desktop/RNASeq_data_analysis_files/scripts/IWGSC_stress_GO.txt", sep=",")
head(all_go)
all_go <- all_go[,c(1,2)]
colnames(all_go) <- c("Gene", "GO_term")
head(all_go)
dim(all_go)

# convert from v1.0 to v1.1 
head(gsub("01G", "02G", all_go$Gene))

all_go$Gene <- (gsub("01G", "02G", all_go$Gene))
head(all_go)
dim(all_go)

all_go_HC <- all_go[!grepl("LC", all_go$Gene),]
head(all_go_HC)
dim(all_go_HC)

length(unique(all_go_HC$Gene)) # number of HC genes with go terms before removing ones which don't match v1.0 to v1.1

# only keep genes which were >99 % ID > 90% coverage from v1.0 to v1.1 
genes_to_transfer <- read.table("~/Desktop/RNASeq_data_analysis_files/scripts/genes_to_transfer_qcov90_pident99_same_ID.txt", sep=",")
head(genes_to_transfer)
colnames(genes_to_transfer) <- c("gene_v1.1", "gene_v1.0")

all_go <- all_go[all_go$Gene %in% genes_to_transfer$gene_v1.1,]
head(all_go)
dim(all_go)

length(unique(all_go$Gene)) # number of genes with go terms

# select only genes which were used for DESeq2 in contrast 2
setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results/DESeq2_results_Final_2_newfilter/")
genes_0.5tpm_ABAandUS <- read.csv("genes_0.5tpm_ABAandUS.csv")

head(genes_0.5tpm_ABAandUS)
dim(genes_0.5tpm_ABAandUS)

all_go <- subset(all_go, Gene %in% genes_0.5tpm_ABAandUS$X)
dim(all_go)

length(unique(genes_0.5tpm_ABAandUS$X)) # number of genes expressed

length(unique(all_go$Gene)) #  number of genes with go terms which were expressed

#create vector for gene_lengths

# need to get lengths of genes not of transcripts
setwd("~/Desktop/RNASeq_data_analysis_files/kallisto_results")
lengths <- read.csv("Anther_gene_length.csv", header=T)
head(lengths)
colnames(lengths) <- c("gene", "length")
head(lengths)

t1 <- subset(lengths, gene %in% genes_0.5tpm_ABAandUS$X)
head(t1)
dim(t1)

# turn into a vector called gene.lens to use with GOSeq
gene.lens <- as.numeric(t1$length)
names(gene.lens) = t1$gene
head(gene.lens)
length(gene.lens)


####### Do GO term enrichment ####

assayed.genes <- as.vector(t1$gene)
length(assayed.genes)

# do the GO term enrichment for contrast 2
setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results/DESeq2_results_Final_2_newfilter/Contrast2_ABAvsUS/")

my_data_con2 <- read.csv("ABA_vs_US_results_0.05.csv")
head(my_data_con2)

upreg2fold_0.05 <- (my_data_con2[my_data_con2$padj < 0.05 & my_data_con2$log2FoldChange > 1,])
downreg2fold_0.05 <- (my_data_con2[my_data_con2$padj < 0.05 & my_data_con2$log2FoldChange < -1,])

# Select separate table each time (for up and downreg. genes) to do GO analysis on each and save separate files  
genes_for_GO <- upreg2fold_0.05

genes_for_GO <- downreg2fold_0.05

head(genes_for_GO)
dim(genes_for_GO)

#now do GO stats analysis on the genes expressed in each pattern compared to all genes expressed
#create a named binary vector for genes where one means differentially expressed and 0 means not differentially expressed
de.genes <- genes_for_GO$X
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
head(gene.vector)

#now carry out the GOseq analysis
pwf = nullp(gene.vector, bias.data = gene.lens, plot.fit = TRUE)
GO.wall = goseq(pwf, gene2cat = all_go) #this gave table with p-values

# add new column with over represented GO terms padj
GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BH")

GO_enriched_contrast <- GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",]
head(GO_enriched_contrast)
dim(GO_enriched_contrast)

setwd("~/Desktop/RNASeq_data_analysis_files/GO_analysis_files_newfilter/GO_analysis_con2/")

### create individuals files for 0.05 up & down regulated genes in this contrast 
write.table(GO.wall[GO.wall$over_rep_padj <0.05,], file = paste0("con2_upreg2fold0.05_GO.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
write.table(GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",], file = paste0("con2_upreg2fold0.05_GO_BP.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)

write.table(GO.wall[GO.wall$over_rep_padj <0.05,], file = paste0("con2_downreg2fold0.05_GO.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
write.table(GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",], file = paste0("con2_downreg2fold0.05_GO_BP.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)


## Start GO analysis again from the start for contrast 3 (VG vs US)

#### read in information about lengths and GO terms #########
# read in GO terms
all_go <- read.table("~/Desktop/RNASeq_data_analysis_files/scripts/IWGSC_stress_GO.txt", sep=",")
head(all_go)
all_go <- all_go[,c(1,2)]
colnames(all_go) <- c("Gene", "GO_term")
head(all_go)
dim(all_go)

# convert from v1.0 to v1.1 
head(gsub("01G", "02G", all_go$Gene))

all_go$Gene <- (gsub("01G", "02G", all_go$Gene))
head(all_go)
dim(all_go)

all_go_HC <- all_go[!grepl("LC", all_go$Gene),]
head(all_go_HC)
dim(all_go_HC)

length(unique(all_go_HC$Gene)) # number of HC genes with go terms before removing ones which don't match v1.0 to v1.1

# only keep genes which were >99 % ID > 90% coverage from v1.0 to v1.1 
genes_to_transfer <- read.table("~/Desktop/RNASeq_data_analysis_files/scripts/genes_to_transfer_qcov90_pident99_same_ID.txt", sep=",")
head(genes_to_transfer)
colnames(genes_to_transfer) <- c("gene_v1.1", "gene_v1.0")

all_go <- all_go[all_go$Gene %in% genes_to_transfer$gene_v1.1,]
head(all_go)
dim(all_go)

length(unique(all_go$Gene)) # number of genes with go terms

# select only genes which were used for DESeq2 in contrast 3
setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results/DESeq2_results_Final_2_newfilter/")
genes_0.5tpm_VGandUS <- read.csv("genes_0.5tpm_VGandUS.csv")

head(genes_0.5tpm_VGandUS)
dim(genes_0.5tpm_VGandUS)

all_go <- subset(all_go, Gene %in% genes_0.5tpm_VGandUS$X)
dim(all_go)

length(unique(genes_0.5tpm_VGandUS$X)) # number of genes expressed

length(unique(all_go$Gene)) #  number of genes with go terms which were expressed

#create vector for gene_lengths

# need to get lengths of genes not of transcripts
setwd("~/Desktop/RNASeq_data_analysis_files/kallisto_results")
lengths <- read.csv("Anther_gene_length.csv", header=T)
head(lengths)
colnames(lengths) <- c("gene", "length")
head(lengths)

t1 <- subset(lengths, gene %in% genes_0.5tpm_VGandUS$X)
head(t1)
dim(t1)

# turn into a vector called gene.lens to use with GOSeq
gene.lens <- as.numeric(t1$length)
names(gene.lens) = t1$gene
head(gene.lens)
length(gene.lens)


####### Do GO term enrichment ####

assayed.genes <- as.vector(t1$gene)
length(assayed.genes)

# do the GO term enrichment for contrast 3
setwd("~/Desktop/RNASeq_data_analysis_files/DESeq2_results/DESeq2_results_Final_2_newfilter/Contrast3_VGvsUS/")

my_data_con3 <- read.csv("VG_vs_US_results_0.05.csv")
head(my_data_con3)

upreg2fold_0.05 <- (my_data_con3[my_data_con3$padj < 0.05 & my_data_con3$log2FoldChange > 1,])
downreg2fold_0.05 <- (my_data_con3[my_data_con3$padj < 0.05 & my_data_con3$log2FoldChange < -1,])

# Select separate table each time (for up and downreg. genes) to do GO analysis on each and save separate files 
genes_for_GO <- upreg2fold_0.05

genes_for_GO <- downreg2fold_0.05

head(genes_for_GO)
dim(genes_for_GO)

#now do GO stats analysis on the genes expressed in each pattern compared to all genes expressed
#create a named binary vector for genes where one means differentially expressed and 0 means not differentially expressed
de.genes <- genes_for_GO$X
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
head(gene.vector)

#now carry out the GOseq analysis
pwf = nullp(gene.vector, bias.data = gene.lens, plot.fit = TRUE)
GO.wall = goseq(pwf, gene2cat = all_go) #this gave table with p-values

# add new column with over represented GO terms padj
GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BH")

GO_enriched_contrast <- GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",]
head(GO_enriched_contrast)
dim(GO_enriched_contrast)

setwd("~/Desktop/RNASeq_data_analysis_files/GO_analysis_files_newfilter/GO_analysis_con3/")

### create individuals files for 0.05 up & down regulated genes in this contrast 
write.table(GO.wall[GO.wall$over_rep_padj <0.05,], file = paste0("con3_upreg2fold0.05_GO.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
write.table(GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",], file = paste0("con3_upreg2fold0.05_GO_BP.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)

write.table(GO.wall[GO.wall$over_rep_padj <0.05,], file = paste0("con3_downreg2fold0.05_GO.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
write.table(GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",], file = paste0("con3_downreg2fold0.05_GO_BP.tsv", sep =""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
