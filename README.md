## Wheat_Anther_Drought_Antitranspirants

## Scripts
### Processing of raw fastq files:

#### 01_fastq_raw_data.sh
To assess quality of raw fastq files.
#### 02_trim_script_used_for_removing_adapters_edited_individually_for_each_sample.sh
To remove adapters from the raw fastq files.
#### 03_fastqc_trimmed_files.sh
Assessing quality of trimmed files.

### Mapping and combining into one dataframe:

#### 04_kallisto_script_used_for_mapping_to_RefSeqv1.1_edited_individually_for_each_sample.sh
Mapping of anther samples to RefSeqv1.1 transcriptome.
#### 05_tximport_counts_tpm_per_gene.R
Summarise data at gene level to get counts per gene and tpm per gene data into single tables for all samples.

### Identifying differentially expressed genes and enriched GO terms:

#### 06_differential_expression_DESeq2.R
Identifying differentially expressed genes in different contrast comparisons. 
#### 07_GO_term_enrichment_analysis.R
Identifying GO enriched terms in different contrast comparions.

## Data files used in scripts

TruSeq3-PE-2.fa - Adapter file used in the trimming step

List of transcripts to gene conversion for v1.1 can be found here: https://github.com/Borrill-Lab/WheatFlagLeafSenescence/blob/master/data/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt

v1.0 GO terms can be found here: https://github.com/Borrill-Lab/WheatFlagLeafSenescence/blob/master/data/IWGSC_stress_GO.csv

List of transcripts for which GO terms were transfered from v1.0 to v1.1 can be found here: https://github.com/Borrill-Lab/WheatFlagLeafSenescence/blob/master/data/genes_to_transfer_qcov90_pident99_same_ID.csv


