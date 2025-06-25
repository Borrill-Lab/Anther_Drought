#!/bin/bash
#SBATCH -p nbi-medium
#SBATCH --mem=30000
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Misbah/kallisto_results/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Misbah/kallisto_results/slurm_output/%x.%N.%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=jek23cen@nbi.ac.uk

cd /jic/scratch/groups/Philippa-Borrill/Misbah/trim_results

source package /nbi/software/testing/bin/kallisto-0.46.1
source package aeee87c4-1923-4732-aca2-f2aff23580cc
kallisto quant -i /jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_ALL_20170706_transcripts.fasta_index -o /jic/scratch/groups/Philippa-Borrill/Misbah/kallisto_results//Sample_ABA2 -t 8 --pseudobam Sample_ABA2_1.paired.fq.gz Sample_ABA2_2.paired.fq.gz
samtools sort -T /jic/scratch/groups/Philippa-Borrill/Misbah/kallisto_results//Sample_ABA2/Sample_ABA2.bam -O bam -o /jic/scratch/groups/Philippa-Borrill/Misbah/kallisto_results//Sample_ABA2.sorted.bam -@ 8 /jic/scratch/groups/Philippa-Borrill/Misbah/kallisto_results//Sample_ABA2/pseudoalignments.bam
samtools index /jic/scratch/groups/Philippa-Borrill/Misbah/kallisto_results//Sample_ABA2.sorted.bam /jic/scratch/groups/Philippa-Borrill/Misbah/kallisto_results//Sample_ABA2.sorted.bai
rm /jic/scratch/groups/Philippa-Borrill/Misbah/kallisto_results//Sample_ABA2/pseudoalignments.bam
