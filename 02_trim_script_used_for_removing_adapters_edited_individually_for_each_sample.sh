#!/bin/bash
#SBATCH -p nbi-medium
#SBATCH --mem=120000
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Misbah/trim_results/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Misbah/trim_results/slurm_output/%x.%N.%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=jek23cen@nbi.ac.uk

cd /jic/research-groups/Philippa-Borrill/raw_data/Misbah/40-794093267/00_fastq

source package 50fcf79b-73a3-4f94-9553-5ed917823423
trimmomatic PE -phred33 -threads 4 ABA5_R1_001.fastq.gz ABA5_R2_001.fastq.gz /jic/scratch/groups/Philippa-Borrill/Misbah/trim_results//Sample_ABA5_1.paired.fq.gz /jic/scratch/groups/Philippa-Borrill/Misbah/trim_results//Sample_ABA5_1.unpaired.fq.gz /jic/scratch/groups/Philippa-Borrill/Misbah/trim_results//Sample_ABA5_2.paired.fq.gz /jic/scratch/groups/Philippa-Borrill/Misbah/trim_results//Sample_ABA5_2.unpaired.fq.gz  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:80
