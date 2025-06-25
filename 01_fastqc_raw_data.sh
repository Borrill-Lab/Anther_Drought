#!/bin/bash
#SBATCH -p nbi-medium
#SBATCH --mem=1G
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Misbah/fastqc_results/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Misbah/fastqc_results/slurm_output/%x.%N.%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=jek23cen@nbi.ac.uk

source fastqc-0.11.8

cd /jic/research-groups/Philippa-Borrill/raw_data/Misbah/40-794093267/00_fastq

fastqc -o /jic/scratch/groups/Philippa-Borrill/Misbah/fastqc_results *.fastq.gz