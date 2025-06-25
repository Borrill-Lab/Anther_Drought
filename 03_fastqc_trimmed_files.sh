#!/bin/bash
#SBATCH -p nbi-medium
#SBATCH --mem=120000
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Misbah/fastqc_results_trim/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Misbah/fastqc_results_trim/slurm_output/%x.%N.%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=jek23cen@nbi.ac.uk

source fastqc-0.11.8

cd /jic/scratch/groups/Philippa-Borrill/Misbah/trim_results

fastqc -o /jic/scratch/groups/Philippa-Borrill/Misbah/fastqc_results_trim *.fq.gz