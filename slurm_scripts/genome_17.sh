#!/bin/bash
#SBATCH --partition=shared
#SBATCH -A jhu152
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --output=slurm_output_17.out
#SBATCH --job-name=genome_17
#SBATCH --time=24:00:00
conda activate mwas
Rscript /expanse/lustre/projects/jhu152/naglemi/mwas/CpGWAS/scripts/deploy_gwas2mwas_batch.R --genome_file_index=17

