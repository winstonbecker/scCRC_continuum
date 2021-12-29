#!/bin/bash
#
# Example bash Script for Running Cell Ranger with pre-mRNA reference using slurm job scheduler 
#################
# set a job name  
#SBATCH --job-name=RNA_A001-C-014
#################  
# filename for job output
#SBATCH --output=RNA_A001-C-014.out.txt
#################
# filename for job errors
#SBATCH --error=RNA_A001-C-014.err
#################
# time requested in minutes in this case, hh:mm:ss
#SBATCH --time=24:00:00
#################
# number of nodes requested
#SBATCH --nodes=1
#################
# partition
#SBATCH --partition=nnn
#################
# for local parallelism
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#################
export PATH=./cellranger-3.1.0:$PATH
srun cellranger count --id=A001-C-014 \
                   --transcriptome=./GRCh38-3.0.0_premrna/ \
                   --fastqs=./fastq_path/ \
                   --sample=A001-C-014  \
                   --expect-cells=9000
