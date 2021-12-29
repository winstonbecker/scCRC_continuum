#!/bin/bash
#
# Example bash Script for Running Cell Ranger ATAC using slurm job scheduler 
#################
# set a job name  
#SBATCH --job-name=ATAC_CRC-1-8810
#################  
# filename for job output
#SBATCH --output=ATAC_CRC-1-8810.out.txt
#################
# filename for job errors
#SBATCH --error=ATAC_CRC-1-8810.err
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
ml biology
export PATH=./cellranger-atac-1.2.0:$PATH
srun cellranger-atac count --id=CRC-1-8810-D \
                   --reference=./refdata-cellranger-atac-GRCh38-1.2.0 \
                   --fastqs=./fastq_path/CRC-1-8810-D/ \
                   --sample=CRC-1-8810-D
