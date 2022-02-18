#Aaron Horning
#!/bin/bash
#SBATCH --job-name=aaron
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 #change to 8?
#SBATCH --time=10-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=ahorning@stanford.edu
#SBATCH --account=mpsnyder

CHR=$1

######## to retrieve info before run:
date
echo "This job ran mutect2 on A001 samples and chromosome $CHR"
echo "I ran on host: $(hostname -s)"

echo "SGE Environment is:"
env | grep "SGE" | sort

echo "My limits are:"
ulimit -a


#this step is for muTect scanning

module load gatk4
cd /home/ahorning/DNAseq/Bulk_Batch123/Mutect2/withWES
# /scg/apps/data/gatk4/Mutect2/af-only-gnomad.hg38.vcf.gz

# GATK Mutect2 Somatic Variant calling for Multiple Samples
# A001
gatk Mutect2 \
     -R /home/ahorning/DNAseq/Reference/Homo_sapiens_assembly38.fasta \
     -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001_blood/02_MAPPING/A001_blood.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C004/02_MAPPING/A001C004.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C005/02_MAPPING/A001C005.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C007/02_MAPPING/A001C007.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C021/02_MAPPING/A001C021.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C102/02_MAPPING/A001C102.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C107/02_MAPPING/A001C107.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C122/02_MAPPING/A001C122.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C124/02_MAPPING/A001C124.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C210/02_MAPPING/A001C210.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C214/02_MAPPING/A001C214.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C219/02_MAPPING/A001C219.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_Batch123/root/A001C222/02_MAPPING/A001C222.sorted.ir.br.rmDup.md.bam -I /home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WES/A001C004/02_MAPPING/A001C004.sorted.ir.br.rmDup.md.rn.bam -I /home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WES/A001C005/02_MAPPING/A001C005.sorted.ir.br.rmDup.md.rn.bam -I /home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WES/A001C021/02_MAPPING/A001C021.sorted.ir.br.rmDup.md.rn.bam -I /home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WES/A001C102/02_MAPPING/A001C102.sorted.ir.br.rmDup.md.rn.bam -I /home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WES/A001C107/02_MAPPING/A001C107.sorted.ir.br.rmDup.md.rn.bam -I /home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WES/A001C122/02_MAPPING/A001C122.sorted.ir.br.rmDup.md.rn.bam -I /home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WES/A001C210/02_MAPPING/A001C210.sorted.ir.br.rmDup.md.rn.bam -I /home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WES/A001C214/02_MAPPING/A001C214.sorted.ir.br.rmDup.md.rn.bam -I /home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WES/A001C219/02_MAPPING/A001C219.sorted.ir.br.rmDup.md.rn.bam -I /home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WES/A001C222/02_MAPPING/A001C222.sorted.ir.br.rmDup.md.rn.bam \
     -normal A001_blood\
     --germline-resource /scg/apps/data/gatk4/Mutect2/af-only-gnomad.hg38.vcf.gz \
     -L $CHR \
     --f1r2-tar-gz /home/ahorning/DNAseq/Bulk_Batch123/Mutect2/withWES/A001-${CHR}-f1r2.tar.gz \
     -O /home/ahorning/DNAseq/Bulk_Batch123/Mutect2/withWES/A001.${CHR}.mutect2.somatic.unfiltered.vcf.gz

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
