#!/bin/bash
#SBATCH --job-name=aaron
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 #change to 8?
#SBATCH --time=10-00:00:00 # big job
#SBATCH --mem=170G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=ahorning@stanford.edu
#SBATCH --account=mpsnyder

# Run this script this way:
# sbatch postMutect2_FilterMutations_perPT.sh A001

PT=$1

# can be run oonce the 67 chromosome intervals from Mutect2 are completed
# run this in the same directory as the results of the GATK_mutect2_${PT}.sh (for example) results.

module load gatk4

# # https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
# # GATK Orientation Bias Model Learning for Filtering later on

# # took 3 hours
# gatk --java-options "-Xmx50g" LearnReadOrientationModel -I ${PT}-chr1:1-50000000-f1r2.tar.gz -I ${PT}-chr1:50000001-100000000-f1r2.tar.gz -I ${PT}-chr1:100000001-150000000-f1r2.tar.gz -I ${PT}-chr1:150000001-200000000-f1r2.tar.gz -I ${PT}-chr1:200000001-248956422-f1r2.tar.gz -I ${PT}-chr2:1-50000000-f1r2.tar.gz -I ${PT}-chr2:50000001-100000000-f1r2.tar.gz -I ${PT}-chr2:100000001-150000000-f1r2.tar.gz -I ${PT}-chr2:150000001-200000000-f1r2.tar.gz -I ${PT}-chr2:200000001-242193529-f1r2.tar.gz -I ${PT}-chr3:1-50000000-f1r2.tar.gz -I ${PT}-chr3:50000001-100000000-f1r2.tar.gz -I ${PT}-chr3:100000001-150000000-f1r2.tar.gz -I ${PT}-chr3:150000001-198295559-f1r2.tar.gz -I ${PT}-chr4:1-50000000-f1r2.tar.gz -I ${PT}-chr4:50000001-100000000-f1r2.tar.gz -I ${PT}-chr4:100000001-150000000-f1r2.tar.gz -I ${PT}-chr4:150000001-190214555-f1r2.tar.gz -I ${PT}-chr5:1-50000000-f1r2.tar.gz -I ${PT}-chr5:50000001-100000000-f1r2.tar.gz -I ${PT}-chr5:100000001-150000000-f1r2.tar.gz -I ${PT}-chr5:150000001-181538259-f1r2.tar.gz -I ${PT}-chr6:1-50000000-f1r2.tar.gz -I ${PT}-chr6:50000001-100000000-f1r2.tar.gz -I ${PT}-chr6:100000001-150000000-f1r2.tar.gz -I ${PT}-chr6:150000001-170805979-f1r2.tar.gz -I ${PT}-chr7:1-50000000-f1r2.tar.gz -I ${PT}-chr7:50000001-100000000-f1r2.tar.gz -I ${PT}-chr7:100000001-150000000-f1r2.tar.gz -I ${PT}-chr7:150000001-159345973-f1r2.tar.gz -I ${PT}-chr8:1-50000000-f1r2.tar.gz -I ${PT}-chr8:50000001-100000000-f1r2.tar.gz -I ${PT}-chr8:100000001-145138636-f1r2.tar.gz -I ${PT}-chr9:1-50000000-f1r2.tar.gz -I ${PT}-chr9:50000001-100000000-f1r2.tar.gz -I ${PT}-chr9:100000001-138394717-f1r2.tar.gz -I ${PT}-chr10:1-50000000-f1r2.tar.gz -I ${PT}-chr10:50000001-100000000-f1r2.tar.gz -I ${PT}-chr10:100000001-133797422-f1r2.tar.gz -I ${PT}-chr11:1-50000000-f1r2.tar.gz -I ${PT}-chr11:50000001-100000000-f1r2.tar.gz -I ${PT}-chr11:100000001-135086622-f1r2.tar.gz -I ${PT}-chr12:1-50000000-f1r2.tar.gz -I ${PT}-chr12:50000001-100000000-f1r2.tar.gz -I ${PT}-chr12:100000001-133275309-f1r2.tar.gz -I ${PT}-chr13:1-50000000-f1r2.tar.gz -I ${PT}-chr13:50000001-100000000-f1r2.tar.gz -I ${PT}-chr13:100000001-114364328-f1r2.tar.gz -I ${PT}-chr14:1-50000000-f1r2.tar.gz -I ${PT}-chr14:50000001-100000000-f1r2.tar.gz -I ${PT}-chr14:100000001-107043718-f1r2.tar.gz -I ${PT}-chr15:1-50000000-f1r2.tar.gz -I ${PT}-chr15:50000001-100000000-f1r2.tar.gz -I ${PT}-chr15:100000001-101991189-f1r2.tar.gz -I ${PT}-chr16:1-50000000-f1r2.tar.gz -I ${PT}-chr16:50000001-90338345-f1r2.tar.gz -I ${PT}-chr17:1-50000000-f1r2.tar.gz -I ${PT}-chr17:50000001-83257441-f1r2.tar.gz -I ${PT}-chr18:1-50000000-f1r2.tar.gz -I ${PT}-chr18:50000001-80373285-f1r2.tar.gz -I ${PT}-chr19:1-50000000-f1r2.tar.gz -I ${PT}-chr19:50000001-58617616-f1r2.tar.gz -I ${PT}-chr20:1-50000000-f1r2.tar.gz -I ${PT}-chr20:50000001-64444167-f1r2.tar.gz -I ${PT}-chr21:1-46709983-f1r2.tar.gz -I ${PT}-chr22:1-50000000-f1r2.tar.gz -I ${PT}-chr22:50000001-50818468-f1r2.tar.gz \
#  -O ${PT}-read-orientation-model.tar.gz
    
# # # Merge all the stats files together too #10 seconds
#   gatk MergeMutectStats -stats ${PT}.chr10:100000001-133797422.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr10:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr10:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr1:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr11:100000001-135086622.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr11:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr1:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr11:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr1:150000001-200000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr1:200000001-248956422.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr12:100000001-133275309.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr12:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr12:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr13:100000001-114364328.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr13:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr13:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr14:100000001-107043718.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr14:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr14:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr1:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr15:100000001-101991189.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr15:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr15:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr16:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr16:50000001-90338345.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr17:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr17:50000001-83257441.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr18:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr18:50000001-80373285.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr19:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr19:50000001-58617616.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr20:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr20:50000001-64444167.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr2:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr21:1-46709983.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr2:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr2:150000001-200000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr2:200000001-242193529.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr22:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr22:50000001-50818468.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr2:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr3:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr3:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr3:150000001-198295559.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr3:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr4:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr4:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr4:150000001-190214555.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr4:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr5:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr5:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr5:150000001-181538259.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr5:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr6:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr6:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr6:150000001-170805979.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr6:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr7:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr7:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr7:150000001-159345973.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr7:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr8:100000001-145138636.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr8:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr8:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr9:100000001-138394717.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr9:1-50000000.mutect2.somatic.unfiltered.vcf.gz.stats -stats ${PT}.chr9:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz.stats \
#          -O ${PT}.merged.stats

# # # Gather VCFs and create and index file # 3 minutes
#  gatk GatherVcfs -I ${PT}.chr1:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr1:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr1:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr1:150000001-200000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr1:200000001-248956422.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr2:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr2:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr2:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr2:150000001-200000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr2:200000001-242193529.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr3:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr3:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr3:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr3:150000001-198295559.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr4:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr4:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr4:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr4:150000001-190214555.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr5:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr5:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr5:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr5:150000001-181538259.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr6:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr6:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr6:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr6:150000001-170805979.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr7:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr7:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr7:100000001-150000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr7:150000001-159345973.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr8:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr8:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr8:100000001-145138636.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr9:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr9:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr9:100000001-138394717.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr10:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr10:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr10:100000001-133797422.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr11:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr11:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr11:100000001-135086622.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr12:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr12:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr12:100000001-133275309.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr13:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr13:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr13:100000001-114364328.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr14:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr14:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr14:100000001-107043718.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr15:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr15:50000001-100000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr15:100000001-101991189.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr16:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr16:50000001-90338345.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr17:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr17:50000001-83257441.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr18:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr18:50000001-80373285.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr19:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr19:50000001-58617616.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr20:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr20:50000001-64444167.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr21:1-46709983.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr22:1-50000000.mutect2.somatic.unfiltered.vcf.gz -I ${PT}.chr22:50000001-50818468.mutect2.somatic.unfiltered.vcf.gz \
#  -RI true \
#  -O ${PT}.chr.Combined.mutect2.somatic.unfiltered.vcf.gz

# # # # Make and Index (.tbi)
#  gatk IndexFeatureFile -I ${PT}.chr.Combined.mutect2.somatic.unfiltered.vcf.gz

# #####################################################################################
# # bl0cklist from ENCODE: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/metadata/GRCh38.metadata.controls.txt.gz
gatk FilterMutectCalls \
   -R /home/ahorning/DNAseq/Reference/Homo_sapiens_assembly38.fasta \
   -V ${PT}.chr.Combined.mutect2.somatic.unfiltered.vcf.gz \
   -stats ${PT}.merged.stats \
   --orientation-bias-artifact-priors ${PT}-read-orientation-model.tar.gz \
   --max-events-in-region 30 \
   --max-alt-allele-count 2 \
   -XL /home/ahorning/DNAseq/Reference/hg38-blacklist.v2.bed \
   -XL /home/ahorning/DNAseq/Reference/ENCFF356LFX.bed \
   -O ${PT}.chr.Combined.mutect2.somatic.filtered.vcf.gz

# #Consider  Not run yet
# # 1 minute
gatk LeftAlignAndTrimVariants \
   -R /home/ahorning/DNAseq/Reference/Homo_sapiens_assembly38.fasta \
   -V ${PT}.chr.Combined.mutect2.somatic.filtered.vcf.gz \
   -O ${PT}.chr.Combined.mutect2.somatic.filtered.left.vcf.gz \
   --split-multi-allelics

# #started at 11:02AM - ~12:10
module load annovar
table_annovar.pl ${PT}.chr.Combined.mutect2.somatic.filtered.left.vcf.gz /scg/apps/software/annovar/2020-06-07/humandb/ \
                   -buildver hg38 \
                   -out ${PT}.chr.Combined.mutect2.somatic.filtered.left.anno.vcf.gz \
                   -remove \
                   -protocol refGene,clinvar_20170905,cosmic70,dbnsfp33a,cytoBand,avsnp150,gnomad_genome \
                   -operation gx,f,f,f,r,f,f \
                   -nastring .  \
                   -polish  \
                   -xref /scg/apps/software/annovar/2020-06-07/example/gene_xref.txt \
                   -vcfinput 

# #  gzip and Index output
module load htslib
bgzip -c ${PT}.chr.Combined.mutect2.somatic.filtered.left.anno.vcf.gz.hg38_multianno.vcf > ${PT}.chr.Combined.mutect2.somatic.filtered.left.anno.vcf.gz.hg38_multianno.vcf.gz
gatk IndexFeatureFile -I ${PT}.chr.Combined.mutect2.somatic.filtered.left.anno.vcf.gz.hg38_multianno.vcf.gz

# # Remove some large intermediate files
rm ${PT}.chr.Combined.mutect2.somatic.filtered.left.vcf.gz
rm ${PT}.chr.Combined.mutect2.somatic.filtered.left.vcf.gz.tbi
rm ${PT}.chr.Combined.mutect2.somatic.filtered.vcf.gz.tbi
rm ${PT}.chr.Combined.mutect2.somatic.filtered.vcf.gz
rm ${PT}.chr.Combined.mutect2.somatic.filtered.left.anno.vcf.gz.hg38_multianno.vcf
rm ${PT}.chr.Combined.mutect2.somatic.filtered.left.anno.vcf.gz.hg38_multianno.txt
rm ${PT}.chr.Combined.mutect2.somatic.filtered.left.anno.vcf.gz.avinput


# # Turn into more readable Table format
gatk VariantsToTable \
      -V ${PT}.chr.Combined.mutect2.somatic.filtered.left.anno.vcf.gz.hg38_multianno.vcf.gz \
      -F CHROM -F POS -F avsnp150 \
      -F REF -F ALT -F FILTER \
      -F TYPE \
      -F TRANSITION \
      -F Func.refGene \
      -F Gene.refGene \
      -F GeneDetail.refGene \
      -F ExonicFunc.refGene \
      -F AAChange.refGene \
      -F CLINSIG \
      -F cosmic70 \
      -F CADD_phred \
      -F GERP++_RS \
      -F SIFT_score \
      -F Polyphen2_HVAR_pred \
      -F cytoBand \
      -F gnomAD_genome_ALL \
      -GF AD \
      -GF DP \
      -O ${PT}.chr.Combined.mutect2.somatic.filtered.left.anno.table.txt
      
# #       --show-filtered \