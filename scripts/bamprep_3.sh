#!/bin/bash
# specify BASH shell
#$ -S /bin/bash
# pass environment variables to job, e.g. LD_LIBRARY_PATH
#$ -v LD_LIBRARY_PATH
# run job in the current working directory where qsub is executed from
#$ -cwd
#  specify that the job requires 16GB of memory
#$ -l m_mem_free=16G
 
# run commands and application
pwd
date

#bam to vcf (bnb)
sample_name="yurolTM"
SAMPLE="yurolTM"

homedir="/sonas-hs/powers/hpc/home/xuwang/work_dir3"    #working directory
bamdir="/sonas-hs/powers/hpc/home/xuwang/YUROL"

ref_path="/sonas-hs/powers/hpc/home/xuwang/reference/hg19"
snpEff_path="/sonas-hs/powers/hpc/home/xuwang/software/snpEff"
GATK_path="/sonas-hs/powers/hpc/home/xuwang/software"
picard_path="/sonas-hs/powers/hpc/home/xuwang/software/picard-tools-1.119"

cd $homedir
#sort, clean and index
#samtools sort $bamdir/${sample_name}.bam $bamdir/${sample_name}.sorted
#samtools view -F 4 -b $bamdir/${sample_name}.sorted.bam  > $bamdir/${sample_name}.clean.bam
#samtools index $bamdir/${sample_name}.clean.bam
#rm $bamdir/${sample_name}.sorted.bam

samtools index $bamdir/${sample_name}.bam
samtools view -bh $bamdir/${sample_name}.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM > ${sample_name}.clean.bam
samtools sort ${sample_name}.clean.bam ${sample_name}.clean.sorted
samtools view -h ${SAMPLE}.clean.sorted.bam | sed '/gl00/d' | sed '/hap/d'> ${SAMPLE}.clean.sorted.sam
samtools view -bS ${SAMPLE}.clean.sorted.sam > ${SAMPLE}.clean2.bam

##order according to the reference genome
java -jar ${picard_path}/ReorderSam.jar INPUT=${SAMPLE}.clean2.bam OUTPUT=${SAMPLE}.reorder.bam REFERENCE=${ref_path}/hg19.fa VALIDATION_STRINGENCY=SILENT
##index
samtools index ${SAMPLE}.reorder.bam

date
