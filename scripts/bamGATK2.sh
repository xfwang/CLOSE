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
SAMPLE="yurolN"

homedir="/sonas-hs/powers/hpc/home/xuwang/work_dir"    #working directory
bamdir="/sonas-hs/powers/hpc/home/xuwang/YUROL"

ref_path="/sonas-hs/powers/hpc/home/xuwang/reference/hg19_orc"
ref_fa="hg19.fa"
snpEff_path="/sonas-hs/powers/hpc/home/xuwang/software/snpEff"
GATK_path="/sonas-hs/powers/hpc/home/xuwang/software"
picard_path="/sonas-hs/powers/hpc/home/xuwang/software/picard-tools-1.119"

##bamGATK.sh
#cd $bamdir
#samtools index ${SAMPLE}.bam #index
#for C in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M
#        do
#        samtools view ${SAMPLE}.bam chr${C} -b > ${SAMPLE}.${C}.chr.bam
#        done
#samtools merge ${SAMPLE}.merged.bam  *.chr.bam
#rm *.chr.bam

#samtools view -H ${SAMPLE}.merged.bam | sed '/gl00/d' | sed '/hap/d'| samtools reheader - ${SAMPLE}.merged.bam > ${SAMPLE}.merged_new.bam

##bamGATK2.sh
cd $bamdir
rm ${SAMPLE}.merged.bam

##sort by coordiante
#java -jar ${picard_path}/SortSam.jar INPUT=${SAMPLE}.merged_new.bam OUTPUT=${SAMPLE}.sort.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
samtools sort ${SAMPLE}.merged_new.bam ${SAMPLE}.sort

##order according to the refernen genome
java -jar ${picard_path}/ReorderSam.jar INPUT=${SAMPLE}.sort.bam OUTPUT=${SAMPLE}.reorder.bam REFERENCE=${ref_path}/${ref_fa} VALIDATION_STRINGENCY=SILENT
##add index file
samtools index ${SAMPLE}.reorder.bam    #.-----> final output file




date
