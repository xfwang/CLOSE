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

cd /sonas-hs/powers/hpc/home/xuwang/reference/hg19
for C in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M
do
	cat chr${C}.fa >> hg19.fa
done
#creat .dict
java -jar ~/software/picard-tools-1.119/CreateSequenceDictionary.jar R= hg19.fa O= hg19.dict
#creat .fai
samtools faidx hg19.fa

date
