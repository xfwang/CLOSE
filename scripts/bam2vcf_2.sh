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
sample_name="yurolTM.reorder"

homedir="/sonas-hs/powers/hpc/home/xuwang/work_dir3"    #working directory
bamdir="/sonas-hs/powers/hpc/home/xuwang/YUROL"

ref_path="/sonas-hs/powers/hpc/home/xuwang/reference/hg19"
snpEff_path="/sonas-hs/powers/hpc/home/xuwang/software/snpEff"
GATK_path="/sonas-hs/powers/hpc/home/xuwang/software"
picard_path="/sonas-hs/powers/hpc/home/xuwang/software/picard-tools-1.119"

#sort, clean and index

cd $homedir
freebayes -f $ref_path/hg19.fa ${sample_name}.bam | vcffilter -f "QUAL > 30" > $homedir/${sample_name}.freebayes.vcf

cd $homedir; cat ${sample_name}.freebayes.vcf | java -Xmx4G -Xms4G -jar $snpEff_path/snpEff.jar eff -motif -nextProt -lof -v -c $snpEff_path/snpEff.config -hgvs -i vcf -o gatk GRCh37.71 > ${sample_name}.freebayes.ann.3.vcf 2> snpEff.3.log

cd $homedir; java -Xmx3G -Xms3G -jar $GATK_path/GenomeAnalysisTK.jar -R $ref_path/hg19.fa -T VariantAnnotator -I ${sample_name}.bam -o ${sample_name}.var_ann.3.vcf --variant ${sample_name}.freebayes.ann.3.vcf -L ${sample_name}.freebayes.ann.3.vcf -A Coverage -A QualByDepth -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A MappingQualityZero -A LowMQ -A RMSMappingQuality -A BaseQualityRankSumTest -rf BadCigar > gatk_ann.3.log 2>&1

cat ${sample_name}.var_ann.3.vcf | vcffilter -f "DP > 7" | python ~/scripts/freebayes_filter.py 40 40 30 4 0 > ${sample_name}.annotated.filtered3.vcf 2> fb_filter3.log


date
