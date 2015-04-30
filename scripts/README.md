Following software need to be installed before using this pipeline: 
Freebayes; GATK; Picard Tools, SnpEff, falcon

**hg19_2.sh**: bash file for generating hg19 reference genome files

**bamprep_3.sh**: sort (by coordinates), clean and index bam files, and then reorder according to the reference genome to be readable by freebayes and GATK. 

**bam2vcf_2.sh**: file to convert BAM to VCF files, assuming bam file is well sorted and cleaned.

**ascn.r**: R-based pipeline to call ASCN from vcf files. change the input file and sample name in the first three lines.The outputs are the summary text file, the ascn results and the plots.

**freebayes_filter.py**:other supporting scripts need to put in appropriate directory, please check bam2vcf.sh

**bamGATK.sh** and **bamGATK_2.sh**: archive scripts, when bamprep_3.sh fails, one can use these two files sequentially. 

