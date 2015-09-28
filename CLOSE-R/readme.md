CLOSE-R 
====

**Title**  A toolkit for **C**NA/**LO**H analysis with **Se**quencing data in R

**Data**  2015-09-28

**Author**  Xuefeng Wang

**Maintainer**  Xuefeng Wang <xuefeng.wang@stonybrook.edu>

**Description**   CLOSE-R is a toolkit for **C**NA and **LO**H analysis (as well as **CLO**nality analysis) with **SE**quencing data implemented in R. Current pipeline majorly facilitates the analysis on paired tumor and normal samples. This pipeline conssits of three major compartments: (1) ASCN (allel-sepcific copy number) estimation using model-free approach (distance-based Chinese Restaurant Process) or model-based approach (MAP, Maximum a posteriori); (3) global purity and ploidy estimation; (3) Genome-wide ASCN visulization

**Depends** R (>= 3.22) DPpackage,grid,ggplot2,VariantAnnotation
<br><br><br>


## Main function of CLOSE-R

#### [CRP.R](https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/segCN.R)

**Description**

  * Estimate segments copy number status using model-free CRP approach (Chinese Restaurant Process)
 
**Usage**

  * CRP(Input, codeDir, outDir, sampleName)
 
**Arguments**

  * **Input**
    * input data matrix, LAF (lesser allale frequency) and LRR (log2 of the read depth ratio tumor/normal) at segment-level. This matrix is required to have five columns: chromosome, start position, end position of the segments, LAF (summerized LAF of the segment, usually mean or median), and LRR (summerized LRR of the segment, usually mean or median). See [example.input](https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/example.input.txt) for example 

| chromosome |     start |       end |         LAF |          LRR |
|-----------:|----------:|----------:|------------:|-------------:|
|          1 |     69223 |  4834523 |  0.445030335 |  -0.038606278 |
|          1 |  6100667 | 14142993 | 0.4375  | 0.042693104|
|          1 | 15287258 | 29652142 |  0.436363636 | 0.015784217 |

  * **codeDir**
    * direcotry where the sub-function script subFunc.R is saved (e.g., /home/CLOSER_code/)
 
  * **outDir**
    * desired location of output files (e.g., /home/CLOSER_output/)
 
  * **sampleName**
    * output prefix; all output files created by `<CRP>` will have this prefix (e.g., .CNstatus.txt, .plotCNR.pdf, etc.). If this option is not provided the default output prefix being used is: "closer"
 
**Output**
 * `<sampleName>`.CNstatus.txt (See [example.CNstatus.txt](https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/example.CNstatus.txt)) containing the following columns:
    * chromosome: chromosome of segments
    * start: start location of segments
    * end：end location of segments
    * LAF: summarized LAF (lesser allele frequency) of segments
    * LRR: summarized LRR of segments
    * minCNR: estimated minor allele copy number ratio of segments
    * majCNR: estimated major allele copy number ratio of segments
    * cluster: index of cluster this segment belongs to (based on a distance-based Chinese Restaurant Process)
    * status_cluster: copy number status estimated based on clusters (Normal, Gain, Loss, CN Neutral LOH, and Ambiguous)
    * status_seg: copy number status estimated based on segments (Normal, Gain, Loss, CN Neutral LOH, and Ambiguous)
 <br>

 * `<sampleName>`.plotCNR.pdf plots the minCNR/majCNR estimated by segCN 
    * **Assessment of copy number status**
    ![majCNR vs. minCNR] (https://github.com/xfwang/CLOSE/blob/master/instr/image/majCNR.vs.minCNR.png)
     Assessment of copy number status of one sample based on distance based modified Chinese restaurant
process (CRP). The relative ASCN (allele specific copy number) estimates can be calculated based on
LRR and BAF. The center of each circle in the plot depictes the minCN and majCN ratios (or relative
ASCN) estiamted from each genomic segment. The size of the circle indicates the length of one segment.
The colors of the circles indicate the grouping results from the Chinese restaurant process (CRP). As
shown in the figure, CRP enables the partition of the genome-wide CNA profile into blocks corresponding
to different CNA status. The main advantage of the algorithm is that it allows unknown number of clusters
without need to model number of clusters directly. The cluster (pink cluster) that is closest to the baseline
point (minCN=1 and majCN=1) corresponds to the normal regions, and other cluster status can be
inferred accordingly, as described in the supplementary method.
<br>

    * **Purity estimation**
    ![purity] (https://github.com/xfwang/CLOSE/blob/master/instr/image/purity.png)
    Minor allelic CN is fitted by Dirchlet Process
<br><br>

#### [MAP.R] (https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/MAP.R)

**Description**
  * Estimate global purity and ploidy using MAP (Maximum a posteriori)
<br><br> 

#### [CLOSE-R.R] (https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/CLOSE-R.R)

**Description**
  * Estimate local copy number given ploidy and purity (based on segmental ASCN coordinates)

**Output**
  * Paramter trajectory plots (canonical lines) of different copy number status along with purity
![canonical] (https://github.com/xfwang/CLOSE/blob/master/instr/image/canonical.png)
<br><br><br>


## Copy number visulization
### Copy number status by chromosome
##### Usage
  * [plotCNstatus.chr(CNstatus, BAF, LRR, sampleName)] (https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/subFunc.R)

##### Arguments
  * CNstatus：`<sampleName>`.CNstatus.txt generated by `<CRP>`  (See [example.CNstatus.txt](https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/example.CNstatus.txt))
  * BAF: BAF (B allele freqency) values for all SNP sites in tumor sample; a matrix with three columns (chr, position, BAF); user may get BAF information of SNP sites from normal and tumor .vcf files using VCFprep.R inside of subFunc.R script
  * LRR: LRR values for all SNP sites; a matrix with three columns (chr, position, BAF); user may get LRR information of SNP sites from normal and tumor .vcf files using VCFprep.R inside of subFunc.R script
  * sampleName: output prefix

#####Output
![chr] (https://github.com/xfwang/CLOSE/blob/master/instr/image/visulization.chr.png)
<br><br>

### Global Copy number status
#####Usage
  * [plotCNstatus.WG(CNstatus, BAF, LRR, sampleName)] (https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/subFunc.R)

#####Arguments
  * see `<plotCNstatus.chr>`

#####Output
![chr] (https://github.com/xfwang/CLOSE/blob/master/instr/image/visulization.WG.png)
<br><br>

### Comparison of copy number calls from WES and SNP array 

#####Usage
  * [compareToArray(CNstatus, SNParray, sampleName)] (https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/subFunc.R)

#####Arguments
  * CNstatus：`<sampleName>`.CNstatus.txt generated by `<CRP>`
  * SNParray: copy number calls from WES and SNP arrays; 4 columns (1)chr, (2)start of region, (3)end of region, (4)copy number estimates 

#####Output
![TCGA] (https://github.com/xfwang/CLOSE/blob/master/instr/image/TCGA.png)
Comparison of global copy number calls from WES and SNP array (Affeymetrix SNP6). The segmented copy number estimates (Segment_mean) based on SNP6 are downloaded from broad firehose website (http://gdac.broadinstitute.org/runs/stddata__2015_04_02).
