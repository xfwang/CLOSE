CLOSE-R 
====

**Title**  A toolkit for **C**NA/**LO**H analysis with **Se**quencing data in R

**Data**  2015-09-28

**Author**  Xuefeng Wang

**Maintainer**  Xuefeng Wang <xuefeng.wang@stonybrook.edu>

**Description**   CLOSE-R is a toolkit for **C**NA and **LO**H analysis (as well as **CLO**nality analysis) with **SE**quencing data implemented in R. Current pipeline majorly facilitates the analysis on paired tumor and normal samples. This pipeline conssits of three major compartments: (1) ASCN (allel-sepcific copy number) estimation using a distance-based Chinese Restaurant Process; (3) global purity and ploidy estimation; (3) Genome-wide ASCN visulization

**Depends** R (>= 3.22) DPpackage,grid,ggplot2,VariantAnnotation
<br><br><br>


## Main function of CLOSE-R

**Description**

  * Perform CLOSE-R analysis starting from a segment-level input
 
**Usage**

  * CLOSER(Input, codeDir, outDir, sampleName)
 
**Arguments**

  * **Input**
    * input data matrix, LAF (lesser allale frequency) and LRR (log2 of the read depth ratio tumor/normal) at segment-level. This matrix is required to have five columns: chromosome, start position, end position of the segments, LAF (summerized LAF of the segment, usually mean or median), and LRR (summerized LRR of the segment, usually mean or median). See [example.input](https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/example.input.txt) for example 
 
  * **codeDir**
    * direcotry where the sub-function script subFunc.R is saved (e.g., /home/CLOSER_code/)
 
  * **outDir**
    * desired location of output files (e.g., /home/CLOSER_output/)
 
  * **sampleName**
    * output prefix; all output files created by CLOSER will have this prefix (e.g., .CNstatus.txt, .plotCNR.pdf, etc.). If this option is not provided the default output prefix being used is: "closer"
 
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
    * status_cluster: copy number status estimated based on clusters
    * status_seg: copy number status estimated based on segments
 <br>

 * `<sampleName>`.plotCNR.pdf plots the minCNR/majCNR estimated by CLOSE-R
    * Assessment of copy number status
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

    * purit estimation
    ![purity] (https://github.com/xfwang/CLOSE/blob/master/instr/image/purity.png)
    
