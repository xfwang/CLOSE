CLOSE-R 
====

**Title**  A toolkit for **C**NA/**LO**H analysis with **Se**quencing data in R

**Data**  2015-09-28

**Author**  Xuefeng Wang

**Maintainer**  Xuefeng Wang <xuefeng.wang@stonybrook.edu>

**Description**   CLOSE-R is a toolkit for **C**NA and **LO**H analysis (as well as **CLO**nality analysis) with **SE**quencing data implemented in R. Current pipeline majorly facilitates the analysis on paired tumor and normal samples. This pipeline conssits of three major compartments: (1) ASCN (allel-sepcific copy number) estimation using a distance-based Chinese Restaurant Process; (3) global purity and ploidy estimation; (3) Genome-wide ASCN visulization

**Depends** R (>= 3.22) DPpackage,grid,ggplot2,VariantAnnotation

_____________________________________________________________________________________________________________

``CLOSER          Main function of CLOSE-R``
_____________________________________________________________________________________________________________


**Description**

  * Perform CLOSE-R analysis starting from a segment-level input
 
**Usage**

  * CLOSER(Input, codeDir, output, sampleName)
 
**Arguments**

  * Input
    * input data matrix, LAF (lesser allale frequency) and LRR (log2 of the read depth ratio tumor/normal) at segment-level. This matrix is required to have five columns: chromosome, start position, end position of the segments, LAF (summerized LAF of the segment, usually mean or median), and LRR (summerized LRR of the segment, usually mean or median). See [example.input](https://github.com/xfwang/CLOSE/blob/master/CLOSE-R/example.input.txt) for example 
 
  * codeDir
    * direcotry where the sub-function script subFunc.R is saved
 
  * output
    * desired location of output files
 
  * sampleName
    * output prefix; all output files created by CLOSER will have this prefix (e.g. .CNstatus.txt, .plotCNR.pdf, etc.). If this option is not provided the default output prefix being used is: "closer"
 

