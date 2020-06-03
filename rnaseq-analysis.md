RNAseq analysis
================
Eduardo González-Grandío
June 2nd, 2020

Data preprocessing (done remotely in usegalaxy.org)

1- Trim low quality bases (Trimmomatic)

2- Align paired-end reads (HISAT2)

3- Remove reads that do not map, map to several regions or are unpaired
(Filter BAM)

4- Count reads per annotated gene (FeatureCounts): outputs
“Feature\_count.txt” and “feature\_length.tabular” files

5- Compare expression between samples (edgeR):

Based on
<https://ucdavis-bioinformatics-training.github.io/2019_March_UCSF_mRNAseq_Workshop/>

Load required libraries:

``` r
library(edgeR)
library(biomaRt)
library(data.table)
library(org.At.tair.db)
library(dendextend)
library(WebGestaltR)
library(ComplexHeatmap)
library(pals)
library(ggsci)
library(tidyverse)
library(splitstackshape)
library(stringr)
library(fgsea)
library(gridtext)
```
