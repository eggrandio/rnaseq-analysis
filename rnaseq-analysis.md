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

5- Compare expression between samples (DESeq2):

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
library(magrittr)
library(splitstackshape)
library(stringr)
library(fgsea)
library(gridtext)
```

Calculate DGEList object (d) and calculate RPKM (outputs “RPKM.txt”
file):

``` r
library(edgeR)
input="./data/Feature_count.txt"
feature.length="./data/Feature_length.txt"

counts = read.delim(input, row.names = 1) #read output from FeatureCounts
head(counts)
```

    ##            NT1  NT2   NT3  NT4  NT5   WT1  WT2  WT3  WT4   WT5  DT1  DT2  DT3
    ## AT1G01010   48   51    97   24   55   178   81   74   57   122  376  336  293
    ## AT1G01020  204  230   172  180  227   147  125  122  132   188  101   90   76
    ## AT1G01030   91   89    59   79  115    87   85   50   42    58   28   46   34
    ## AT1G01040 2897 2954  2915 2113 2762  1942 1630 1263 1679  2171 1329 1202  977
    ## AT1G01050 1011 1030  1140 1018 1137   914  703  718  869  1270  867  913  607
    ## AT1G01060 8784 9078 10806 8453 9749 10315 6535 7330 8135 10937 5318 6376 4365
    ##            DT4  DT5  RT1  RT2  RT3  RT4  RT5
    ## AT1G01010  314  342  153   79   78   98   71
    ## AT1G01020  124  120  128  124  111   78  182
    ## AT1G01030   60   59   91   44   51   47   36
    ## AT1G01040 1484 1365 1246 1329 1132 1204 1809
    ## AT1G01050  837  868  583  646  626  596 1002
    ## AT1G01060 6398 6445 5707 6205 5307 4789 8054

``` r
GeneLength = (read.delim(feature.length, row.names = 1)) #read file with length of genes to calculate RPKM, Make sure they are in the same order as the counts file

snames = colnames(counts)
treatment = as.factor(substr(snames, 1, 2))
treatment = factor(treatment, levels=c("NT","WT","RT","DT"))

#Create Differential Gene Expression List Object (DGEList) object
d0 = DGEList(counts = counts, 
             genes = data.frame(Length=GeneLength),
             group = treatment)

d0 = calcNormFactors(d0) #Calculate normalization factors based on library size

RPKM = rpkm(d0) #Calculate reads per kilobase per million reads using normalization factors

RPKM_output = RPKM %>% as.data.frame %>% rownames_to_column("GeneID")
  
write.table(RPKM_output, file = paste0("./output/RPKM.txt"), row.names=F, sep="\t", quote=F)

d0 = d0[filterByExpr(d0), ,keep.lib.sizes=FALSE] #remove genes with low counts and recalculate library size

d0 = calcNormFactors(d0) #Re-calculate normalization factors after removing genes with low amount of reads
```
