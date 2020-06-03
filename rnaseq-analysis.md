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
```

    ## Loading required package: limma

``` r
library(biomaRt)
library(data.table)
library(org.At.tair.db)
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colMeans,
    ##     colnames, colSums, dirname, do.call, duplicated, eval, evalq,
    ##     Filter, Find, get, grep, grepl, intersect, is.unsorted, lapply,
    ##     lengths, Map, mapply, match, mget, order, paste, pmax, pmax.int,
    ##     pmin, pmin.int, Position, rank, rbind, Reduce, rowMeans, rownames,
    ##     rowSums, sapply, setdiff, sort, table, tapply, union, unique,
    ##     unsplit, which, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     first, second

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:data.table':
    ## 
    ##     shift

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## 

``` r
library(dendextend)
```

    ## 
    ## ---------------------
    ## Welcome to dendextend version 1.13.4
    ## Type citation('dendextend') for how to cite the package.
    ## 
    ## Type browseVignettes(package = 'dendextend') for the package vignette.
    ## The github page is: https://github.com/talgalili/dendextend/
    ## 
    ## Suggestions and bug-reports can be submitted at: https://github.com/talgalili/dendextend/issues
    ## Or contact: <tal.galili@gmail.com>
    ## 
    ##  To suppress this message use:  suppressPackageStartupMessages(library(dendextend))
    ## ---------------------

    ## 
    ## Attaching package: 'dendextend'

    ## The following object is masked from 'package:data.table':
    ## 
    ##     set

    ## The following object is masked from 'package:stats':
    ## 
    ##     cutree

``` r
library(WebGestaltR)
```

    ## ******************************************

    ## *                                        *

    ## *          Welcome to WebGestaltR !      *

    ## *                                        *

    ## ******************************************

``` r
library(ComplexHeatmap)
```

    ## Loading required package: grid

    ## ========================================
    ## ComplexHeatmap version 2.3.4
    ## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    ## Github page: https://github.com/jokergoo/ComplexHeatmap
    ## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    ## 
    ## If you use it in published research, please cite:
    ## Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    ##   genomic data. Bioinformatics 2016.
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(ComplexHeatmap))
    ## ========================================

``` r
library(pals)
library(ggsci)
library(tidyverse)
```

    ## -- Attaching packages -------------------------------------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.0     v purrr   0.3.4
    ## v tibble  3.0.1     v dplyr   0.8.5
    ## v tidyr   1.0.2     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts ----------------------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::between()    masks data.table::between()
    ## x dplyr::collapse()   masks IRanges::collapse()
    ## x dplyr::combine()    masks Biobase::combine(), BiocGenerics::combine()
    ## x dplyr::desc()       masks IRanges::desc()
    ## x tidyr::expand()     masks S4Vectors::expand()
    ## x dplyr::filter()     masks stats::filter()
    ## x dplyr::first()      masks S4Vectors::first(), data.table::first()
    ## x dplyr::lag()        masks stats::lag()
    ## x dplyr::last()       masks data.table::last()
    ## x ggplot2::Position() masks BiocGenerics::Position(), base::Position()
    ## x purrr::reduce()     masks IRanges::reduce()
    ## x dplyr::rename()     masks S4Vectors::rename()
    ## x dplyr::select()     masks AnnotationDbi::select(), biomaRt::select()
    ## x dplyr::slice()      masks IRanges::slice()
    ## x purrr::transpose()  masks data.table::transpose()

``` r
library(splitstackshape)
library(stringr)
library(fgsea)
```

    ## Loading required package: Rcpp

``` r
library(gridtext)
```
