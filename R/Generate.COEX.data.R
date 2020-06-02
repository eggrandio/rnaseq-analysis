library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # to set parameters of dendrogram and to compare two dendrograms
library(NbClust)
library(biomaRt)
library(data.table)
library(org.At.tair.db)

ensembl = useMart("plants_mart",host="plants.ensembl.org")
ensembl = useDataset("athaliana_eg_gene",mart=ensembl)

geneIDs = getBM(attributes = c("entrezgene_id", "ensembl_gene_id"), mart = ensembl)

#Generation of Coexpression dataset
#Data downloaded from http://atted.jp/download.shtml

setwd("D:/Box Sync/Landry lab/Landry lab (eggrandio@berkeley.edu)/RNAseq/Ath rB.v17 12.G22760 S2120.combat.mrgeo.d_1/Ath-rB.v17-12.G22760-S2120.combat.mrgeo.d/")

#Make sure only files from ATTED-II are in the folder:

temp = list.files()
temp = setNames(object = temp, nm = temp) #name each element in the list with the file name
temp = temp[order(as.numeric(names(temp)))]

V1 = fread(temp[1], head=F) %>% arrange(V1) %>% dplyr::select(-V2)

data = lapply(temp, function(x) {
  fread(x, head=F) %>% arrange(V1) %>% dplyr::select(-V1)
})
test = cbind(V1, do.call(cbind, data))


row.names(test) = test[,1]
test[,1] = NULL
colnames(test) = temp
test = test[order(as.numeric(row.names(test))), order(as.numeric(names(test)))]

genes_coexRNA = length(colnames(test)) # number of genes with coexpression data in RNA dataset
coexRNA = test/genes_coexRNA #range-standardized MR index
logit_coexRNA = log(coexRNA/(1-coexRNA)) #logit(MR) converted dataset (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5914358/ )


#saveRDS(test, "D:/Box Sync/Landry lab/Landry lab (eggrandio@berkeley.edu)/RNAseq/COEX_table_Entrezid_mARRAY.Rds")

saveRDS(logit_coexRNA, "D:/Box Sync/Landry lab/Landry lab (eggrandio@berkeley.edu)/RNAseq/COEX_table_Entrezid_RNAseq_logit.Rds")
