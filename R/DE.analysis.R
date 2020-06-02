DE.analysis = function(input="./data/Feature_count.txt",feature.lenght="./data/feature_length.tabular"){
  library(edgeR)
  counts = read.delim(input, row.names = 1) #read output from FeatureCounts
  GeneLength = (read.delim(feature.lenght, row.names = 1)) #read file with length of genes to calculate RPKM
  
  #Create Differential Gene Expression List Object (DGEList) object
  d0 = DGEList(counts=counts,genes=data.frame(Length=GeneLength)) 
  d0 = calcNormFactors(d0) #Calculate normalization factors based on library size
  RPKM = rpkm(d0) #Calculate reads per kilobase per million reads using normalization factors
  
  #Process data and generate PCA plot of samples based on the transcriptomic profile
  cutoff = 1
  drop = which(apply(cpm(d0), 1, max) < cutoff)
  d = d0[-drop,] 
  snames = colnames(counts)
  treatment = as.factor(substr(snames, 1, 2))
  postscript(file="./output/MDSplot.ps")
  plotMDS(d, col = as.numeric(treatment))
  suppressMessages(dev.off())
  
}