GO_hm = function(genelist,ORA_GO_genes){
  
  GO_cats = ORA_GO_genes$geneSet #select enriched GOs only
  GO_heatmap = as.data.frame(matrix(ncol = length(GO_cats), nrow = length(genelist), dimnames=list(names(genelist),GO_cats)))
  #data_GO = mapIds(org.At.tair.db, names(genelist), "GOALL", "TAIR",multiVals = "list") #get all GOs for each gene
  genes = ORA_GO_genes$userId
  
  data_GO = vector(mode = "list", length = length(GO_cats))
  names(data_GO) = GO_cats
  
  for(i in seq(data_GO)){
    data_GO[i] = strsplit(genes[i],";")
  }
  
  for(i in seq(colnames(GO_heatmap))){
    GO = colnames(GO_heatmap)[i]
    for(j in 1:length(rownames(GO_heatmap))){
      GO_heatmap[j,i] = ifelse(rownames(GO_heatmap)[j] %in% unlist(data_GO[GO]),1,0)
    }
  }
  colnames(GO_heatmap) = ORA_GO_genes$description
  return(GO_heatmap)
}