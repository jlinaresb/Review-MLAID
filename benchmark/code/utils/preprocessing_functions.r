# Match patients with condition
# ===
match.condition = function(datapath, omic = 'transcriptomic'){
  
  # Arguments:
  #  datapath
  #  omic: transcriptomic, epigenomic or metagenomic
  # 
  # Value:
  #  return data with condition column. Genes in columns and samples in rows. 
  
  if (omic == 'transcriptomic') {
    meta = read.delim2(
      'benchmark/extdata/single-omic/transcriptomic/metadata.tsv',
      header = T)
  } else if (omic == 'epigenomic'){
    meta = read.delim2(
      'benchmark/extdata/single-omic/epigenomic/metadata.tsv',
      header = T)
  } else if (omic == 'metagenomic'){
    meta = NULL
  }
  
  data = data.table::fread(datapath, header = T)
  data = as.data.frame(data)
  rownames(data) = data$gene
  data = data[, -1]
  
  pats = colnames(data)
  meta = meta[match(pats, meta$Sample), ]
  
  data = as.data.frame(t(data))
  data$Condition = meta$Condition
  
  print('done!')
  
  return(data)
  
}


# NA removing
# ===
na.delete = function(data){
  
  pre = ncol(data)
  data = data[, apply(data, 2, function(x) !any(is.na(x)))]
  post = ncol(data)
  
  print(paste0(pre-post, ' features were removed'))
  
  return(data)
  
}


# Delete constant features
# ===
remove.constant = function(data){
  
  pre = ncol(data)
  data = data[,apply(data, 2, var) >= 0.001]
  post = ncol(data)
  
  print(paste0(pre-post, ' features were removed'))
  return(data)
}



# PCA
# ===
pca_fun = function(data){
  
  require(PCAtools)
  fit = PCAtools::pca(data,
                      scale = T,
                      center = T,
                      removeVar = 0.1, 
                      transposed = T)
  # Select number of PCA?s
  elbow = findElbowPoint(fit$variance)
  
  data = as.matrix(fit$rotated[, 1:elbow])
  
  return(data)
  
}



# Most variables features
# ===
MostVariables = function(data, ntopGenes = 1000){
  
  mads = apply(data, 2, mad)
  data = data[, rev(order(mads))[1:ntopGenes]]
  
  print(paste0('The ', ntopGenes, ' most variable features were selected!'))
  
  return(data)  
}
 

levene = function(data, condition = 'Condition', alpha = 0.05){
  
  # Arguments 
  #   - data: matrix or data.frame with genes in columns and samples in rows. One
  #           columns must be condition.
  #   - condition: character string with the name of variable name of condition.
  #                 This variable codes control and sample patients.
  #   - alpha: level of significance
  # 
  # Return
  #     vector with padj of each gene 
  
  stopifnot(condition %in% colnames(data))
  
  # Convert group to factor
  data[, condition] = as.factor(data[, condition])
  
  # Select genes
  genes = colnames(data)
  genes = genes[- grep(condition, genes)]
  
  # Leven test (bonferroni posthoc)
  pvalues = list()
  for (i in seq_along(genes)) {
    
    pvalues[[i]] = car::leveneTest(data[,i] ~ data[, condition], data)$`Pr(>F)`[1]
    print(paste0('Levene test calculated for: ', genes[i], ' ',  i, ' / ', length(genes)))
    
  }
  names(pvalues) = genes
  pvalues = unlist(pvalues)
  
  # Post hoc test
  padj = p.adjust(pvalues, method = 'bonferroni')
  
  # Sort and select by significance
  padj.f = padj[which(padj < alpha)]
  
  if (length(padj.f) < 100) {
    padj.f = sort(padj)[1:100]
  }
  
  # Select these genes in data
  data = data[, match(names(padj.f), colnames(data))]
  
  return(data)
  
}

