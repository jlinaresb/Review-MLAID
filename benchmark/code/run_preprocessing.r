source('~/git/Review-MLAID/benchmark/code/preprocessing_functions.r')

# Transcriptomics
# ===
setwd("~/git/Review-MLAID/benchmark/extdata/single-omic/transcriptomic/")

datapath = list.files(pattern = 'GSE')
for (d in seq_along(datapath)) {
  data = match.condition(datapath = datapath[d], omic = 'transcriptomic')
  res = levene(data)
  outpath = 'filtered/'
  filename = gsub('.tsv', '.rds', datapath[d])
  saveRDS(res, file = paste0(outpath, filename))
}


# Epigenomics
# ===
setwd("~/git/Review-MLAID/benchmark/extdata/single-omic/epigenomic/")

datapath = list.files(pattern = 'GSE')
for (d in seq_along(datapath)) {
  data = match.condition(datapath = datapath[d], omic = 'epigenomic')
  
  y = data$Condition
  x = subset(data, select = -c(Condition))
  x = as.matrix(x)
  
  data = na.delete(x)
  data = remove.constant(data)
  data = lumi::beta2m(data)
  res = MostVariables(data, ntopGenes = 1000)
  res = as.data.frame(res)
  res$Condition = y
  
  outpath = 'filtered/'
  filename = gsub('.tsv', '.rds', datapath[d])
  saveRDS(res, file = paste0(outpath, filename))
}
