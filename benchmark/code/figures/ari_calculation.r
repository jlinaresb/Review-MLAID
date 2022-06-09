setwd('~/git/Review-MLAID/benchmark/results/')

files = list.files()
time = list()
l = list()
res = list()
# f = 1
for (f in seq_along(files)) {
  data = readRDS(files[f])
  data = data[-c(3, 4, 8)]
  
  
  for (i in seq_along(data)) {
    l[[i]] = data.frame(
      ids = names(data[[i]]$labels),
      labels = data[[i]]$labels
    )
    colnames(l[[i]])[2] = names(data)[i]
  }
  
  x = as.data.frame(l)
  
  res[[f]] = x[,-grep('id', colnames(x))]
  
} 

names(res) = files

# ARI calculation
# ===

require(mclust)
require(corrplot)
require(viridis)
kk = res$metagenomic_morgan.rds
corrplots = list()
for (i in seq_along(res)) {
  kk = res[[i]]
  name = names(res)[i]
  
  n <- ncol(kk)
  algs = colnames(kk)
  
  ari <- function(i, j, data) {adjustedRandIndex(data[,i], data[,j])}
  corp <- Vectorize(ari, vectorize.args = list("i","j"))
  ari_res = outer(1:n, 1:n, corp, data=kk)
  
  rownames(ari_res) = algs
  colnames(ari_res) = algs
  
  corrplots[[i]] = corrplot(ari_res,
                            type = 'lower',
                            col = viridis(10),
                            tl.col = 'black',
                            title = name)
}
names(corrplots) = names(res)
