setwd('~/git/Review-MLAID/benchmark/results/preciseads/')

# Single-omic clustering
# ===
files = list.files(pattern = 'single')
ari.s = list()
# i = 1
for (i in seq_along(files)) {
  res = readRDS(files[i])
  
  labels = list(
    Consensus.H  = res[[1]]$labels, 
    Consensus.KM = res[[2]]$labels,
    SOM          = res[[3]]$labels,
    MClust       = res[[4]]$labels,
    NMF          = res[[5]]$labels
  )
  
  ari = data.frame(
    hierarchical = labels$Consensus.H,
    kmeans = labels$Consensus.KM,
    som = labels$SOM,
    mclust = labels$MClust,
    nmf = labels$NMF,
    row.names = names(labels[[1]])
  )
  
  # ARI calculation
  require(mclust)
  require(corrplot)
  
  n = ncol(ari)
  algs = colnames(ari)
  
  ari.fun = function(i, j, data) {adjustedRandIndex(data[, i], data[, j])}
  corp = Vectorize(ari.fun, vectorize.args = list('i', 'j'))
  ari.s[[i]] = outer(1:n, 1:n, corp, data = ari)
  rownames(ari.s[[i]]) = algs 
  colnames(ari.s[[i]]) = algs
  
  names(ari.s)[i] = files[i]
  
}


# Multi-omic clustering
# ===
files = list.files(pattern = 'multi')
files = files[-grep('coca', files)]
ari.m = list()
# i = 27
for (i in seq_along(files)) {
  res = readRDS(files[i])
  
  ari.m[[i]] = data.frame(
    labels = res$labels,
    ids = names(res$labels),
    algorithm = sapply(strsplit(files[i], '_'), '[[', 2),
    features = sapply(strsplit(files[i], '_'), '[[', 3),
    row.names = 1:length(res$labels)
  )
}
ari.m = data.table::rbindlist(ari.m)


head(ari.m)
table(ari.m$features)
# Split by experiment
# ==
experiment = unique(ari.m$features)

exp = list()
# i = 1
for (i in seq_along(experiment)) {
  byexp = ari.m[which(ari.m$features == experiment[i]),]
  
  exp[[i]] = reshape2::dcast(dat = byexp,
                             formula = ids ~ algorithm,
                             fun.aggregate = sum,
                             value.var = 'labels')
  
}
dim(exp[[i]])