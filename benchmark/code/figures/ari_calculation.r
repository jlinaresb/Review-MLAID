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
table(ari.m$algorithm)
# Split by experiment
# ==
experiment = unique(ari.m$features)

exp = list()
i = 1
for (i in seq_along(experiment)) {
  byexp = ari.m[which(ari.m$features == experiment[i]),]
  
  exp[[i]] = reshape2::dcast(dat = byexp,
                             formula = ids ~ algorithm,
                             fun.aggregate = sum,
                             value.var = 'labels')
  exp[[i]] = exp[[i]][,2:9]
  n = ncol(exp[[i]])
  algs = colnames(exp[[i]])
  
  ari.fun = function(i, j, data) {adjustedRandIndex(data[, i], data[, j])}
  corp = Vectorize(ari.fun, vectorize.args = list('i', 'j'))
  exp[[i]] = outer(1:n, 1:n, corp, data = exp[[i]])
  rownames(exp[[i]]) = algs 
  colnames(exp[[i]]) = algs
  
}
names(exp) = experiment
ari.m = exp


# Plotting
# ===
# ari.s
# ari.m

require(corrplot)
require(viridis)
# corrplot(ari.s[[f]],
#          method = 'square',
#          type = 'lower', is.corr = T,
#          col = rev(viridis(10)),
#          title = names(ari.m)[f],
#          tl.col = 'black')
# 
# names(ari.s)
f = 3
ggcorrplot::ggcorrplot(ari.s[[f]],
                       hc.order = F, 
                       type = 'lower',
                       show.diag = T,
                       lab = F,
                       legend.title = 'ARI',
                       method = 'square') +
  scale_fill_continuous(type = 'viridis') +
  ggtitle(gsub('.rds', '', names(ari.s)[f]))



# Probe

all = c(ari.s, ari.m)

plotting = list()
names = names(all)
for (i in seq_along(names)) {
  plotting[[i]] = all[[i]]
  plotting[[i]][upper.tri(plotting[[i]])] = NA
  plotting[[i]] = melt(plotting[[i]], na.rm = T)
  plotting[[i]]$experiment = gsub('.rds', '', names[i])
}
plotting = data.table::rbindlist(plotting)

ggplot(data = plotting, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_continuous(type = 'viridis') +
  theme_light() +
  facet_wrap(~experiment, scales = 'free') +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank())
  
  
