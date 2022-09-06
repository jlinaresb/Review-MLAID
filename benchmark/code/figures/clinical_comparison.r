# Clinical comparisons
# ===
require(reshape2)
require(ggpubr)

setwd('~/git/Review-MLAID/benchmark/')

clinical = read.delim2('extdata/PRECISEADS/Metadata.tsv', 
                       header = T, sep = '\t', row.names = 1)
rownames(clinical) = make.names(rownames(clinical))
clinical[,c(11, 17:22)] = apply(clinical[,c(11, 17:22)], 2, function(x) as.numeric(x))

# Single omic
# ===
path = 'results/preciseads/'
files = list.files(path = path, pattern = 'single_')
compare = list()
# i = 1
for (i in seq_along(files)) {
  res = readRDS(paste0(path, files[i]))
  
  l = list(
    Consensus.H  = res[[1]]$labels, 
    Consensus.KM = res[[2]]$labels,
    SOM          = res[[3]]$labels,
    MClust       = res[[4]]$labels,
    NMF          = res[[5]]$labels
  )
  
  labels = data.frame(
    hierarchical = l$Consensus.H,
    kmeans = l$Consensus.KM,
    som = l$SOM,
    mclust = l$MClust,
    nmf = l$NMF,
    row.names = names(l[[1]])
  )
  
  info = clinical[match(rownames(labels), rownames(clinical)), 17:22]
  labels = cbind.data.frame(labels, info)
  
  compare[[i]] = labels
}
names(compare) = gsub('.rds', '', files)

toPlot = compare[[1]]

algs = names(toPlot)[1:5]
toPlot = melt(toPlot, 
              measure.vars = c('Monocytes', 'Neutrophils', 'BCells',
                          'NKCells', 'CD4TCells', 'CD8TCells'))

# ggboxplot(data = toPlot, 
#           x = 'hierarchical', 
#           y = 'value',
#           facet.by = 'variable') + stat_compare_means()
# 
# ggboxplot(data = toPlot, 
#           x = 'kmeans', 
#           y = 'value',
#           facet.by = 'variable') + stat_compare_means()
# 
# ggboxplot(data = toPlot, 
#           x = 'som', 
#           y = 'value',
#           facet.by = 'variable') + stat_compare_means()
# 
# ggboxplot(data = toPlot, 
#           x = 'mclust', 
#           y = 'value',
#           facet.by = 'variable') + stat_compare_means()
# 
# ggboxplot(data = toPlot, 
#           x = 'nmf', 
#           y = 'value',
#           facet.by = 'variable') + stat_compare_means()


# Multi-omic
files = list.files(path = path, pattern = 'multi_')
files = files[-grep('coca', files)]
compare = list()
# i = 27
for (i in seq_along(files)) {
  res = readRDS(paste0(path, files[i]))
  
  compare[[i]] = data.frame(
    labels = res$labels,
    ids = names(res$labels),
    algorithm = sapply(strsplit(files[i], '_'), '[[', 2),
    features = sapply(strsplit(files[i], '_'), '[[', 3),
    row.names = 1:length(res$labels)
  )
}
compare = data.table::rbindlist(compare)
head(compare)


# Split by experiment
# ==
experiment = unique(compare$features)

exp = list()
# i = 1
for (i in seq_along(experiment)) {
  byexp = compare[which(compare$features == experiment[i]),]
  
  exp[[i]] = reshape2::dcast(dat = byexp,
                             formula = ids ~ algorithm,
                             fun.aggregate = sum,
                             value.var = 'labels')
  
  info = clinical[match(exp[[i]]$id, rownames(clinical)), 17:22]
  exp[[i]] = cbind.data.frame(exp[[i]], info)

}
names(exp) = experiment

toPlot = exp[[1]]
algs = names(toPlot)[2:7]
toPlot = melt(toPlot, 
              measure.vars = c('Monocytes', 'Neutrophils', 'BCells',
                               'NKCells', 'CD4TCells', 'CD8TCells'))


head(toPlot)
ggboxplot(data = toPlot, 
          x = 'hierarchical',
          y = 'value',
          facet.by = 'variable') + stat_compare_means()

ggboxplot(data = toPlot,
          x = 'kmeans',
          y = 'value',
          facet.by = 'variable') + stat_compare_means()

ggboxplot(data = toPlot,
          x = 'som',
          y = 'value',
          facet.by = 'variable') + stat_compare_means()

ggboxplot(data = toPlot,
          x = 'mclust',
          y = 'value',
          facet.by = 'variable') + stat_compare_means()

ggboxplot(data = toPlot,
          x = 'nmf',
          y = 'value',
          facet.by = 'variable') + stat_compare_means()

ggboxplot(data = toPlot,
          x = 'snf',
          y = 'value',
          facet.by = 'variable') + stat_compare_means()

ggboxplot(data = toPlot,
          x = 'icluster',
          y = 'value',
          facet.by = 'variable') + stat_compare_means()

