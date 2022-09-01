setwd('~/git/Review-MLAID/benchmark/results/preciseads/')

# Single-omic clustering
# ===
files = list.files(pattern = 'single')
nclust1 = list()
for (i in seq_along(files)) {
  res = readRDS(files[i])
  
  nclust1[[i]] = data.frame(
    Consensus.H  = res[[1]]$nclust, 
    Consensus.KM = res[[2]]$nclust,
    SOM          = res[[3]]$nclust,
    MClust       = res[[4]]$nclust,
    NMF          = res[[5]]$nclust,
    omic         = sapply(strsplit(files[i], '_'), '[[', 2),
    features     = sapply(strsplit(files[i], '_'), '[[', 3)
  )
}

nclust1 = data.table::rbindlist(nclust1)
nclust1 = reshape2::melt(nclust1,
                         id.vars = c('omic', 'features'), 
                         value.name = 'nclust')
names(nclust1)[3] = 'algorithm'


# Multi-omic clustering
# ===
files = list.files(pattern = 'multi')
nclust2 = list()
for (i in seq_along(files)) {
  res = readRDS(files[i])
  
  nclust2[[i]] = data.frame(
    nclust = res$nclust,
    algorithm = sapply(strsplit(files[i], '_'), '[[', 2),
    features = sapply(strsplit(files[i], '_'), '[[', 3),
    omic = 'multi'
  )
}
nclust2 = data.table::rbindlist(nclust2)

# Create data to plot
nclust1 = nclust1[, match(names(nclust2), names(nclust1))]
toPlot = rbind.data.frame(nclust1, nclust2)

toPlot$algorithm[which(toPlot$algorithm == 'Consensus.H')] = 'hierarchical'
toPlot$algorithm[which(toPlot$algorithm == 'Consensus.KM')] = 'kmeans'
toPlot$algorithm[which(toPlot$algorithm == 'MClust')] = 'mclust'
toPlot$algorithm[which(toPlot$algorithm == 'NMF')] = 'nmf'
toPlot$algorithm[which(toPlot$algorithm == 'SOM')] = 'som'

# Plotting
# ===
require(ggpubr)
require(viridis)
plot = ggplot(toPlot, aes(x=algorithm, y=nclust, color=algorithm)) + 
  geom_violin(trim = T) +
  geom_jitter(position=position_jitter(.5)) +
  scale_color_manual(values = viridis(8)) +
  facet_wrap(~omic, scales = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = 'top') +
  scale_y_continuous(breaks=seq(0, 11, 1), limits=c(0, 11))
plot

# ggsave(plot = plot_time, 
#        filename = 'time_singleomic.pdf',
#        device = 'pdf',
#        path = '~/git/Review-MLAID/benchmark/plots/',
#        height = 5,
#        width = 7)
