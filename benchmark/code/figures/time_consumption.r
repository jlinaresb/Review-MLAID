setwd('~/git/Review-MLAID/benchmark/results/preciseads/')

# Single-omic clustering
# ===
files = list.files(pattern = 'single')
time = list()
for (i in seq_along(files)) {
  res = readRDS(files[i])
  
  time[[i]] = data.frame(
    Consensus.H  = res[[1]]$time, 
    Consensus.KM = res[[2]]$time,
    SOM          = res[[3]]$time,
    MClust       = res[[4]]$time,
    NMF          = res[[5]]$time,
    fs           = files[i]
  )
}

time = data.table::rbindlist(time)
nfeatures = c(43, 43, 279, 279, 479, 101, 101, 3775, 3775, 3775)
time$nfeatures = nfeatures

time = time[-grep('noctrls', time$fs), ]
time = time[-6,]

require(reshape2)
toPlot = melt(time, id.vars = c('nfeatures', 'fs'))
toPlot$value = log10(as.numeric(toPlot$value))

require(ggpubr)
require(viridis)
plot_time = ggplot(toPlot, aes(x=nfeatures, y=value, color = variable)) + 
  geom_line() +
  scale_color_manual(values = viridis(5)) +
  # geom_smooth(method = lm, se = F) +
  # facet_wrap(~omic, scales = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = 'top')
plot_time






# Multi-omic clustering
# ===
files = list.files(pattern = 'multi')

time = list()
for (i in seq_along(files)) {
  res = readRDS(files[i])
  
  # features type
  features = gsub('.rds', '', sapply(strsplit(files[i], '_'), '[[', 3))
  
  # nfeatures
  if (features == 'fcbf-noctrls') {
    nfeatures = 144
  } else if (features == 'fcbf'){
    nfeatures = 144
  } else if (features == 'levene-kruskal-noctrls'){
    nfeatures = 4054
  } else if (features == 'levene-kruskal'){
    nfeatures = 4054
  } else if (features == 'levene'){
    nfeatures = 4254
  }
  
  time[[i]] = data.frame(
    time = res$time,
    algorithm = sapply(strsplit(files[i], '_'), '[[', 2),
    features = features,
    nfeatures = nfeatures
  )
}
time = data.table::rbindlist(time)
time$time = as.numeric(time$time)
time = as.data.frame(time)
time = time[-grep('noctrls', time$features), ]
time$time = log10(time$time)


toPlot = time

plot_time = ggplot(toPlot, aes(x=nfeatures, y=time, color = algorithm)) + 
  geom_line() +
  scale_color_manual(values = viridis(8)) +
  # geom_smooth(method = lm, se = F) +
  # facet_wrap(~omic, scales = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = 'top')
plot_time


# ggsave(plot = plot_time, 
#        filename = 'time_singleomic.pdf',
#        device = 'pdf',
#        path = '~/git/Review-MLAID/benchmark/plots/',
#        height = 5,
#        width = 7)
