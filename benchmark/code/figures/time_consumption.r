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
toPlot.s = melt(time, id.vars = c('nfeatures', 'fs'))
toPlot.s$value = log10(as.numeric(toPlot.s$value))

require(ggpubr)
# require(viridis)
# plot_time.s = ggplot(toPlot, aes(x=nfeatures, y=value, color = variable)) + 
#   geom_line() +
#   scale_color_manual(values = viridis(5)) +
#   # geom_smooth(method = lm, se = F) +
#   # facet_wrap(~omic, scales = 'free_x') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(legend.position = 'top')

toPlot.s$omic = 'single'





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


toPlot.m = time
toPlot.m$omic = 'multi'

# plot_time.m = ggplot(toPlot, aes(x=nfeatures, y=time, color = algorithm)) + 
#   geom_line() +
#   scale_color_manual(values = viridis(8)) +
#   # geom_smooth(method = lm, se = F) +
#   # facet_wrap(~omic, scales = 'free_x') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(legend.position = 'top')


head(toPlot.s)
head(toPlot.m)

names(toPlot.s)[2] = 'features'
names(toPlot.s)[3] = 'algorithm'
names(toPlot.s)[4] = 'time'
toPlot.s = toPlot.s[,match(names(toPlot.m), names(toPlot.s))]

toPlot = rbind.data.frame(toPlot.s, toPlot.m)
toPlot$algorithm = as.character(toPlot$algorithm)

toPlot$algorithm[which(toPlot$algorithm == 'Consensus.H')] = 'hierarchical'
toPlot$algorithm[which(toPlot$algorithm == 'Consensus.KM')] = 'kmeans'
toPlot$algorithm[which(toPlot$algorithm == 'MClust')] = 'mclust'
toPlot$algorithm[which(toPlot$algorithm == 'NMF')] = 'nmf'
toPlot$algorithm[which(toPlot$algorithm == 'SOM')] = 'som'

table(toPlot$algorithm)
str(toPlot)
toPlot$time = abs(toPlot$time)

ggplot(toPlot, aes(x=nfeatures, 
                   y=time, 
                   color = algorithm)) +
  geom_line() +
  scale_color_manual(values = viridis(8)) +
  # geom_smooth(method = lm, se = F) +
  facet_wrap(~omic, scales = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = 'top')




# ggsave(plot = plot_time, 
#        filename = 'time_singleomic.pdf',
#        device = 'pdf',
#        path = '~/git/Review-MLAID/benchmark/plots/',
#        height = 5,
#        width = 7)
