setwd('~/git/Review-MLAID/benchmark/results/')


files = list.files()
time = list()
res = list()
# f = 1
# i = 1
for (f in seq_along(files)) {
  data = readRDS(files[f])
  # data = data[-c(3, 4, 8)]
  for (i in seq_along(data)) {
    tt = data[[i]]$time 
    if (i == 7) {
      time[[i]] = as.numeric(tt) * 60
    } else{
      time[[i]] = as.numeric(tt)
    }
  }
  
  res[[f]] = data.frame(
    cohort = files[f],
    time = unlist(time),
    samples = length(data[[1]]$labels),
    algorithm = names(data)
  )
  
}

res = data.table::rbindlist(res)
str(res)
res$omic = NA
res$omic[grep('metagenomic', res$cohort)] = 'metagenomic'
res$omic[grep('transcriptomic', res$cohort)] = 'transcriptomic'
res$omic[grep('epigenomic', res$cohort)] = 'epigenomic'

require(ggpubr)

res$time.log10 = log10(res$time)

plot_time = ggplot(res, aes(x=samples, y=time.log10, color = algorithm)) + 
  geom_line() +
  scale_color_manual(values = viridis(5)) +
  # geom_smooth(method = lm, se = F) +
  facet_wrap(~omic, scales = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_time
