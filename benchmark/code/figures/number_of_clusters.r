setwd('~/git/Review-MLAID/benchmark/results/single-omic/')

files = list.files()
nclust = list()
l = list()
results = list()
# f = 1
# i = 1
for (f in seq_along(files)) {
  data = readRDS(files[f])
  # data = data[-c(3)]
  

  for (i in seq_along(data)) {
    l[[i]] = data[[i]]$nclust
  }
  names(l) = names(data)
  res = as.data.frame(l)
  res$data = files[f]
  
  results[[f]] = res
} 

results = data.table::rbindlist(results)
results = as.data.frame(results)
rownames(results) = results$data
results = results[, -6]


# Plotting
# ===
require(tidyverse)


datasets = rownames(results)
algorithms = colnames(results)

data <- expand.grid(datasets = datasets, algorithms = algorithms) %>% bind_cols(results = unlist(results))
data = data[which(data$algorithms != 'data'),]
data$results = as.numeric(data$results)

## plot data
require(ggplot2)
require(viridis)
plot_nclust = ggplot(data, aes(datasets, algorithms)) +
  geom_tile(aes(fill = results)) + 
  geom_text(aes(label = results)) +
  scale_fill_gradient(low = viridis(1), high = viridis(8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_nclust



# ggarrange(plot_time, plot_nclust, nrow = 2, ncol = 1)
