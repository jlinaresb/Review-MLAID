setwd('~/git/Review-MLAID/benchmark/results/multi-omic/')
source('~/git/Review-MLAID/benchmark/code/utils.r')

# GSE117931
exp = readRDS('GSE117931_expression.rds')
met = readRDS('GSE117931_methylation.rds')
multi = readRDS('multiomic_GSE117931.rds')

expression = c(exp$Consensus.H$nclust,
exp$Consensus.KM$nclust,
exp$SOM$nclust,
exp$MClust$nclust,
exp$NMF$nclust,
NA, NA)

methylation = c(met$Consensus.H$nclust,
met$Consensus.KM$nclust,
met$SOM$nclust,
met$MClust$nclust,
met$NMF$nclust,
NA, NA)

multiomic = c(multi$Consensus.H$nclust,
multi$Consensus.KM$nclust,
multi$SOM$nclust,
multi$MClust$nclust,
multi$NMF$nclust,
multi$SNF$nclust,
select_k_icluster(multi$iCluster$fit)$nclust)

res = data.frame(
  'Expression' = expression,
  'Methylation' = methylation,
  'Exp_Meth' = multiomic,
  row.names = c('Hierarchical', 'Kmeans', 'SOM', 'MClust', 'NMF', 'SNF', 'iCluster')
)

datatoPlot = as.data.frame(t(res))
require(tidyverse)

datasets = rownames(datatoPlot)
algorithms = colnames(datatoPlot)

data <- expand.grid(datasets = datasets, algorithms = algorithms) %>% bind_cols(results = unlist(datatoPlot))
data = data[which(data$algorithms != 'data'),]
data$results = as.numeric(data$results)
data$cohort = 'GSE117931'
data1 = data

# ## plot data
# require(ggplot2)
# require(viridis)
# plot_nclust = ggplot(data, aes(algorithms, datasets)) +
#   geom_tile(aes(fill = results)) + 
#   geom_text(aes(label = results)) +
#   scale_fill_gradient(low = viridis(1), high = viridis(8)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# plot_nclust
# 
# 
# ggsave(plot = plot_nclust, 
#        filename = 'GSE117931_nclust_multiomic.pdf',
#        device = 'pdf',
#        path = '~/git/Review-MLAID/benchmark/plots/',
#        height = 5,
#        width = 7)



rm(list = setdiff(ls(), c('data1', 'select_k_icluster')))





# GSE82221
exp = readRDS('GSE82221_expression.rds')
met = readRDS('GSE82221_methylation.rds')
multi = readRDS('multiomic_GSE82221.rds')

expression = c(exp$Consensus.H$nclust,
               exp$Consensus.KM$nclust,
               exp$SOM$nclust,
               exp$MClust$nclust,
               exp$NMF$nclust,
               NA, NA)

methylation = c(met$Consensus.H$nclust,
                met$Consensus.KM$nclust,
                met$SOM$nclust,
                met$MClust$nclust,
                met$NMF$nclust,
                NA, NA)

multiomic = c(multi$Consensus.H$nclust,
              multi$Consensus.KM$nclust,
              multi$SOM$nclust,
              multi$MClust$nclust,
              multi$NMF$nclust,
              multi$SNF$nclust,
              select_k_icluster(multi$iCluster$fit)$nclust)

res = data.frame(
  'Expression' = expression,
  'Methylation' = methylation,
  'Exp_Meth' = multiomic,
  row.names = c('Hierarchical', 'Kmeans', 'SOM', 'MClust', 'NMF', 'SNF', 'iCluster')
)

datatoPlot = as.data.frame(t(res))
require(tidyverse)

datasets = rownames(datatoPlot)
algorithms = colnames(datatoPlot)

data <- expand.grid(datasets = datasets, algorithms = algorithms) %>% bind_cols(results = unlist(datatoPlot))
data = data[which(data$algorithms != 'data'),]
data$results = as.numeric(data$results)
data$cohort = 'GSE82221'
data2 = data

rm(list = setdiff(ls(), c('data1', 'data2')))

## plot data
require(ggplot2)
require(viridis)
data = rbind.data.frame(data1, data2)
data$datasets = factor(data$datasets, levels = c('Exp_Meth', 'Methylation', 'Expression'))
plot_nclust = ggplot(data, aes(algorithms, datasets)) +
  geom_tile(aes(fill = results)) +
  geom_text(aes(label = results)) +
  facet_wrap(~cohort) +
  theme(strip.text.x =  element_text(angle = 45)) +
  scale_fill_gradient(low = viridis(10)[1], high = viridis(10)[10]) +
  theme(
    axis.text = element_text(size = 7)) +
  theme_bw()

plot_nclust


ggsave(plot = plot_nclust,
       filename = 'nclust_multiomic.pdf',
       device = 'pdf',
       path = '~/git/Review-MLAID/benchmark/plots/',
       height = 2.5,
       width = 7)

# ggsave(plot = plot_nclust,
#        filename = 'nclust_multiomic.png',
#        device = 'png',
#        path = '~/git/Review-MLAID/benchmark/plots/',
#        height = 5,
#        width = 7)
