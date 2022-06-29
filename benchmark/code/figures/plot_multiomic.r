setwd('~/git/Review-MLAID/benchmark/results/multi-omic/')
source('~/git/Review-MLAID/benchmark/code/utils.r')


# GSE117931
exp = readRDS('GSE117931_expression.rds')
met = readRDS('GSE117931_methylation.rds')
multi = readRDS('multiomic_GSE117931.rds')

expression = c(exp$Consensus.H$nclust,
exp$Consensus.KM$nclust,
# exp$SOM$nclust,
# exp$MClust$nclust,
# exp$NMF$nclust,
NA, NA)

methylation = c(met$Consensus.H$nclust,
met$Consensus.KM$nclust,
# met$SOM$nclust,
# met$MClust$nclust,
# met$NMF$nclust,
NA, NA)

multi = c(multi$Consensus.H$nclust,
multi$Consensus.KM$nclust,
multi$SNF$nclust,
select_k_icluster(multi$iCluster$fit)$nclust)

res = data.frame(
  'Expression' = expression,
  'Methylation' = methylation,
  'Exp_Meth' = multi,
  row.names = c('Hierarchical', 'Kmeans', 'SNF', 'iCluster')
)

t(res)


# GSE82221
exp = readRDS('GSE82221_expression.rds')
met = readRDS('GSE82221_methylation.rds')
multi = readRDS('multiomic_GSE82221.rds')

expression = c(exp$Consensus.H$nclust,
               exp$Consensus.KM$nclust,
               # exp$SOM$nclust,
               # exp$MClust$nclust,
               # exp$NMF$nclust,
               NA, NA)

methylation = c(met$Consensus.H$nclust,
                met$Consensus.KM$nclust,
                # met$SOM$nclust,
                # met$MClust$nclust,
                # met$NMF$nclust,
                NA, NA)

multi = c(multi$Consensus.H$nclust,
          multi$Consensus.KM$nclust,
          multi$SNF$nclust,
          select_k_icluster(multi$iCluster$fit)$nclust)

res = data.frame(
  'Expression' = expression,
  'Methylation' = methylation,
  'Exp_Meth' = multi,
  row.names = c('Hierarchical', 'Kmeans', 'SNF', 'iCluster')
)

t(res)


