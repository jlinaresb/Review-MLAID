# Run

cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else{
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/utils.r')
source('benchmark/code/multi-omic/cluster_algorithms_multiomics.r')

exp = readRDS('benchmark/data/preciseads/expression_filtered_levene.rds')
met = readRDS('benchmark/data/preciseads/methylation_filtered_levene.rds')

conCluster.h(data1 = exp,
              data2 = met,
              maxK = 10,
              outPath = 'benchmark/results/preciseads/',
              file = 'hierarchical_preciseads_multiOmic',
              return = F,
              save = T)